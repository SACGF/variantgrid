import zipfile
from abc import ABC, abstractmethod
from datetime import datetime
from io import StringIO
from typing import Optional, List, Iterator, Tuple, Any

from django.conf import settings
from django.http.response import HttpResponseBase
from django.http import HttpResponse, StreamingHttpResponse
from more_itertools import peekable
from stream_zip import NO_COMPRESSION_64, stream_zip, NO_COMPRESSION_32, ZIP_64, ZIP_32
from threadlocals.threadlocals import get_current_request
from library.guardian_utils import bot_group
from library.log_utils import NotificationBuilder
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter


class ClassificationExportFormatter2(ABC):
    """
    Extend this class to export classification data into different formats
    """

    def __init__(self, classification_filter: ClassificationFilter):
        self.classification_filter = classification_filter
        self.row_count = 0
        self.file_count = 0
        self.started = datetime.utcnow()

    @property
    def is_genome_build_relevant(self) -> bool:
        return True

    def filename(self, part: Optional[int] = None, extension_override: Optional[str] = None) -> str:
        """
        Generate a filename
        :param part: If the data is being split, what index is this file (otherwise None)
        :param extension_override: If creating a wrapper file, e.g. "zip"
        :return: The appropriate filename
        """
        filename_parts: List[str] = [self.classification_filter.file_prefix]
        if self.classification_filter.file_include_date:
            filename_parts.append(self.classification_filter.date_str)

        if self.is_genome_build_relevant:
            filename_parts.append(str(self.classification_filter.genome_build))

        if part is not None:
            filename_parts.append(f"part_{part:02}")

        filename = "_".join(filename_parts)

        return f"{filename}.{extension_override or self.extension()}"

    def serve(self) -> HttpResponseBase:
        """
        Start generating the data and return it in a HTTP Response
        """
        if self.classification_filter.rows_per_file:
            response = StreamingHttpResponse(stream_zip(self._yield_streaming_files()), content_type='application/zip')
            response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
            return response
        else:
            self.file_count = 1
            # can stream in single file
            response = StreamingHttpResponse(self._yield_file(), self.content_type())
            response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
            return response

    @abstractmethod
    def content_type(self) -> str:
        """
        :return: Http content type
        """
        pass

    def _yield_streaming_entry(self, source: peekable) -> Iterator[bytes]:
        # source should be a peekable of List[str]

        def byte_me(rows: List[str]) -> bytes:
            return "\n".join(rows).encode()

        if header := self.header():
            yield byte_me(header)

        this_entry_row_count = 0
        while True:
            # peek the next entry to see if it's big enough that we should split the file
            # if it's not too big (or if we've reached the end where False will be returned)
            # then add the footer and this iteration is done
            if next_rows := source.peek(False):
                if this_entry_row_count + len(next_rows) > self.classification_filter.rows_per_file:
                    break
                # now that we're sure we want the rows, call next (which will give us the same data)
                # but progress the iterator
                next_rows = next(source)
                next_rows_count = len(next_rows)
                this_entry_row_count += next_rows_count
                self.row_count += next_rows_count
                yield byte_me(next_rows)
            else:
                break

        if footer := self.footer():
            yield byte_me(footer)

    def _yield_streaming_files(self) -> Iterator[Tuple[str, datetime, int, Any, bytes]]:
        def row_iterator():
            for allele_data in self.classification_filter.allele_data_filtered():
                if rows := self.row(allele_data):
                    yield rows

        modified_at = datetime.now()
        perms = 0o600

        data_peek = peekable(row_iterator())
        while data_peek.peek(False):
            self.file_count += 1
            yield self.filename(part=self.file_count), modified_at, perms, NO_COMPRESSION_64, self._yield_streaming_entry(data_peek)


    def _yield_file(self) -> Iterator[str]:
        """
        :return: An iterator for a single streaming file, call either this or yield_file
        """
        try:
            for header in self.header():
                yield header
            for allele_data in self.classification_filter.allele_data_filtered():

                row_data = self.row(allele_data)
                for row in row_data:
                    self.row_count += 1
                    yield row
            for footer in self.footer():
                yield footer
        except:
            yield "An error occurred generating the file"
            raise

        self.send_stats()

    @abstractmethod
    def extension(self) -> str:
        """
        Filename extension, e.g. csv, json
        """
        pass

    def header(self) -> List[str]:
        """
        Return rows to start each file, typically 0 to 1 line
        :return: A list of rows to be \n at the top of each file
        """
        return list()

    @abstractmethod
    def row(self, allele_data: AlleleData) -> List[str]:
        """
        Return the row or rows that represent this AlleleData
        :param allele_data: All the valid classifications for an Allele
        """
        pass

    def footer(self) -> List[str]:
        """
        Return rows to end each file, typically 0
        :return: A list of rows to be \n at the bottom of each file
        """
        return list()

    def send_stats(self):
        """
        Alerts Slack of download
        """

        # don't report bots downloading
        user = self.classification_filter.user
        if user.groups.filter(name=bot_group().name):
            return
        end = datetime.utcnow()
        row_count = self.row_count

        body_parts = [f":simple_smile: {user.username}"]
        if request := get_current_request():
            body_parts.append(f"URL : `{request.path_info}`")
        body_parts.append(f"Rows Downloaded : *{row_count}*")
        if self.file_count > 1:
            body_parts.append(f"File Count : *{self.file_count}*")

        nb = NotificationBuilder(message="Classification Download")\
            .add_header(":arrow_down: Classification Download Completed")\
            .add_markdown("\n".join(body_parts), indented=True)
        if request := get_current_request():
            for key, value in request.GET.items():
                nb.add_field(key, value)
        nb.add_field("Duration", str((end - self.started).seconds) + " seconds")
        nb.send()
