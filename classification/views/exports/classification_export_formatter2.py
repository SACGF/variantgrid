import zipfile
from abc import ABC, abstractmethod
from datetime import datetime
from io import StringIO
from typing import Optional, List, Iterator, Tuple, Any

from django.http import HttpResponse, StreamingHttpResponse
from django.http.response import HttpResponseBase
from django.shortcuts import render
from more_itertools import peekable
from stream_zip import stream_zip, ZIP_64
from threadlocals.threadlocals import get_current_request

from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter
from library.guardian_utils import bot_group
from library.log_utils import NotificationBuilder, report_exc_info


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

        # set by the decorator
        filename_parts.append(self.format_type)

        if part is not None:
            filename_parts.append(f"part_{part:02}")

        filename = "_".join(filename_parts)

        return f"{filename}.{extension_override or self.extension()}"

    def serve(self) -> HttpResponseBase:
        """
        Start generating the data and return it in a HTTP Response
        """
        if self.classification_filter.benchmarking:
            for count, row in enumerate(self._yield_single_file()):
                if count > 1000:
                    break
            return render(get_current_request(), "snpdb/benchmark.html", {"content": "TODO"})

        if self.classification_filter.rows_per_file:
            # Had subtle issues with stream_zip, maybe try again after a version increase
            # response = StreamingHttpResponse(stream_zip(self._yield_streaming_zip_entries()), content_type='application/zip')
            # response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
            # return response
            return self._non_streaming_zip()
        else:
            self.file_count = 1
            # can stream in single file
            response = StreamingHttpResponse(self._yield_single_file(), self.content_type())
            response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
            return response

    @abstractmethod
    def content_type(self) -> str:
        """
        :return: Http content type
        """
        pass

    def _peekable_data(self) -> peekable:  # peekable[List[str]]
        def row_iterator():
            for allele_data in self.classification_filter.allele_data_filtered_pre_processed():
                if rows := self.row(allele_data):
                    yield rows
        return peekable(row_iterator())

    def _yield_streaming_entry_str(self, source: peekable) -> Iterator[List[str]]:
        # source should be a peekable of List[str]
        # yield's a file's worth of data (in several chunks)

        if header := self.header():
            yield header

        this_entry_row_count = 0
        while True:
            # peek the next entry to see if it's big enough that we should split the file
            # if it's not too big (or if we've reached the end where False will be returned)
            # then add the footer and this iteration is done
            if next_rows := source.peek(default=False):
                next_rows_count = len(next_rows)
                # code somewhat assumes header count = 1
                if not this_entry_row_count == 0 and this_entry_row_count + next_rows_count >= self.classification_filter.rows_per_file:
                    break
                # now that we're sure we want the rows, call next (which will give us the same data)
                # and progress the iterator
                next(source)

                this_entry_row_count += next_rows_count
                self.row_count += next_rows_count
                yield next_rows
            else:
                break

        if footer := self.footer():
            yield footer

    def _yield_streaming_entry_bytes(self, source: peekable) -> Iterator[bytes]:
        for entry in self._yield_streaming_entry_str(source):
            yield "\n".join(entry).encode()

    def _yield_streaming_zip_entries(self) -> Iterator[Tuple[str, datetime, int, Any, bytes]]:
        modified_at = datetime.now()
        perms = 0o600
        try:
            data_peek = self._peekable_data()
            while data_peek.peek(default=False):
                self.file_count += 1
                yield self.filename(part=self.file_count), modified_at, perms, ZIP_64, self._yield_streaming_entry_bytes(data_peek)
        except:
            report_exc_info()
            def yield_error_bytes():
                yield "An error occurred generating the file".encode()
            yield "error.txt", modified_at, perms, ZIP_64, yield_error_bytes()
            raise

    def _streaming_zip(self) -> StreamingHttpResponse:
        # Had some issues with stream_zip telling macOS couldn't extract file, but then being able to manually extract the downloaded file fine
        # not sure if the error is on macOS, stream_zip or my implementation, so using non streaming version for now
        response = StreamingHttpResponse(stream_zip(self._yield_streaming_zip_entries()), content_type='application/zip')
        response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
        return response

    def _non_streaming_zip(self) -> HttpResponse:
        response = HttpResponse(content_type='application/zip')
        with zipfile.ZipFile(response, 'w') as zf:
            data_peek = self._peekable_data()

            data_peek = self._peekable_data()
            while data_peek.peek(default=False):
                self.file_count += 1

                def next_file() -> str:
                    str_buffer = StringIO()
                    for str_data in self._yield_streaming_entry_str(source=data_peek):
                        str_buffer.writelines(str_data)
                    str_buffer.flush()
                    return str_buffer.getvalue()

                zf.writestr(self.filename(part=self.file_count), next_file())

        response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
        self.send_stats()
        return response

    def _yield_single_file(self) -> Iterator[str]:
        """
        :return: An iterator for a single streaming file, call either this or yield_file
        """
        try:
            for header in self.header():
                yield header
            for allele_data in self.classification_filter.allele_data_filtered_pre_processed():

                row_data = self.row(allele_data)
                for row in row_data:
                    self.row_count += 1
                    yield row
            for footer in self.footer():
                yield footer
        except:
            report_exc_info()
            yield "An error occurred generating the file"
            raise
        finally:
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
