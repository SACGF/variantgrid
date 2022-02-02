import zipfile
from abc import ABC, abstractmethod
from datetime import datetime
from io import StringIO
from typing import Optional, List, Iterator

from django.conf import settings
from django.http.response import HttpResponseBase
from django.http import HttpResponse, StreamingHttpResponse
from threadlocals.threadlocals import get_current_request
from library.guardian_utils import bot_group
from library.log_utils import NotificationBuilder
from classification.views.exports.classification_export_filter import AlleleData, ClassificationFilter


class FileWriter:
    """
    In memory file writer, used if export create a zip file with the data potentially split across
    multiple entries
    :var file: StringIO data to write to
    :var row_count: How many lines have been written to this file
    """

    def __init__(self):
        self.file = StringIO()
        self.row_count = 0

    def __bool__(self):
        return self.row_count > 0

    def write(self, rows: List[str], count: bool = True):
        self.file.writelines(rows)
        if count:
            self.row_count += len(rows)

    def finish(self):
        content = self.file.getvalue()
        self.file.close()
        return content


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
            filename_parts.append(f"part_{part+1:02}")

        filename = "_".join(filename_parts)

        return f"{filename}.{extension_override or self.extension()}"

    def serve(self) -> HttpResponseBase:
        """
        Start generating the data and return it in a HTTP Response
        """
        if self.classification_filter.rows_per_file:
            response = HttpResponse(content_type='application/zip')
            response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'

            with zipfile.ZipFile(response, 'w') as zf:
                # zf.writestr("readme.txt", f"This file was downloaded from {settings.SITE_NAME}\nThis is buffer" + ("0") * 16384)
                # response.flush()

                for index, entry in enumerate(self._yield_files()):
                    self.file_count += 1
                    # can we be more efficient than converting StringIO to a string to be converted back into bytes?
                    zf.writestr(self.filename(part=index), str(entry.file.getvalue()))

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

    def _yield_files(self) -> Iterator[FileWriter]:
        """
        :return: An iterator of chunked file data (each value could be made into a ZipEntry)
        """
        fw: Optional[FileWriter] = None

        for allele_data in self.classification_filter.allele_data_filtered():
            to_rows = self.row(allele_data)
            self.row_count += len(to_rows)
            if not fw or (self.classification_filter.rows_per_file and fw.row_count + len(to_rows) > self.classification_filter.rows_per_file):
                if fw:
                    fw.write(self.footer(), count=False)
                    yield fw
                fw = FileWriter()
                fw.write(self.header(), count=False)

            fw.write(to_rows, count=True)

        if fw:
            fw.write(self.footer(), count=False)
            yield fw
        self.send_stats()

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
