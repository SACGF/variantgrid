import zipfile
from datetime import datetime
from io import StringIO
from typing import Optional, Iterator, Any

from django.http import HttpResponseBase, StreamingHttpResponse, HttpResponse
from django.shortcuts import render
from more_itertools.more import peekable
from stream_zip import stream_zip, ZIP_64
from threadlocals.threadlocals import get_current_request
from classification.views.exports_grouping.classification_grouping_export_filter import \
    ClassificationGroupingExportFormat, ClassificationGroupingExportFileSettings
from library.log_utils import report_exc_info


class ClassificationGroupingExportProcess:

    def __init__(self, classification_export_format: ClassificationGroupingExportFormat, export_settings: ClassificationGroupingExportFileSettings):
        self.classification_export_format = classification_export_format
        self.export_settings = export_settings
        self.format_properties = classification_export_format.format_properties()

        self.row_count = 0
        self.file_count = 0

    def filename(self, part: Optional[int] = None, extension_override: Optional[str] = None) -> str:
        """
        Generate a filename
        :param part: If the data is being split, what index is this file (otherwise None)
        :param extension_override: If creating a wrapper file, e.g. "zip"
        :return: The appropriate filename
        """
        filename_parts: list[str] = []
        if file_prefix := self.format_properties.filename_suffix:
            filename_parts.append(file_prefix)

        # if self.classification_filter.allele_origin_filter != AlleleOriginFilterDefault.SHOW_ALL:
        #     filename_parts.append(self.classification_filter.allele_origin_filter.label.lower())

        # FIXME add date_str back
        # filename_parts.append(self.classification_filter.date_str)

        if self.format_properties.is_genome_build_relevant:
            filename_parts.append(str(self.classification_export_format.classification_grouping_filter.genome_build))

        # set by the decorator
        # noinspection PyUnresolvedReferences
        # filename_parts.append(self.format_type)

        if part is not None:
            if isinstance(part, int):
                filename_parts.append(f"part_{part:02}")
            else:
                filename_parts.append(f"part_{part}")

        filename = "_".join(filename_parts)

        return f"{filename}.{extension_override or self.format_properties.extension}"

    def serve(self) -> HttpResponseBase:
        """
        Start generating the data and return it in an HTTP Response
        """
        benchmarking = False
        if benchmarking:
            for count, row in enumerate(self._yield_single_file()):
                if count > 100:
                    break
            return render(get_current_request(), "snpdb/benchmark.html", {"content": "TODO"})

        if self.export_settings.rows_per_file:
            # Had subtle issues with stream_zip, maybe try again after a version increase
            # response = StreamingHttpResponse(stream_zip(self._yield_streaming_zip_entries()), content_type='application/zip')
            # response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
            # return response
            return self._non_streaming_zip()
        else:
            self.file_count = 1
            # can stream in single file
            response = StreamingHttpResponse(streaming_content=self._yield_single_file(), content_type=self.format_properties.http_content_type)
            response.minify_response = False
            # FIXME
            # response['Last-Modified'] = self.classification_filter.last_modified_header
            response['Content-Disposition'] = f'attachment; filename="{self.filename()}"'
            return response

    def serve_in_memory(self) -> Iterator[str]:
        data_peek = self._peekable_data()

        while data_peek.peek(default=False):
            def next_file() -> str:
                str_buffer = StringIO()
                for str_data in self._yield_streaming_entry_str(source=data_peek):
                    str_buffer.writelines(str_data)
                str_buffer.flush()
                return str_buffer.getvalue()

            yield next_file()

    def _peekable_data(self) -> peekable:  # peekable[list[str]]
        return self.classification_export_format.peekable()

    def with_new_lines(self, data: list[str]) -> list[str]:
        if data:
            return [row + self.format_properties.delimiter_for_row for row in data]
        else:
            return []

    def _yield_streaming_entry_str(self, source: peekable) -> Iterator[list[str]]:
        # source should be a peekable of list[str]
        # yield's a file's worth of data (in several chunks)

        if header := self.with_new_lines(self.classification_export_format.header()):
            yield header

        this_entry_row_count = 0

        is_first_row = True
        while True:
            # peek the next entry to see if it's big enough that we should split the file
            # if it's not too big (or if we've reached the end where False will be returned)
            # then add the footer and this iteration is done
            next_rows = source.peek(default=None)
            if next_rows is not None:
                use_rows = []
                for row in next_rows:
                    if not is_first_row:
                        row = f"{self.format_properties.delimiter_for_row}{row}"
                    else:
                        is_first_row = False
                    use_rows.append(row)

                next_rows_count = len(use_rows)
                # code somewhat assumes header count = 1

                # if we have a rows per file limit, and we're not on row 0 (don't want to get stuck in a scenario where
                # we have infinite files with 0 rows due to some weird configuration)
                if self.export_settings.rows_per_file and not this_entry_row_count == 0:
                    # make sure our existing data, plus the next round of data wont exceed the row limit
                    if this_entry_row_count + next_rows_count >= self.export_settings.rows_per_file:
                        break
                # now that we're sure we want the rows, call next (which will give us the same data)
                # and progress the iterator
                next(source)

                this_entry_row_count += next_rows_count
                self.row_count += next_rows_count
                yield use_rows
            else:
                break

        if footer := self.with_new_lines(self.classification_export_format.footer()):
            yield footer

    def _yield_streaming_entry_bytes(self, source: peekable) -> Iterator[bytes]:
        for entry in self._yield_streaming_entry_str(source):
            yield "\n".join(entry).encode()

    def _yield_streaming_zip_entries(self) -> Iterator[tuple[str, datetime, int, Any, bytes]]:
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
        # FIXME re-establish last_modified_header
        # response['Last-Modified'] = self.classification_filter.last_modified_header
        response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
        return response

    def _non_streaming_zip(self) -> HttpResponse:
        response = HttpResponse(content_type='application/zip')
        with zipfile.ZipFile(response, 'w') as zf:

            data_peek = self._peekable_data()
            while data_peek.peek(default=None) is not None:
                self.file_count += 1

                def next_file() -> str:
                    str_buffer = StringIO()
                    for str_data in self._yield_streaming_entry_str(source=data_peek):
                        str_buffer.writelines(str_data)
                    str_buffer.flush()
                    return str_buffer.getvalue()

                zf.writestr(self.filename(part=self.file_count), next_file())

            # if extra := self.extra_data(as_individual_file=True):
            #     str_buffer = StringIO()
            #     str_buffer.write(extra.content)
            #     zf.writestr(self.filename(part=extra.filename_part), str_buffer.getvalue())

        # FIXME response['Last-Modified'] = self.classification_filter.last_modified_header
        response['Content-Disposition'] = f'attachment; filename="{self.filename(extension_override="zip")}"'
        self.send_stats()
        return response

    def _yield_single_file(self) -> Iterator[str]:
        """
        :return: An iterator for a single streaming file, call either this or yield_file
        """
        try:
            for header in self.with_new_lines(self.classification_export_format.header()):
                yield header

            for rows in self.classification_export_format.row_generator():
                rows = self.with_new_lines(rows)
                for row in rows:
                    self.row_count += 1

                    # if not is_first_row:
                    #     row = f"{self.format_properties.delimiter_for_row}{row}"
                    # else:
                    #     is_first_row = False

                    yield row

            for footer in self.with_new_lines(self.classification_export_format.footer()):
                yield footer

            # if extra_data := self.extra_data(as_individual_file=False):
            #     yield extra_data.content

        except:
            report_exc_info()
            yield "An error occurred generating the file"
            raise
        finally:
            self.send_stats()

    def send_stats(self):
        pass