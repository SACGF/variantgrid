import logging
import os
import shutil
import subprocess
import sys
from tempfile import NamedTemporaryFile

from django.conf import settings
from django.utils.crypto import salted_hmac

from pedigree.models import PedFileFamily


def get_ped_parser_command() -> list[str]:
    """ ped_parser installs a plain script (not a console_scripts entry point or a runnable
        module), so we invoke it with the current interpreter. This avoids relying on the venv's
        bin dir being on PATH, which fails for local installs (e.g. UV) when Django runs as a
        service - see https://github.com/SACGF/variantgrid/issues/1572 """
    ped_parser = os.path.join(os.path.dirname(sys.executable), "ped_parser")
    if not os.path.exists(ped_parser):
        ped_parser = shutil.which("ped_parser") or ped_parser
    return [sys.executable, ped_parser]


class PedigreeChart:
    """ One chart per family: Madeline2 names its output <prefix>_<family>.svg when the input holds
        several families, so only a single family gives the <prefix>.svg we wait for.

        To use this you need to install:
            * ped_parser - python3 -m pip install ped_parser
            * Madeline2 -  https://madeline.med.umich.edu/madeline/
    """
    def __init__(self, ped_file_family_id):
        self.ped_file_family_id = ped_file_family_id

    def get_params_hash(self):
        """ This should be unique for every graph produced (from parameters). Salted so the media URL of a
            family's chart can't be worked out from its id """
        return salted_hmac("pedigree_chart", f"ped_file_family={self.ped_file_family_id}",
                           algorithm="sha256").hexdigest()

    def get_filename(self):
        return self._get_prefix() + self._get_extension()

    def _get_prefix(self):
        return os.path.join(settings.GENERATED_DIR, self.get_name(), self.get_params_hash())

    @staticmethod
    def _get_extension():
        return ".svg"

    def get_name(self):
        return self.__class__.__name__

    def save(self, _filename):
        madeline2 = getattr(settings, "PEDIGREE_MADELINE2_COMMAND", None)
        if madeline2 is None:
            raise ValueError("settings.PEDIGREE_MADELINE2_COMMAND not set!")
        ped_file_family = PedFileFamily.objects.get(pk=self.ped_file_family_id)
        with NamedTemporaryFile() as madeline_file, NamedTemporaryFile(suffix=".ped") as ped_temp_file:
            madeline_filename = madeline_file.name
            ped_file_family.write_ped(ped_temp_file.name)

            # Need to convert to Madeline input file
            convert_command = [
                *get_ped_parser_command(),
                ped_temp_file.name,
                "--to_madeline",
                "-o",
                madeline_filename,
            ]
            subprocess.check_call(convert_command)
            logging.info("Converted .ped to madeline: %s", madeline_filename)

            command = [
                madeline2,
                madeline_filename,
                "--outputprefix",
                self._get_prefix(),
            ]
            logging.info("Running: \n'%s'", " ".join(command))
            output = subprocess.check_output(command)
        output = output.decode()
        success = "Pedigree output file is" in output
        if not success:
            raise OSError(output)
