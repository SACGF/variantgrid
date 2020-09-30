import logging
import os
import subprocess
from tempfile import NamedTemporaryFile

from django.conf import settings

from library.utils import sha1_str
from pedigree.models import PedFile


class PedigreeChart:
    """ To use this you need to install:
            * ped_parser - python3 -m pip install ped_parser
            * Madeline2 -  https://madeline.med.umich.edu/madeline/
    """
    def __init__(self, ped_file_id):
        self.ped_file_id = ped_file_id

    def get_params_hash(self):
        """ This should be unique for every graph produced (from parameters) """
        return sha1_str(str(self.ped_file_id))

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
        ped_file = PedFile.objects.get(pk=self.ped_file_id)
        ped_filename = ped_file.uploadedpedfile.uploaded_file.get_filename()
        temp_file = NamedTemporaryFile(delete=False)
        madeline_filename = temp_file.name

        # Need to convert to Madeline input file
        convert_command = [
            "ped_parser",
            ped_filename,
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
        os.unlink(madeline_filename)
        output = output.decode()
        success = "Pedigree output file is" in output
        if not success:
            raise OSError(output)
