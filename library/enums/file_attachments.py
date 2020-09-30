import os


class AttachmentFileType:
    """ JFU file attached to some other object """
    IMAGE = "I"
    PDF = "P"
    WORD_DOC = "W"
    TEXT = "T"
    SPREADSHEET = "S"
    OTHER = "O"
    CHOICES = (
        (IMAGE, "Image"),
        (PDF, "Image"),
        (WORD_DOC, "Word Doc"),
        (TEXT, "Text"),
        (SPREADSHEET, "Spreadsheet"),
        (OTHER, "Image"),
    )

    _EXTENSIONS = {IMAGE: ['.png', '.gif', '.jpg', '.jpeg'],
                   PDF: ['.pdf'],
                   WORD_DOC: [".doc", ".docx"],
                   TEXT: ['.txt'],
                   SPREADSHEET: ['.csv', '.xls', '.xlsx']}

    @staticmethod
    def get_type_for_file(filename):
        types_by_extension = {}
        for k, v_list in AttachmentFileType._EXTENSIONS.items():
            for ext in v_list:
                types_by_extension[ext] = k

        (_, ext) = os.path.splitext(filename)
        ext = ext.lower()
        return types_by_extension.get(ext, AttachmentFileType.OTHER)
