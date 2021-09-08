import os

from django.db import models
from django.db.models.deletion import CASCADE
from django.urls.base import reverse

from classification.models.classification import Classification
from library.django_utils.django_file_system_storage import PrivateUploadStorage
from library.enums.file_attachments import AttachmentFileType


class ClassificationAttachment(models.Model):
    UPLOAD_PATH = 'classification_attachments'  # in media root

    classification = models.ForeignKey(Classification, on_delete=CASCADE)
    file = models.FileField(upload_to=UPLOAD_PATH,
                            storage=PrivateUploadStorage())
    file_type = models.CharField(max_length=1, choices=AttachmentFileType.CHOICES)
    thumbnail_path = models.TextField(null=True)

    def get_file_dict(self):
        basename = os.path.basename(self.file.path)
        image_url = self.get_absolute_url()

        if self.file_type == AttachmentFileType.IMAGE:
            thumb_url = reverse('view_classification_file_attachment_thumbnail', kwargs={'pk': self.pk})
        else:
            thumb_url = None

        if os.path.exists(self.file.path):
            size = self.file.size
        else:
            size = 0

        return {
            'name': basename,
            'size': size,
            'url': image_url,
            'thumbnailUrl': thumb_url,  # If thumbnail set, JFU displays in gallery
            'deleteUrl': reverse('classification_file_delete', kwargs={'pk': self.pk}),
            'deleteType': 'POST',
        }

    def get_absolute_url(self):
        return reverse('view_classification_file_attachment', kwargs={'pk': self.pk})
