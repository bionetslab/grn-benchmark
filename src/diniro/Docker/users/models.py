from django.db import models
from django.contrib.auth.models import User
from django.conf import settings
from django.core.validators import FileExtensionValidator
from django_mysql.models import ListCharField
from django import forms

# Create your models here.


def upload_location(instance, filename):
    return "%s/%s" %(instance.user,filename)

def my_default():
    return {'foo': 'bar'}

class UserFiles(models.Model):

    user = models.CharField(max_length=1000,null=True)
    main_anndata = models.FileField(upload_to=upload_location, default='anndata.h5ad', validators=[FileExtensionValidator(allowed_extensions=['h5ad'])])
    main_TFs = models.FileField(upload_to=upload_location, default='', blank=True, null=True, validators=[FileExtensionValidator(allowed_extensions=['txt'])])
    points = models.JSONField(default=my_default)
    colors = models.JSONField(default=my_default)
    samples = models.JSONField(default=my_default)
    sessionname = models.CharField(max_length=50, null=True)
    map = models.CharField(max_length=1000, null=True)
    obs = models.CharField(max_length=1000, null=True)
    sample1x_lasso = models.JSONField(default=my_default)
    sample2x_lasso = models.JSONField(default=my_default)
    sample1y_lasso = models.JSONField(default=my_default)
    sample2y_lasso = models.JSONField(default=my_default)
    selectionType = models.JSONField(default={'selectionType': 'None'})
    sample1Name = models.CharField(max_length=10, null=True)
    sample2Name = models.CharField(max_length=10, null=True)





    def __str__(self):
        return f'{self.user} userfile'