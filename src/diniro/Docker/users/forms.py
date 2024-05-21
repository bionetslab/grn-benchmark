from django import forms
from users.models import UserFiles

class UserFilesForm(forms.ModelForm):
    widget = forms.FileInput

    class Meta:
        model = UserFiles
        fields = ('main_anndata','main_TFs','sessionname',)

    def __init__(self, *args, **kwargs):
        super(UserFilesForm, self).__init__(*args, **kwargs)

        # sets the placeholder key/value in the attrs for a widget
        # when the form is instantiated (so the widget already exists)
        self.fields['sessionname'].widget.attrs['placeholder'] = 'Session Name ...'



