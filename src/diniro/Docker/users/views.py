from django.shortcuts import render
import random
import string
from django.contrib.auth.models import User
from django.contrib.auth import login, authenticate,logout
from django.conf import settings
from django.core.files.storage import FileSystemStorage

def creat_user(request):
    username = ''.join(random.SystemRandom().choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for _ in range(32))
    raw_password = ''.join(random.SystemRandom().choice(string.ascii_lowercase + string.ascii_uppercase + string.digits) for _ in range(32))
    user = User(username=username, first_name=raw_password)
    user.set_unusable_password()
    user.save()
    # u = authenticate(user=u)
    login(request, user)
    context = {'data': username}
    #logout(request)
    return render(request, 'users/username.html', context)



