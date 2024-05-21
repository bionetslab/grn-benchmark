from django.contrib.auth.models import User
from django.contrib import messages
from django.shortcuts import render
from users.models import UserFiles


class Common():

    def del_user(request, username):
        print(username)
        try:
            u = User.objects.get(username = username)
            user_files = UserFiles.objects.get(user=request.user.username)
            u.delete()
            user_files.delete()
            messages.success(request, "The user is deleted")

        except User.DoesNotExist:
            messages.error(request, "User doesnot exist")
            return render(request, 'test_interface/HOME.html')

        except Exception as e:
            return render(request, 'test_interface/HOME.html',{'err':e.message})

        return render(request, 'test_interface/HOME.html')