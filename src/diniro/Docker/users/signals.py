# from django.db.models.signals import post_save
# from django.contrib.auth.models import User
# from django.dispatch import receiver
# from users.models import UserFiles
#
#
# @receiver(post_save, sender=User)
# def creat_user(sender, instance, created, **kwargs):
#     if created:
#         UserFiles.objects.create(user=instance)
# #
# @receiver(post_save, sender=User)
# def save_user(sender, instance, **kwargs):
#     instance.userfiles.save()