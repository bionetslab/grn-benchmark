from django.urls import path
from . import views
from users import views as user_viwes
from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('', views.home, name='test-home'),
    # path('home2/', views.home2, name='test-home2'),
    path('upload/', views.upload, name='test-upload'),
    path('about/', views.about, name='test-about'),
    path('nowhere/', views.selection1, name='test-select'),
    path('nowhere2/', views.lasso_selecton, name='test-lasso'),
    path('result/', views.result, name='test-result'),
    path('result2/', views.result_page, name='test-result-page'),
    path('nowhere3/', views.lasso_reset, name='test-lasso-reset'),
    path('nowhere4/', views.selection2, name='test-sample-selection'),
    path('download/', views.download, name='test-download'),

]
