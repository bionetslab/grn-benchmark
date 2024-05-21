from django.urls import reverse
from django import template

register = template.Library()


@register.simple_tag(name='rurl')
def rurl(url, *args, **kwargs):
    return reverse(url)


@register.simple_tag(name='rstatic')
def rstatic(file, *args, **kwargs):
    return "/static/" + file