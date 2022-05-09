from django.urls import path, re_path

from . import views

from django.contrib import admin
from django.urls import path
from roma.views import * 
from RTangle.views import *

urlpatterns = [    
    path('RTangle/', RTangle, name='RTanlge'),
    ]