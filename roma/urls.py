from django.urls import path, re_path

from . import views

from django.contrib import admin
from django.urls import path
from roma.views import * 

urlpatterns = [
    path('prony/', pronyView, name='roma'),
    #re_path(r'clusters/(?P<d>[\w_]+)/$',views.clusters),
    #path('mapmx/', mapmx),
    #path('frequency/', frequency, name='frequency'),
    #path('angle/', angle, name='angle'),
    #path('Home/', Home, name='Home'), 
    #path('Zonemap/', Zonemap, name='Zonemap'),
    #path('hola/', signals),
]

