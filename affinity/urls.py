from django.urls import path, re_path

from . import views

urlpatterns = [
    path('clusters/', views.clusters, name='clusters'),
    re_path(r'clusters/(?P<d>[\w_]+)/$',views.clusters)
]
