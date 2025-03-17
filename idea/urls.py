from django.contrib import admin
from django.urls import path
from idea.views import idea, csv_write, user_guide, latex_code

from django.conf import settings
from django.conf.urls.static import static

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', idea, name='home'),  # Agregar esta línea para que la URL raíz redirija a la vista `idea`
    #path('idea/', idea, name='idea'),
    path('user_guide/', user_guide, name='user_guide'),
    path('csv-write/', csv_write, name='csv_write'),
    path('latex_code/', latex_code, name='latex_code'),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)