

import csv
import json
from datetime import datetime

from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.http import HttpResponse
from django.shortcuts import render

import idea.plotly_graphics as ply_graph
import idea.plt_graphics as plt_graph
from idea.analyze_data import analyze, get_data, get_signals_name

graficos_generados = []
def idea(request):
    global graficos_generados
    dic = {}

    if request.method == 'POST' and 'upload' in request.POST:
        uploaded_file = request.FILES['examinar']
        fs = FileSystemStorage()
        fs.save(uploaded_file.name, uploaded_file)
        dic['file_name'] = uploaded_file.name

    elif request.method == 'POST' and 'analyze' in request.POST:
        dic['file_name'] = request.POST.get('file_name', '')
        dic['method'] = request.POST.get('method', '')
        dic['fs'] = request.POST.get('fs', '')
        dic['lfft'] = request.POST.get('lfft', '')
        dic['tolerance'] = request.POST.get('tolerance', '')
        dic['windows'] = request.POST.get('windows', '')
        dic['hwindows'] = request.POST.get('hwindows', '')
        dic['sigma'] = request.POST.get('sigma', '')

        graficos_generados.clear()  # Limpiar antes de analizar
        graficos_generados = analyze(settings.MEDIA_ROOT, dic)  # Obtener gráficos

        dic['graficos_generados'] = graficos_generados
        request.session['dic'] = dic

    elif request.method == 'POST' and 'clear' in request.POST:
        dic = {}

    return render(request, 'baseidea.html', dic)


def user_guide(request):

    response = HttpResponse(content_type = 'pdf')
    response['Content-Disposition'] = 'attachment; filename = Toolbox_ringdown_python_draft.pdf'
    # CHECK ...
    return response


 
def csv_write(request):

    csv_name = 'idea_table_' + datetime.now().strftime("%d%b%Y_%Hh%Mm%Ss") + '.csv'

    response = HttpResponse(content_type = 'text/csv')
    response['Content-Disposition'] = 'attachment; filename = {}'.format(csv_name)

    writer = csv.writer(response)
    writer.writerow(['Ringdown Oscillations Monitoring and Analytics'])
    writer.writerow(['Date:', datetime.now().strftime("%d%b%Y_%Hh%Mm%Ss")])
    writer.writerow([])
    writer.writerow(['NAME', 'MODE', 'TYPE', 'FREQUENCY', 'AMPLITUDE', 'DAMPING', 'DAMPING RATIO', 'PHASE', 'ENERGY'])

    res_table = request.session['dic']['res_table']
    for i in res_table:
        writer.writerow(i)

    return response


def latex_code(request):

    def list_to_str(list):
        txt = ''
        for i in list:
            txt = txt + ' & ' + str(i)
        return txt[3:]

    text_name = 'idea_latex_' + datetime.now().strftime("%d%b%Y_%Hh%Mm%Ss") + '.txt'

    response = HttpResponse(content_type = 'text')
    response['Content-Disposition'] = 'attachment; filename = {}'.format(text_name)

    response.write(r'\begin{table*}[]' + '\n')
    response.write(r'  \centering' + '\n')
    response.write(r'  \caption{Dynamic parameters by IDEA}' + '\n')
    response.write(r'  \label{idea_table}' + '\n')
    response.write(r'  \begin{tabular}{lllllllll}' + '\n')
    response.write(r'      \hline' + '\n')
    response.write(r'      \hline' + '\n')
    response.write(r'      Name & Mode & Type & Frequency (Hz.) & Amplitude & Damping (1/s) & Damping ratio (\%) & Phase (rad) & Energy (\%)' + r' \\' + '\n')

    iname       = ''
    res_table   = request.session['dic']['res_table']
    for i in res_table:
        if iname != i[0]:
            response.write(r'        \hline' + '\n')
        iname = i[0]
        i = '       ' + list_to_str(i) + r' \\' + ' \n'
        response.write(i)

    response.write(r'      \hline' + '\n')
    response.write(r'      \hline' + '\n')
    response.write(r'  \end{tabular}' + '\n')
    response.write(r'\end{table*}' + '\n')

    return response
