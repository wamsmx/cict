import csv
import roma.plt_graphics        as plt_graph
import roma.plotly_graphics     as ply_graph

from datetime                   import datetime
from django.conf                import settings
from django.http                import HttpResponse
from django.shortcuts           import render
from django.core.files.storage  import FileSystemStorage
from roma.analyze_data          import analyze, get_signals_name, get_data


def roma(request):

    dic = {}


    if (request.method == 'POST') and ('upload' in request.POST):   
        uploaded_file = request.FILES['examinar']
        fs = FileSystemStorage()
        fs.save(uploaded_file.name, uploaded_file)
        dic['file_name']        = uploaded_file.name


    elif (request.method == 'POST') and ('preview' in request.POST):
        dic['file_name']        = request.POST['file_name']
        dic['csv_layout']       = request.POST['csv_layout']
        dic['method']           = request.POST['method']
        dic['isignals_name']    = request.POST['isignals_name']
        dic['modes_energy']     = request.POST['modes_energy']
        dic['downsampling']     = request.POST['downsampling']
        dic['nsignal']          = request.POST['nsignal']
        dic['time_start']       = request.POST['time_start']
        dic['time_end']         = request.POST['time_end']
        dic['fourier_start']    = request.POST['fourier_start']
        dic['fourier_end']      = request.POST['fourier_end']
        
        h_vec, t_vec, y_vec     = get_data(settings.MEDIA_ROOT, dic)

        dic['graph_preview']        = ply_graph.preview(h_vec, t_vec, y_vec, dic['time_start'], dic['time_end'])
        dic['graph_preview_mean']   = ply_graph.preview_mean(h_vec, t_vec, y_vec, dic['time_start'], dic['time_end'])
        

    elif (request.method == 'POST') and ('analyze' in request.POST):
        dic['file_name']        = request.POST['file_name']
        dic['csv_layout']       = request.POST['csv_layout']
        dic['method']           = request.POST['method']
        dic['isignals_name']    = request.POST['isignals_name']
        dic['modes_energy']     = request.POST['modes_energy']
        dic['downsampling']     = request.POST['downsampling']
        dic['nsignal']          = request.POST['nsignal']
        dic['time_start']       = request.POST['time_start']
        dic['time_end']         = request.POST['time_end']
        dic['fourier_start']    = request.POST['fourier_start']
        dic['fourier_end']      = request.POST['fourier_end']

        h_vec, t_vec, y_vec, dt, N, t_aprx, med, can, m_res, m_rts, l_mod, l_svn = analyze(settings.MEDIA_ROOT, dic)

        dic['graph_fft']        = ply_graph.fast_fourier_transform(dic, h_vec, y_vec, can, dt)
        dic['graph_sig']        = ply_graph.signal_reconstruction(h_vec, t_vec, y_vec, N, t_aprx, can, m_res, l_mod)
        dic['graph_pzr']        = ply_graph.poles_zeros(h_vec, can, m_res, m_rts)
        dic['res_table']        = ply_graph.res_table(h_vec, m_res, l_mod)

        if (dic['method'] == '2') or (dic['method'] == '3'):
            dic['graph_svn'] = ply_graph.svn(dic['nsignal'], h_vec, l_svn, l_mod)

        if dic['nsignal'] == '2':
            dic['graph_msh'] = ply_graph.mode_shapes(can, h_vec, m_res)
        
        request.session['dic'] = dic

    elif (request.method == 'POST') and ('clear' in request.POST):
        dic = {}

    return render(request, 'roma/index.html', dic)


def user_guide(request):

    response = HttpResponse(content_type = 'pdf')
    response['Content-Disposition'] = 'attachment; filename = Toolbox_ringdown_python_draft.pdf'
    # CHECK ...
    return response


 
def csv_write(request):

    csv_name = 'roma_table_' + datetime.now().strftime("%d%b%Y_%Hh%Mm%Ss") + '.csv'

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

    text_name = 'roma_latex_' + datetime.now().strftime("%d%b%Y_%Hh%Mm%Ss") + '.txt'

    response = HttpResponse(content_type = 'text')
    response['Content-Disposition'] = 'attachment; filename = {}'.format(text_name)

    response.write(r'\begin{table*}[]' + '\n')
    response.write(r'  \centering' + '\n')
    response.write(r'  \caption{Dynamic parameters by ROMA}' + '\n')
    response.write(r'  \label{roma_table}' + '\n')
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
