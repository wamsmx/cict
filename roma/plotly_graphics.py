
import os
import math

import numpy                    as np
import plotly.graph_objects     as go
import plotly.express           as px

from scipy.fftpack              import fft, fftshift

global font_size
font_size = 13
from matplotlib import cm
import matplotlib.colors as mcolors


def ncolors(n):

    if n == 1: n = 2
    colors=cm.get_cmap('rainbow', n)
    return [mcolors.rgb2hex(c) for c in colors(range(n))]
    #list = []
    #for i in range(n): list.append(i/(n-1))
    #return px.colors.sample_colorscale('Rainbow', list)


def fast_fourier_transform(dic, h_vec, y_vec, can, dt):

    y_z     = np.zeros([2**14 , 1]) 
    medf    = len(y_vec[:,0])

    S0 = np.zeros([len(y_z) + medf , can])
    Es = np.zeros([len(y_z) + medf , can], dtype = complex)

    for i in range(can):
        S0[: , i] = np.concatenate((y_vec[: , i], y_z), axis = None)
        Es[: , i] = fftshift(fft(S0[: , i]))

    fx  = np.linspace(-1.0 / (2.0 * dt), 1.0 / (2.0 * dt), len(S0[: , 0]))

    icolor = ncolors(can)
    fig = go.Figure()
    for i in range(can):
        fig.add_trace(go.Scatter(x = fx, y = abs(Es[: , i]), line={'width': 3,'color': icolor[i]},name = h_vec[i]))

    fig.update_layout(
        title           = '<b>Fourier spectra</b>',
        height          = 550,
        template        = 'simple_white',
        xaxis_title     = 'Frequency (Hz.)',
        yaxis_title     = 'Magnitude',
        legend_title    = 'Signal',
        xaxis_range     = [float(dic['fourier_start']), float(dic['fourier_end'])],
        
        xaxis = dict(
            showgrid=True
            ),

        yaxis = dict(
            showgrid=True,
            ),

        font = dict(
            # family  = 'Liberation Serif',
            size    = font_size,
            color   = 'black',
        )
    )

    return fig.to_html()


def svn(nsignal, h_vec, l_svn, l_mod):

    if nsignal == '1': n = len(h_vec)
    if nsignal == '2': n = 1

    icolor  = ncolors(n)
    fig     = go.Figure()
    for i in range(n):
        fig.add_trace(go.Scatter(
            y           = np.log10(l_svn[i]),
            line        = {'width': 3,'dash': 'solid', 'color': icolor[i]},
            legendgroup = str(i),
            name        = h_vec[i] + ': (' + str(int(round(l_mod[i]/2.0, 0))) + ')'
            )
        )
        fig.add_vline(
            name        = h_vec[i],
            x           = l_mod[i],
            line_width  = 3,
            line_dash   = 'dash',
            line_color  = icolor[i],
            opacity     = 1.0
            ),

    fig.update_layout(
        height          = 550,
        template        = 'simple_white',
        title           = '<b>SVD from Hankel matrix</b>',
        xaxis_title     = 'Singular value number',
        yaxis_title     = 'log10(\u03C3)',
        legend_title    = 'Signal: (max modes)',

        xaxis = dict(
            showgrid=True
            ),

        yaxis = dict(
            showgrid=True,
            ),

        font = dict(
            # family  = 'Liberation Serif',
            size    = font_size,
            color   = 'black',
            ),
        )

    return fig.to_html()


def signal_reconstruction(h_vec, t_vec, y_vec, N, t_aprx, can, m_res, l_mod):

    mag     = m_res[:, :, 0]
    ang     = m_res[:, :, 1]
    damp    = m_res[:, :, 2]
    freq    = m_res[:, :, 3]

    y_aprx = np.zeros([N, can])
    for ican, i in enumerate(l_mod):
        for j in range(i):
            if freq[j][ican] > 0:
                imag    = mag[j][ican]
                iang    = ang[j][ican]
                idamp   = damp[j][ican]
                ifreq   = freq[j][ican]

                y_aprx[: , ican] = y_aprx[: , ican] + (
                    imag * np.exp(idamp * t_aprx) * np.cos((2 * np.pi * ifreq) * t_aprx + iang)
                        )

    icolor = ncolors(can)
    fig     = go.Figure()
    for i in range(len(l_mod)):

        s0      = np.linalg.norm(list(y_vec[: , i]))
        s1      = np.linalg.norm(list(y_aprx[: , i]))
        sfull2  = pow(s0 , 2)
        e2      = pow(s0 - s1 , 2)
        SNR     = round(10 * math.log10(sfull2 / e2) , 2)

        fig.add_trace(go.Scatter(
            x = t_vec, y = y_vec[: , i],
            line={'width': 3,'dash': 'solid', 'color': icolor[i]},
            legendgroup = str(i),
            name = h_vec[i] + ': ' + str(SNR) + 'dB')
        )
        fig.add_trace(go.Scatter(
            x = t_vec, y = y_aprx[: , i],
            line={'width': 3,'dash': 'dash', 'color': icolor[i]},
            legendgroup = str(i),
            name = h_vec[i] + '_aprox',
            showlegend = False)
        )

    fig.update_layout(
        height          = 550,
        template        = 'simple_white',
        title           = '<b>Signal reconstruction</b>',
        xaxis_title     = 'Time (s)',
        yaxis_title     = 'Magnitude',
        legend_title    = 'Signal: Noise Ratio',

        xaxis = dict(
            showgrid=True
            ),

        yaxis = dict(
            showgrid=True,
            ),

        font = dict(
            # family  = 'Liberation Serif',
            size    = font_size,
            color   = 'black',
            ),
    )

    return fig.to_html()

    
def poles_zeros(h_vec, can, m_res, m_rts):

    freq    = m_res[: , : , 3]
    roots   = m_rts[: , :]

    icolor = ncolors(can)
    fig = go.Figure()
    for ican in range(can):
        xdata = []
        ydata = []

        for i in range(freq.shape[0]):
            if (freq[i , ican] != 0) and (freq[i , ican] != -0.0001):
                xdata.append(roots[i , ican].real)
                ydata.append(roots[i , ican].imag)

        fig.add_trace(
            go.Scatter(
                x       = xdata,
                y       = ydata,
                name    = h_vec[ican],
                mode    = 'markers',
                marker  = dict(symbol = 'circle-dot', size = 9, color = 'rgba(0,0,0,0)', 
                    line = dict(width   = 2, color = icolor[ican])
                    )
                )
            )

    fig.add_shape(
        type        = 'circle',
        fillcolor   = 'rgba(0,0,0,0)',
        xref        = 'x',
        yref        = 'y',
        x0          = -1, 
        y0          = -1, 
        x1          = 1, 
        y1          = 1,
        line = dict(width   = 2, color = 'black')
        )

    fig.update_xaxes(
        range           =[-1.2, 1.2],
        zeroline        = True,
        zerolinewidth   = 0.5,
        zerolinecolor   = 'grey',
        )

    fig.update_yaxes(
        range           =[-1.2, 1.2],
        zeroline        = True,
        zerolinewidth   = 0.5,
        zerolinecolor   = 'grey',
        )

    fig.update_layout(
        height          = 600,
        width           = 600,
        template        = 'simple_white',
        title           = '<b>Pole-Zero</b>',
        xaxis_title     = 'Re',
        yaxis_title     = 'Im',
        legend_title    = 'Signal',

        xaxis = dict(showgrid = True),
        yaxis = dict(showgrid = True),

        font = dict(
            # family  = 'Liberation Serif',
            size    = font_size,
            color   = 'black',
            ),
        )

    return fig.to_html()


def mode_shapes(can, h_vec, m_res):

    mag     = m_res[:, :, 0]
    ang     = m_res[:, :, 1]
    freq    = m_res[:, :, 3]

    fig = go.Figure()
    icolor = ncolors(can)

    m = 0
    for j in range(freq.shape[0]):
        if freq[j , 1] > 0:
            m = m + 1
            
            imag = list(mag[j,:])
            iang = list(ang[j,:])

            indx_mag_max = imag.index(max(imag))

            imag_ref = imag / max(imag)
            iang_ref = (iang - iang[indx_mag_max]) * (180.0 / np.pi)

            for i in range(freq.shape[1]):

                ir  = [0, imag_ref[i]]
                ith = [0, iang_ref[i]]

                fig.add_trace(
                    go.Scatterpolar(
                        r = ir, theta = ith,
                        mode = 'lines+markers',
                        line = {'width': 3, 'dash': 'solid', 'color': icolor[i]},
                        legendgroup = str(j),
                        name = h_vec[i] + ': mode ' + str(m)
                    )
                )

    fig.update_layout(
        height          = 600,
        width           = 600,
        template        = 'simple_white',
        title           = '<b>Mode-shapes</b>',
        legend_title    = 'Signal: mode',
        
        polar = dict(
            angularaxis = dict(
                showgrid    = True,
            ),
            radialaxis = dict(
                range           = [0, 1],
                showticklabels  = False,
                showgrid        = True
            ),
        ),

        font = dict(
            # family  = 'Liberation Serif',
            size    = font_size,
            color   = 'black',
        ),
    )

    return fig.to_html()


def res_table(h_vec, m_res, l_mod):

    mag     = m_res[:, :, 0]
    ang     = m_res[:, :, 1]
    damp    = m_res[:, :, 2]
    freq    = m_res[:, :, 3]
    drat    = m_res[:, :, 4]
    enrg    = m_res[:, :, 5]

    res_table = []
    for ican, i in enumerate(l_mod):
        imod = 0
        for j in range(i):
            if freq[j][ican] > 0:
                imod    = imod + 1
                imag    = mag[j][ican]
                iang    = ang[j][ican]
                idamp   = damp[j][ican]
                ifreq   = freq[j][ican]
                idrat   = drat[j][ican]
                ienrg   = enrg[j][ican]

                if ifreq >= 0.10 and ifreq < 0.80:
                    osc_mod = 'Inter-Area'
                elif ifreq >= 0.80 and ifreq < 2.00:
                    osc_mod = 'Local'
                elif ifreq >= 2.00 and ifreq < 3.00:
                    osc_mod = 'Intra-Plant'
                else:
                    osc_mod = '-'

                res_table.append([h_vec[ican], imod, osc_mod, np.round(ifreq, 3), np.round(imag, 3), np.round(idamp, 3), np.round(idrat, 3), np.round(iang, 3), np.round(ienrg,4)])

    return res_table


def preview(h_vec, t_vec, y_vec, time_start, time_end):

    icolor  = ncolors(len(h_vec))
    fig     = go.Figure()
    for i in range(len(h_vec)):

        fig.add_trace(go.Scatter(
            x = t_vec, y = y_vec[: , i],
            line = {'width': 3, 'color': icolor[i]},
            legendgroup = str(i),
            name = h_vec[i]
            )
        )

    fig.add_vrect(x0 = float(time_start),
        x1              = float(time_end),
        line_width      = 2,
        line_dash       = 'dash',
        line_color      = 'black',
        fillcolor       = 'white',

        annotation_text         = 'Window analysis: ' + str(float(time_start)) + ' to ' + str(float(time_end)) + 's.',
        annotation_font_size    = 16,
        # annotation_font_family  = 'Liberation Serif',

        )

    fig.update_layout(
        height          = 550,
        template        = 'simple_white',
        title           = '<b>Input signals preview</b>',
        xaxis_title     = 'Time (s)',
        yaxis_title     = 'Magnitude',
        legend_title    = 'Signal:',

        xaxis = dict(
            showgrid=True
            ),

        yaxis = dict(
            showgrid=True,
            ),

        font = dict(
            # family  = 'Liberation Serif',
            size    = font_size,
            color   = 'black',
            ),
        )

    return fig.to_html()


def preview_mean(h_vec, t_vec, y_vec, time_start, time_end):

    for i in range(y_vec.shape[1]):
        y_vec[: , i] = y_vec[: , i] - np.mean(y_vec[: , i])

    icolor  = ncolors(len(h_vec))
    fig     = go.Figure()
    for i in range(len(h_vec)):

        fig.add_trace(go.Scatter(
            x = t_vec, y = y_vec[: , i],
            line = {'width': 3, 'color': icolor[i]},
            legendgroup = str(i),
            name = h_vec[i]
            )
        )

    fig.add_vrect(x0 = float(time_start),
        x1              = float(time_end),
        line_width      = 2,
        line_dash       = 'dash',
        line_color      = 'black',
        fillcolor       = 'white',

        annotation_text         = 'Window analysis: ' + str(float(time_start)) + ' to ' + str(float(time_end)) + 's.',
        annotation_font_size    = 16,
        # annotation_font_family  = 'Liberation Serif',

        )

    fig.update_layout(
        height          = 550,
        template        = 'simple_white',
        title           = '<b>Input signals preview</b> (without mean value)',
        xaxis_title     = 'Time (s)',
        yaxis_title     = 'Magnitude',
        legend_title    = 'Signal:',

        xaxis = dict(
            showgrid=True
            ),

        yaxis = dict(
            showgrid=True,
            ),

        font = dict(
            # family  = 'Liberation Serif',
            size    = font_size,
            color   = 'black',
            ),
        )

    return fig.to_html()
