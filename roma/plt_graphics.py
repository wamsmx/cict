
import os
import math
import numpy 				as np
import matplotlib.pyplot 	as plt

from io 					import StringIO
from scipy.fftpack          import fft, fftshift


def fast_fourier_transform(dic, h_vec, y_vec, can, dt):

    y_z     = np.zeros([2**14 , 1]) 
    medf    = len(y_vec[:,0])

    S0 = np.zeros([len(y_z) + medf , can])
    Es = np.zeros([len(y_z) + medf , can], dtype = complex)

    for i in range(can):
        S0[: , i] = np.concatenate((y_vec[: , i], y_z), axis = None)
        Es[: , i] = fftshift(fft(S0[: , i]))

    fx  = np.linspace(-1.0 / (2.0 * dt), 1.0 / (2.0 * dt), len(S0[: , 0]))
    fig = plt.Figure(dpi = 100) 
    
    ax = fig.add_subplot(111)
    ax.set_title('Fourier Spectrum')
    ax.set_xlabel('Frequency (Hz.)')
    ax.set_ylabel('Magnitude')
    ax.set_xlim(float(dic['fourier_start']), float(dic['fourier_end']))    

    colors = plt.cm.get_cmap('gist_rainbow', can+1)
    for i in range(can):
        ax.plot(fx, abs(Es[: , i]) / (2**14), color = colors(i), label = h_vec[i], lw = '2')
    
    ax.legend(loc = 'upper right')
    ax.grid(True, color = 'grey',  linestyle = ':', linewidth = 0.75)

    imgdata = StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)
    data = imgdata.getvalue()
    
    return data


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

    fig = plt.Figure(dpi = 100) 
    ax  = fig.add_subplot(111)
    ax.set_title('Signal Recosntruction')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Magnitude')

    for i in range(len(l_mod)):

        colors  = plt.cm.get_cmap('gist_rainbow', can+1)
        s0      = np.linalg.norm(list(y_vec[: , i]))
        s1      = np.linalg.norm(list(y_aprx[: , i]))
        sfull2  = pow(s0 , 2)
        e2      = pow(s0 - s1 , 2)
        SNR     = round(10 * math.log10(sfull2 / e2) , 2)

        ax.plot(t_vec, y_vec[: , i], color = colors(i), label = h_vec[i] + ': ' + str(SNR), lw = '2')
        ax.plot(t_vec, y_aprx[: , i], dashes = [6,4], color = colors(i), lw = '2')


    ax.legend(loc = 'upper right', title = 'Signal Noise Ratio [dB]')
    ax.grid(True, color = 'grey',  linestyle = ':', linewidth = 0.75)

    imgdata = StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)
    data = imgdata.getvalue()
    
    return data


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




