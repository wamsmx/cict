# Paso 1: Configurar el backend de Matplotlib para evitar problemas en entornos sin GUI
import matplotlib
matplotlib.use('Agg')  # Backend sin GUI

# 📌 Librerías estándar de Python
import os
import sys
import csv
import io
import json
import base64
from datetime import datetime, timezone, timedelta

import numpy as np
import pandas as pd
from scipy.fft import fft, fftfreq
from scipy.linalg import hankel, svd, pinv, eig

import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.io import to_json
import mpld3


graficos_generados = []


def get_data(media_path, dic):

	csv_file = open(os.path.join(media_path, dic['file_name']))
	csv_data = csv.reader(csv_file)

	csv_list = []
	for row in csv_data:
		csv_list.append(row)

	h_vec   = csv_list[0][1:]
	t, y    = [], []
	for index, icsv_list in enumerate(csv_list[1:]):
		t.append(float(icsv_list[0]))
		y.append([])
		for jcsv_list in icsv_list[1:]:
			y[index].append(float(jcsv_list))

	y_vec = np.array(y)
	t_vec = np.array(t)
	t_vec = t_vec - t_vec[0]

	if (dic['downsampling'] != '') and (int(dic['downsampling']) > 1):
		t_vec, y_vec = downsampling(int(dic['downsampling']), t_vec, y_vec)

	if (dic['isignals_name'] != ''):
		h_vec, y_vec = select_signals(h_vec, y_vec, dic)

	return h_vec, t_vec, y_vec


def s_prony(Y, N, dt, fct_method):

	modes = int(fct_method)
	if (modes == 0) or (modes > int(Y.shape[0] / 2)):
		modes = int(Y.shape[0] / 2)

	T = np.zeros((N - modes, modes))
	for i in range(modes):
		T[:, i] = Y[modes-1-i : N-1-i]

	b = Y[modes : N]
	a = np.linalg.inv(T.T @ T) @ T.T @ b

	z0 = [1]
	for i in a:
		z0.append(-i)

	z = np.roots(z0)
	Z = np.zeros((N, modes), dtype = complex)

	for i in range(N):
		Z[i, :] = pow(z, i)

	lam = np.log(z)/dt

	B = np.linalg.inv(Z.T @ Z) @ Z.T @ Y

	mag 	= 2.0 * abs(B)
	ang 	= np.angle(B)
	damp 	= lam.real
	freq 	= lam.imag / (2.0 * np.pi)

	freq_zrs = list(np.where(freq == 0.0)[0])
	for fz in freq_zrs:
		freq[fz] = -0.0001

	omga	= 2.0 * np.pi * freq
	damprat	= (damp / omga) * 100.0
	enrgy 	= (1.0 / 2.0) * (omga**2)  * (mag**2) * 100.0
	roots	= z

	return modes, mag, ang, damp, freq, damprat, enrgy, roots, 0


def s_mp(Y, N, dt, fct_method):

	r = int(np.around((N/2.0)-1.0, 0))
	Hankel = hankel(Y)
	H0 = Hankel[0:r, 0:r]
	H1 = Hankel[1:r+1, 0:r]

	u, s0, v0 = np.linalg.svd(H0, full_matrices=False)
	v = v0.T

	sum_s = sum(s0)
	sum_st = 0.0
	for modes, i in enumerate(s0, start = 1):
		sum_st 		= sum_st + i
		pc_sum_st 	= (sum_st / sum_s) * 100.0
		if pc_sum_st > float(fct_method): break

	V1 = np.zeros((r - 1 , modes))
	V2 = np.zeros((r - 1 , modes))

	for i in range(modes):
		V1[: , i] = v[0 : r - 1 , i]
		V2[: , i] = v[1 : r , i]

	Y1 = np.dot(V1.T , V1)
	Y2 = np.dot(V2.T , V1)

	z = np.linalg.inv(Y1)
	z = np.dot(z , Y2)
	z, b = np.linalg.eig(z)

	# MPM finished, continue as Prony method...
	lam = np.log(z) / dt

	Z = np.zeros((N, modes), dtype = complex)
	for i in range(N):
		Z[i, :] = pow(z, i)

	B = np.linalg.inv(Z.T @ Z) @ Z.T @ Y

	mag 	= 2.0 * abs(B)
	ang 	= np.angle(B)
	damp 	= lam.real
	freq 	= lam.imag / (2.0 * np.pi)

	freq_zrs = list(np.where(freq == 0.0)[0])
	for fz in freq_zrs:
		freq[fz] = -0.0001

	omga	= 2.0 * np.pi * freq
	damprat	= (damp / omga) * 100.0
	enrgy 	= (1.0 / 2.0) * (omga**2)  * (mag**2) * 100.0
	roots	= z

	return modes, mag, ang, damp, freq, damprat, enrgy, roots, s0


def s_era(Y, N, dt, fct_method):

	r = int(np.around((N/2.0)-1.0, 0))
	Hankel = hankel(Y)
	H0 = Hankel[0:r, 0:r]
	H1 = Hankel[1:r+1, 0:r]

	u, s0, v0 = np.linalg.svd(H0, full_matrices=False)
	s = np.diag(s0)
	v = v0.T

	sum_s 	= np.sum(s0)
	sum_st 	= 0.0
	for modes, i in enumerate(s0,  start = 1):
		sum_st 		= sum_st + i
		pc_sum_st 	= (sum_st / sum_s) * 100.0
		if pc_sum_st > float(fct_method): break

	U = u[: , 0:modes]
	S = s[0:modes , 0:modes]
	V = v[: , 0:modes]

	Sr = np.zeros((modes , modes))
	for i in range(modes):
		Sr[i,i] = pow(S[i,i], -0.5)

	A1 = Sr @ U.T @ H1 @ V @ Sr

	z, b = np.linalg.eig(A1)

	Z = np.zeros((N, modes), dtype = complex)
	for i in range(N):
		Z[i, :] = pow(z, i)

	# ERA finished, continue as Prony method...
	lam = np.log(z) / dt

	Z = np.zeros((N, modes), dtype = complex)
	for i in range(N):
		Z[i, :] = pow(z, i)

	B = np.linalg.inv(Z.T @ Z) @ Z.T @ Y

	mag 	= 2.0 * abs(B)
	ang 	= np.angle(B)
	damp 	= lam.real
	freq 	= lam.imag / (2.0 * np.pi)

	freq_zrs = list(np.where(freq == 0.0)[0])
	for fz in freq_zrs:
		freq[fz] = -0.0001

	omga	= 2.0 * np.pi * freq
	damprat	= (damp / omga) * 100.0
	enrgy 	= (1.0 / 2.0) * (omga**2)  * (mag**2) * 100.0
	roots	= z

	return modes, mag, ang, damp, freq, damprat, enrgy, roots, s0


def m_era(Y, N, dt, fct_method):

	med 	= Y.shape[0]
	can 	= Y.shape[1]
	r 		= int(np.around((N / 2.0) - 1.0, 0))

	H0 = np.zeros([r * can, r])
	H1 = np.zeros([r * can, r])
	for j in range(r):
		for i in range(r):
			H0[can*i: can*(i+1) , j] = Y[j+i+1 , :]
			H1[can*i: can*(i+1) , j] = Y[j+i+2 , :]

	u, s0, v0 = np.linalg.svd(H0, full_matrices=True)
	s = np.diag(s0)
	v = v0.T

	sum_s = sum(s0)
	sum_st = 0.0
	for modes, i in enumerate(s0, start = 1):
		sum_st 		= sum_st + i
		pc_sum_st 	= (sum_st / sum_s) * 100.0
		if pc_sum_st > float(fct_method): break
	
	U = u[: , 0:modes]
	S = s[0:modes , 0:modes]
	V = v[: , 0:modes]

	Sr = np.zeros((modes , modes))
	for i in range(modes):
		Sr[i,i] = pow(S[i,i], -0.5)

	A = Sr @ U.T @ H1 @ V @ Sr

	eigA, z = np.linalg.eig(A)
	lam = np.log(eigA)/dt

	Z = np.zeros((N, modes), dtype = complex)
	for i in range(N):
		Z[i, :] = pow(eigA, i)

	B = np.dot(np.linalg.pinv(Z), Y)

	# Results ------------------------------------------------->
	roots0 	= eigA
	damp0 	= lam.real
	freq0 	= lam.imag / (2.0 * np.pi)
	mag 	= 2 * abs(B)
	ang 	= np.angle(B)

	# Change any zero in frequency vector by -0.0001
	zrs = list(np.where(freq0 == 0.0)[0])
	for z in zrs: freq0[z] = -0.0001

	omga		= 2.0 * np.pi * freq0
	damprat0	= (damp0 / omga) * 100.0

	enrgy = np.zeros((len(freq0) , can))
	for i in range(len(freq0)):
		enrgy[i , :] = (1.0 / 2.0) * (omga[i]**2)  * (mag[i , :]**2) * 100.0


	freq 	= np.zeros([modes , can])
	damp 	= np.zeros([modes , can])
	damprat = np.zeros([modes , can])
	roots 	= np.zeros([modes , can], dtype = complex)

	for i in range(can):
		damp[: , i] 	= damp0
		freq[: , i] 	= freq0
		damprat[: , i] 	= damprat0
		roots[: , i] 	= roots0

	return modes, mag, ang, damp, freq, damprat, enrgy, roots, s0


def ringdown(h_vec, t_vec, y_vec, dic):

	for i in range(y_vec.shape[1]):
		y_vec[: , i] = y_vec[: , i] - np.mean(y_vec[: , i])

	if dic['method'] == '1':
		fct_method = int(dic['modes_energy'])
	else:
		fct_method = float(dic['modes_energy'])

	t_start 		= float(dic['time_start'])
	t_end 			= float(dic['time_end'])

	dt = float(t_vec[1] - t_vec[0])
	pa = int(round(t_start / dt, 0))
	pb = int(round(t_end / dt, 0)) + int(1)

	t_vec 	= t_vec[pa : pb]
	y_vec 	= y_vec[pa : pb , :]
	N 		= len(t_vec)
	t_aprx 	= np.linspace(0, t_end - t_start, N)
	med 	= y_vec.shape[0]
	can 	= y_vec.shape[1]

	L_res 	= []
	l_mod 	= []
	l_svn 	= []
	l_isvn 	= []
	if dic['nsignal'] == '1':

		for c in range(can):

			Y = np.zeros([y_vec.shape[0], 1])
			Y = (y_vec.T)[c]

			if dic['method'] == '1': # PRONY
				modes, mag, ang, damp, freq, damprat, enrgy, roots, svn = s_prony(Y, N, dt, fct_method)

			if dic['method'] == '2': # MATRIX PENCIL
				modes, mag, ang, damp, freq, damprat, enrgy, roots, svn = s_mp(Y, N, dt, fct_method)

			if dic['method'] == '3': # ERA
				modes, mag, ang, damp, freq, damprat, enrgy, roots, svn = s_era(Y, N, dt, fct_method)

			l_svn.append(svn)
			L_res.append([mag, ang, damp, freq, damprat, enrgy, roots])
			l_mod.append(modes)

	if dic['nsignal'] == '2':

		Y = y_vec
		if dic['method'] == '3': # ERA
			modes, mag, ang, damp, freq, damprat, enrgy, roots, svn = m_era(Y, N, dt, fct_method)

			l_svn.append(svn)
			L_res.append([mag, ang, damp, freq, damprat, enrgy, roots])
			l_mod.append(modes)
			l_mod = l_mod * can

	max_mod = max(l_mod)
	m_res 	= np.zeros([max_mod, can, 7])
	m_rts 	= np.zeros([max_mod, can], dtype = complex)

	if dic['nsignal'] == '1':
		for c in range(can):
			for m in range(l_mod[c]):
				m_res[m , c , 0] 	= L_res[c][0][m]	# mag
				m_res[m , c , 1] 	= L_res[c][1][m]	# ang
				m_res[m , c , 2] 	= L_res[c][2][m]	# damp
				m_res[m , c , 3] 	= L_res[c][3][m]	# freq
				m_res[m , c , 4] 	= L_res[c][4][m]	# damprat
				m_res[m , c , 5] 	= L_res[c][5][m]	# enrgy
				m_rts[m , c] 		= L_res[c][6][m]	# roots

	if dic['nsignal'] == '2':
		m_res[: , : , 0] 	= mag
		m_res[: , : , 1] 	= ang
		m_res[: , : , 2] 	= damp
		m_res[: , : , 3] 	= freq
		m_res[: , : , 4] 	= damprat
		m_res[: , : , 5] 	= enrgy
		m_rts[: , :] 		= roots

	return h_vec, t_vec, y_vec, dt, N, t_aprx, med, can, m_res, m_rts, l_mod, l_svn

def downsampling(m, t_vec, y_vec):

	t, y = [], []

	for i in range(0 , t_vec.shape[0] , m):
		t.append(t_vec[i])
		y.append(y_vec[i])

	return np.array(t), np.array(y)


def str_to_list(str_signals):

	str_signals = str_signals.replace(' ', '')

	if len(str_signals) == 0: return False

	if str_signals[0] 	== ',': str_signals = str_signals[1:]
	if str_signals[-1] 	== ',': str_signals = str_signals[:-1]

	list_signals = []

	a = 0
	for i in range(len(str_signals)):
		if str_signals[i] == ',':
			list_signals.append(str_signals[a : i])
			a = i + 1

	list_signals.append(str_signals[a : len(str_signals)])

	return list_signals


def select_signals(h_vec, y_vec, dic):

	list_signals = str_to_list(dic['isignals_name'])

	for i in range(len(h_vec), 0, -1):
		if h_vec[i-1] not in list_signals:
			y_vec = np.delete(y_vec, (i - 1), axis = 1)

	return list_signals, y_vec


def extract_time_and_frequency(media_path, dic):
    # Leer el archivo Excel
	csv_file = open(os.path.join(media_path, dic['file_name']))
	df = pd.read_csv(csv_file)
    # Extraer la primera columna (Tiempo) excepto el primer valor
	time_values = df.iloc[1:, 0].tolist()
    
    # Extraer las señales de frecuencia (todas las demás columnas)
	freq_values = [df.iloc[1:, i].tolist() for i in range(1, df.shape[1])]
    
    # Crear v_time con la misma cantidad de listas que señales de frecuencia
	v_time = [time_values] * len(freq_values)
	v_freq = freq_values
    
	return v_time, v_freq

def create_data_windows(freqarre, window, hwindow):
    """
    Crea ventanas deslizantes a partir de una serie de datos.
    
    Parámetros:
    - freqarre: Array de datos de frecuencia.
    - window: Tamaño de la ventana.
    - hwindow: Paso de desplazamiento.

    Retorna:
    - Array de ventanas de datos.
    """
    cols = len(freqarre)
    data_windows = [freqarre[j:j + window] for j in range(0, cols - window + 1, hwindow)]
    return np.array(data_windows)

def detect_valid_windows(data_windows, method, sigma, Fs, LFFT, Tolerance, iteracion):
    """
    Detecta ventanas válidas usando tres métodos de análisis.
    
    Parámetros:
    - data_windows: Matriz con ventanas de datos.

    Retorna:
    - Lista con índices de ventanas válidas.
    """
    if method == 0:
        ventanas_fuera, _ = matrix_pencil(data_windows, sigma, Fs, Tolerance, method, iteracion)
    elif method == 1:
        ventanas_fuera, _ = fft_events(data_windows, sigma, Fs, LFFT, method, iteracion)
    elif method == 2:
        ventanas_fuera, _ = min_max(data_windows, sigma, method, iteracion)
    #return sorted(set(ventanas_fuera_mm) & set(ventanas_fuera_fft) & set(ventanas_fuera_mp))
    return sorted(ventanas_fuera)

def filter_non_consecutive(ventanas_validas):
    """
    Filtra índices consecutivos, dejando solo los principales eventos detectados.

    Parámetros:
    - ventanas_validas: Lista de índices de ventanas.

    Retorna:
    - Lista de índices filtrados.
    """
    return [ventanas_validas[i] for i in range(len(ventanas_validas)) if i == 0 or ventanas_validas[i] != ventanas_validas[i - 1] + 1]

import matplotlib
matplotlib.use('Agg')  # Para evitar problemas en entornos sin GUI

import matplotlib.pyplot as plt
import numpy as np
import io
import base64

def plot_results1(timet, freqarre, ventanas_validas, resultado, windows, hwindows, iteracion):
    """
    Genera un gráfico profesional de frecuencia vs tiempo, destacando solo el primer evento detectado.

    Parámetros:
    - timet: Lista de tiempos.
    - freqarre: Lista de frecuencias.
    - ventanas_validas: Índices de ventanas válidas.
    - resultado: Índices de eventos detectados.
    - windows: Tamaño de la ventana de análisis.
    - hwindows: Separación entre ventanas.
    - iteracion: Número de iteración actual.
    """

    # 📌 Validar si hay eventos detectados
    if not resultado or len(resultado) == 0:
        print("⚠️ No hay eventos detectados, no se generará gráfico.")
        return

    # 📌 Tomar solo el primer evento detectado
    first_event_idx = resultado[0]
    if first_event_idx + (windows - 1) >= len(timet):
        print("⚠️ El primer evento está fuera del rango de la señal. No se generará gráfico.")
        return

    event_x = timet[first_event_idx * hwindows + int(windows/2) - 1]  # Tiempo exacto del evento
    event_time_str = f"{event_x:.2f} s"  # Formato para mostrar en el gráfico

    # 📌 Crear figura
    fig, ax = plt.subplots(figsize=(10, 5), dpi=100)

    # 📌 Graficar la señal completa pero enfocándose en el evento
    ax.plot(timet, freqarre, label="Frecuencia", color="royalblue", linewidth=2)

    # 📌 Resaltar solo el primer evento detectado
    ax.axvline(x=event_x, color='red', linestyle='--', linewidth=2, label="Evento detectado")
    ax.annotate(f"Evento en {event_time_str}",
                xy=(event_x, max(freqarre)), 
                xytext=(event_x, max(freqarre) + 0.2), 
                arrowprops=dict(facecolor='red', shrink=0.05), 
                fontsize=10, color="red")

    # 📌 Ajustar el eje X para centrarse en el evento detectado
    x_margin = 5  # Margen de 5 segundos antes y después del evento
    x_min = max(min(timet), event_x - x_margin)
    x_max = min(max(timet), event_x + x_margin)
    ax.set_xlim(x_min, x_max)

    # 📌 Configurar ejes y estilos
    ax.set_xlabel('Tiempo (s)', fontsize=12)
    ax.set_ylabel('Frecuencia (Hz)', fontsize=12)
    ax.set_title(f'Primer Evento Analizado - Iteración {iteracion}', fontsize=14, fontweight='bold')
    ax.grid(True, linestyle=":", linewidth=0.75, alpha=0.7)

    # 📌 Mejorar leyenda
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = dict(zip(labels, handles))
    ax.legend(unique_labels.values(), unique_labels.keys(), loc="upper right", fontsize=10, frameon=True)

    # 📌 Ajuste de márgenes para mejor visualización
    plt.tight_layout()

    # 📌 Guardar imagen en formato base64 para Django
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', bbox_inches="tight", dpi=100)
    buffer.seek(0)
    image_png = buffer.getvalue()
    buffer.close()
    graphic = base64.b64encode(image_png).decode('utf-8')

    # 📌 Guardar en la lista para visualización en Django
    graficos_generados.append((f'Evento {iteracion}', graphic))

    # 📌 Cerrar figura para liberar memoria
    plt.close(fig)
	
def plot_results2(amplitudes_maximas, limite_inferior, limite_superior, method, iteracion):
    """
    Genera gráficos de amplitudes con límites y diferentes métodos de análisis.

    Parámetros:
    - amplitudes_maximas: Datos de amplitudes máximas.
    - limite_inferior: Límite inferior de detección.
    - limite_superior: Límite superior de detección.
    - method: Método de análisis (0 = MP, 1 = FFT, 2 = Min Max).
    - iteracion: Número de iteración actual.
    """

    # Definir colores y etiquetas válidas
    colors = ['green', 'blue', 'black']
    labels = ['Matrix Pencil', 'FFT', 'Min Max']

    # Validar que method esté en el rango correcto
    if not isinstance(method, int) or method < 0 or method >= len(labels):
        method = 0  # Asigna un método por defecto en caso de error


    # Crear figura
    fig, ax = plt.subplots(figsize=(10, 5), dpi=100)

    # Graficar amplitudes con el método correspondiente
    ax.scatter(range(len(amplitudes_maximas)), amplitudes_maximas, color=colors[method], label=labels[method], s=50)

    # Agregar líneas de límites
    ax.axhline(y=limite_superior, color=colors[method], linestyle='--', linewidth=2, label=f'{labels[method]} ±3σ')
    ax.axhline(y=limite_inferior, color=colors[method], linestyle='--', linewidth=2)

    # Configurar ejes
    ax.set_xlabel('Índice de Ventana', fontsize=12)
    ax.set_ylabel('Amplitud / Frecuencia', fontsize=12)
    ax.set_title(f'Análisis - Método {labels[method]} - Iteración {iteracion}', fontsize=14, fontweight='bold')

    # Mejorar visualización
    ax.grid(True, linestyle=":", linewidth=0.75, alpha=0.7)

    # Optimizar leyenda sin duplicados
    handles, labels_leg = ax.get_legend_handles_labels()
    unique_labels = dict(zip(labels_leg, handles))
    ax.legend(unique_labels.values(), unique_labels.keys(), loc="upper right", fontsize=10, frameon=True)

    # Ajuste de márgenes
    plt.tight_layout()

    # Convertir a Base64 para Django
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', bbox_inches="tight", dpi=100)
    buffer.seek(0)
    image_png = buffer.getvalue()
    buffer.close()
    graphic = base64.b64encode(image_png).decode('utf-8')

    # Verificar que graficos_generados está inicializado
    global graficos_generados

    # Guardar en la lista
    graficos_generados.append((f'Método {labels[method]} - {iteracion}', graphic))

    # Cerrar la figura para liberar memoria
    plt.close(fig)

def voltage_angle_difference(timet1, freqarre1, dic):
    """
    Función principal que detecta eventos de diferencia de ángulo de voltaje.

    Parámetros:
    - timet1: Lista de listas de tiempos.
    - freqarre1: Lista de listas de frecuencias.

    Retorna:
    - event_detect: Lista de eventos detectados.
    - ventanas_freq1: Lista de ventanas de frecuencia relevantes.
    - ventanas_time1: Lista de ventanas de tiempo relevantes.
    """
    try:
        Fs = int(dic['fs'])
    except (ValueError, KeyError) as e:
        print(f"Error al convertir 'fs': {e}. Asignando valor predeterminado 0.")
        Fs = 60

    try:
        LFFT = int(dic['lfft'])
    except (ValueError, KeyError) as e:
        print(f"Error al convertir 'lfft': {e}. Asignando valor predeterminado 0.")
        LFFT = 1024

    try:
        Tolerance = float(dic['tolerance'])
    except (ValueError, KeyError) as e:
        print(f"Error al convertir 'tolerance': {e}. Asignando valor predeterminado 0.0.")
        Tolerance = 1e-6

    event_detect = []
    ventanas_freq1 = []
    ventanas_time1 = []

    # Asegurándote de que 'windows', 'hwindows', 'sigma', y 'method' sean enteros
    try:
        windows = int(dic['windows'])
    except (ValueError, KeyError) as e:
        print(f"Error al convertir 'windows': {e}. Asignando valor predeterminado 0.")
        windows = 10

    try:
        hwindow = int(dic['hwindows'])
    except (ValueError, KeyError) as e:
        print(f"Error al convertir 'hwindows': {e}. Asignando valor predeterminado 0.")
        hwindow = 1

    try:
        sigma = int(dic['sigma'])
    except (ValueError, KeyError) as e:
        print(f"Error al convertir 'sigma': {e}. Asignando valor predeterminado 0.")
        sigma = 3

    try:
        method = int(dic['method'])
    except (ValueError, KeyError) as e:
        print(f"Error al convertir 'method': {e}. Asignando valor predeterminado 0.")
        method = 0


    #for timet, freqarre in zip(timet1, freqarre1):
    for iteracion, (timet, freqarre) in enumerate(zip(timet1, freqarre1)):
        data_windows = create_data_windows(freqarre, windows, hwindow)
        ventanas_validas = detect_valid_windows(data_windows, method, sigma, Fs, LFFT, Tolerance, iteracion)
        resultado = filter_non_consecutive(ventanas_validas)

        #ventanas_freq, ventanas_time = extract_event_windows(timet, freqarre, resultado)

        #ventanas_freq1.append(ventanas_freq)
        #ventanas_time1.append(ventanas_time)

        # Visualización (comentado para entorno web)
        plot_results1(timet, freqarre, ventanas_validas, resultado, windows, hwindow, iteracion)

    return event_detect, ventanas_freq1, ventanas_time1

def deep_tuple(lst):
    return tuple(deep_tuple(x) if isinstance(x, list) else x for x in lst)

def min_max(data_windows, sigma, method, iteracion):
    amplitudes_maximas = np.ptp(data_windows, axis=1)    
    desviacion_estandar_m = np.mean(amplitudes_maximas)
    limite_superior_mm = desviacion_estandar_m + sigma * np.std(amplitudes_maximas)
    limite_inferior_mm = desviacion_estandar_m - sigma * np.std(amplitudes_maximas)
    ventanas_fuera_mm = np.where((amplitudes_maximas > limite_superior_mm) | (amplitudes_maximas < limite_inferior_mm))
    ventanas_fuera_mm = ventanas_fuera_mm[0]
    
    plot_results2(amplitudes_maximas,limite_inferior_mm,limite_superior_mm,method, iteracion)
    return ventanas_fuera_mm, amplitudes_maximas

def matrix_pencil(data_windows, sigma, Fs, Tolerance, method, iteracion):
    T = 1 / Fs
    # Aplicar Matrix Pencil sin bucles
    L = data_windows.shape[1] // 2
    def apply_matrix_pencil(signal):
        amplitudes, _, _, _, _ = matrix_pencil_analysis(signal, L, T, Tolerance)
        return np.max(amplitudes)
    peak_amplitudes_mp = np.array([apply_matrix_pencil(sig) for sig in data_windows])
    desviacion_estandar_mp_amplitudes = np.mean(peak_amplitudes_mp)
    peak_amplitudes_mp = np.array([apply_matrix_pencil(sig) for sig in data_windows])
    limite_superior_mp = desviacion_estandar_mp_amplitudes + sigma * np.std(peak_amplitudes_mp)
    limite_inferior_mp = desviacion_estandar_mp_amplitudes - sigma * np.std(peak_amplitudes_mp)
    ventanas_fuera_mp = np.where((peak_amplitudes_mp > limite_superior_mp) | (peak_amplitudes_mp < limite_inferior_mp))
    ventanas_fuera_mp = ventanas_fuera_mp[0]
    
    plot_results2(peak_amplitudes_mp,limite_inferior_mp,limite_superior_mp,method,iteracion)
    return ventanas_fuera_mp, peak_amplitudes_mp

def fft_events(data_windows, sigma, Fs, LFFT, method, iteracion):
    """
    Aplica FFT a las ventanas de datos y grafica correctamente la parte positiva del espectro.

    Parámetros:
    - data_windows: Matriz con ventanas de datos.
    - sigma: Umbral para detección de eventos.
    - Fs: Frecuencia de muestreo.
    - LFFT: Longitud de la FFT.
    - method: Método de análisis.
    - iteracion: Iteración actual del análisis.

    Retorna:
    - ventanas_fuera_fft: Índices de ventanas detectadas fuera del umbral.
    - amplitudes_maximasfft: Valores máximos de amplitud para cada ventana.
    """
    # Definir constantes
    T = 1 / Fs  # Período de muestreo

    # Calcular la FFT de todas las ventanas de datos
    FFTResult = np.fft.fft(data_windows, n=LFFT, axis=1)

    # Obtener la escala de frecuencias correctamente
    freqs = np.fft.fftfreq(LFFT, d=T)

    # Filtrar solo la parte positiva del espectro
    positive_freqs = freqs[:LFFT // 2]
    P2 = np.abs(FFTResult / LFFT)  # Normalización correcta
    P1 = P2[:, :LFFT // 2]  # Tomar solo la parte positiva

    # Obtener las amplitudes máximas en la parte positiva del espectro
    amplitudes_maximasfft = np.max(P1, axis=1)

    #Filtrar solo frecuencias entre 0 y 2 Hz
    mask = (positive_freqs >= 0) & (positive_freqs <= 2)
    freqs_filtradas = positive_freqs[mask]
    P1_filtrado = P1[:, mask]

    # Crear la figura de la FFT con un tamaño adecuado
    fig, ax = plt.subplots(figsize=(10, 5), dpi=100)

    # Definir colores dinámicos para cada ventana
    colormap = plt.cm.viridis
    num_series = P1_filtrado.shape[0]
    colors = [colormap(i / num_series) for i in range(num_series)]

    # Graficar la FFT para cada ventana
    for i in range(num_series):
        ax.plot(freqs_filtradas, P1_filtrado[i, :], label=f'Ventana {i+1}', color=colors[i], linewidth=1.5, alpha=0.8)

    # Anotar los 3 picos más altos dentro del rango
    num_peaks = min(3, num_series)
    for i in range(num_peaks):
        if len(freqs_filtradas) > 0:  # Asegurar que hay datos en el rango
            peak_index = np.argmax(P1_filtrado[i, :])  # Índice del pico máximo
            peak_freq = freqs_filtradas[peak_index]
            peak_amp = P1_filtrado[i, peak_index]

            ax.annotate(f'{peak_freq:.2f} Hz', 
                        xy=(peak_freq, peak_amp), 
                        xytext=(peak_freq, peak_amp * 1.1), 
                        arrowprops=dict(facecolor='red', arrowstyle='->'),
                        fontsize=10, color='black')

    # Configuración del gráfico
    ax.set_xlabel("Frecuencia (Hz)", fontsize=12)
    ax.set_ylabel("Amplitud |P1(f)|", fontsize=12)
    ax.set_title(f'Espectro FFT (0 - 2 Hz) - Iteración {iteracion}', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 2)  # 🔹 Fija los límites del eje X entre 0 y 2 Hz
    ax.grid(True, linestyle=':', linewidth=0.75, alpha=0.7)

    # Optimizar la leyenda (mostrar solo si hay ≤10 ventanas)
    if num_series <= 10:
        ax.legend(loc="upper right", fontsize=10, frameon=True)

    # Ajuste de márgenes
    plt.tight_layout()

    # Guardar imagen en formato base64
    buffer = io.BytesIO()
    fig.savefig(buffer, format='png', bbox_inches='tight', dpi=100)
    buffer.seek(0)
    image_png = buffer.getvalue()
    buffer.close()
    
    # Convertir imagen a base64 para mostrar en la web
    graphic = base64.b64encode(image_png).decode('utf-8')

    # Evitar acumulación de memoria
    plt.close(fig)

    # Guardar gráfico en la lista para visualización en la web
    graficos_generados.append((f'FFT (0-2 Hz) - Iteración {iteracion}', graphic))

    # Calcular límites para detección de eventos
    desviacion_estandar_amplitudes = np.mean(amplitudes_maximasfft)
    limite_superior_fft = desviacion_estandar_amplitudes + sigma * np.std(amplitudes_maximasfft)
    limite_inferior_fft = desviacion_estandar_amplitudes - sigma * np.std(amplitudes_maximasfft)

    # Identificar ventanas fuera de ±3σ
    ventanas_fuera_fft = np.where(
        (amplitudes_maximasfft > limite_superior_fft) | (amplitudes_maximasfft < limite_inferior_fft)
    )[0]

    # Generar gráfico de detección de eventos
    plot_results2(amplitudes_maximasfft, limite_inferior_fft, limite_superior_fft, method, iteracion)

    return ventanas_fuera_fft, amplitudes_maximasfft
def matrix_pencil_analysis(signal, L, dt, tol):
    """
    Apply the Matrix Pencil Method to analyze the signal.
    
    Args:
    - signal: The input signal (1D or 2D array). If 2D, each row is analyzed separately.
    - L: The pencil parameter.
    - dt: The sampling interval.
    - tol: Tolerance for rank determination.
    
    Returns:
    - amplitudes: Estimated amplitudes of the frequencies.
    - frequencies: Estimated frequencies of the components (Hz).
    - damping_ratios: Estimated damping ratios (dimensionless).
    - phases: Estimated initial phases (radians).0
    - y_hat: Reconstructed signal using the estimated components.
    """
    
    # Si la señal es 2D, procesamos cada fila de forma independiente
    if signal.ndim == 2:
        results = [matrix_pencil_analysis(row, L, dt, tol) for row in signal]
        return results  # Devuelve una lista de resultados para cada fila
    
    N = len(signal)
    
    if N <= L:
        raise ValueError("L must be smaller than the signal length N.")
    
    # Construcción de la matriz Y
    Y = np.column_stack([signal[ii:N - L + ii] for ii in range(L + 1)])
    
    # Descomposición en valores singulares (SVD)
    U, S, Vt = svd(Y, full_matrices=False)
    
    # Determinación del rango usando la tolerancia
    n = np.argmax(S[1:] / S[0] <= tol) if np.any(S[1:] / S[0] <= tol) else len(S) - 1

    # Seleccionar componentes significativos
    S = S[:n + 1]
    V_prime = Vt.T[:, :n + 1]

    # Construcción de las matrices V1 y V2
    V1_prime, V2_prime = V_prime[:-1, :], V_prime[1:, :]

    # Resolver sistema usando pseudoinversa
    A = pinv(V1_prime) @ V2_prime
    z = eig(A, right=False)  # Valores propios

    # Matriz Z con exponenciales
    Z = np.array([z**qt for qt in range(N)], dtype=complex)

    # Manejar NaNs e infinitos
    Z = np.nan_to_num(Z, nan=0.0, posinf=0.0, neginf=0.0)
    
    # Resolver sistema para obtener coeficientes R
    R = pinv(Z) @ signal

    # Cálculo de amplitudes
    amplitudes = np.abs(R)
    
    # Cálculo de frecuencias (Hz)
    frequencies = np.angle(z) / (2 * np.pi * dt)

    # Cálculo de damping ratios (\(\zeta\))
    damping_ratios = -np.log(np.abs(z)) / dt
    
    # Cálculo de fases iniciales
    phases = np.angle(R)

    # Ordenar los resultados por amplitud
    sorted_indices = np.argsort(amplitudes)[::-1]
    amplitudes, frequencies, damping_ratios, phases = (
        amplitudes[sorted_indices], frequencies[sorted_indices],
        damping_ratios[sorted_indices], phases[sorted_indices]
    )

    # Reconstrucción de la señal
    y_hat = Z @ R

    return amplitudes, frequencies, damping_ratios, phases, y_hat


def analyze(media_path,dic):

	v_time, v_freq = extract_time_and_frequency(media_path,dic)
	event, señales, tiempo = voltage_angle_difference(v_time, v_freq, dic)

	return graficos_generados

def get_signals_name(media_path, dic):

	csv_file = open(os.path.join(media_path, dic['file_name']))
	csv_data = csv.reader(csv_file)


	for list_signal in csv_data: break

	signals_name = ''
	for i in list_signal[1:]:
		signals_name = signals_name + ',' + i

	return signals_name[1:]

	







	