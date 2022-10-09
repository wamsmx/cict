
import os
import csv
import numpy 				as np
import matplotlib.pyplot 	as plt

from scipy.linalg 			import hankel


def get_data(media_path, dic):

	csv_file = open(os.path.join(media_path, dic['file_name']))
	csv_data = csv.reader(csv_file)
	print(dic)
	csv_list = []
	for row in csv_data:
		csv_list.append(row)

	if dic['csv_layout'] == 'vcsv':

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

	if dic['csv_layout'] == 'hcsv':

		h_vec   = []
		y_vec   = []
		t, y    = [], []

		for index, icsv_list in enumerate(csv_list):
            
			h_vec.append(icsv_list[0])

			if index == 0:
				t = icsv_list[1:]
				t = list(map(float, t))
			else:
				iy = icsv_list[1:]
				iy = list(map(float, iy))
				y.append(iy)
				
		h_vec = h_vec[1:]
		t_vec = np.array(t)
		t_vec = t_vec - t_vec[0]
		y_vec = np.array(y)
		y_vec = y_vec.T

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


def analyze(media_path, dic):

	h_vec, t_vec, y_vec = get_data(media_path, dic)

	return ringdown(h_vec, t_vec, y_vec, dic)


def get_signals_name(media_path, dic):

	csv_file = open(os.path.join(media_path, dic['file_name']))
	csv_data = csv.reader(csv_file)

	if dic['csv_layout'] == 'vcsv':

		for list_signal in csv_data: break

		signals_name = ''
		for i in list_signal[1:]:
		    signals_name = signals_name + ',' + i

		return signals_name[1:]

	if dic['csv_layout'] == 'hcsv':

		signals_name = ''
		for index, icsv_list in enumerate(csv_data):
			signals_name = signals_name + ',' + icsv_list[0]

		return signals_name[1:]







	
