import numpy as np
#def Prony(self, Y, N, dt, fct_method):
def Prony(Y, N, dt, fct_method):
    modes = int(fct_method)
    if (modes == 0) or (modes > int(N / 2)):
        modes = int(N / 2)
        #self.val_W1ent7.set(modes)
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

    mag = 2.0 * abs(B)
    ang = np.angle(B)
    damp = lam.real
    freq = lam.imag / (2.0 * np.pi)

    freq_zrs = list(np.where(freq == 0.0)[0])
    for fz in freq_zrs:
        freq[fz] = -0.0001

    omga = 2.0 * np.pi * freq
    damprat = (damp / omga) * 100.0
    enrgy = (1.0 / 2.0) * (omga**2)  * (mag**2) * 100.0
    roots = z

    return modes, mag, ang, damp, freq, damprat, enrgy, roots
