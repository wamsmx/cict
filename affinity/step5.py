
import numpy as np
import scipy.io as sio

def step5(mean1,sigma,tauMean,cte):
	n=mean1.shape[0]
	muR=np.zeros(mean1.shape[0])

	for i in range(n):
		for j in range(n):
			if(i!=j):
				mu_i=mean1[i][:]
				mu_j=mean1[j][:]
				sigma_i=sigma[i]
				tau_i=tauMean[i]
				tau_j=tauMean[j]
				condicao1=np.linalg.norm((mu_i-mu_j))<=(cte*sigma_i)
				condicao2=tau_i<tau_j
				if (condicao1&condicao2):
					muR[i]=1
	output={"muR":muR}
	return output
