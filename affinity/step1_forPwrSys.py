import numpy as np
from scipy.stats import pearsonr

#step1_forPwrSys

def calculaN(xi,x):
	diferenca=x-xi
	saida=1/(sum(np.sqrt(np.sum(pow(diferenca,2),axis=1))))
	return saida

def calculaqN(xi,x):
	xi=(np.transpose(xi))
	diferenca=np.zeros(x.shape[0])
	for i in range(x.shape[0]):
		diferenca[i],_=pearsonr(x[:,i],xi);

	saida=sum((pow(diferenca,2)));
	return (saida,diferenca)


#step1_forPwrSys
def step1_forPwrSys(datasetX):
	cN=np.zeros(datasetX.shape[0])
	qN=np.zeros(datasetX.shape[0])
	rho=np.zeros([datasetX.shape[0],datasetX.shape[0]])
	print(">>>>>>>>",datasetX)
	for i in range(datasetX.shape[0]):
		cN[i]=calculaN(datasetX[i][:],datasetX)
		(qN[i],rho[:][i])=calculaqN(datasetX[i][:],datasetX);
	
	sigmaN=2*qN/sum(qN)

	epsN=(datasetX.shape[0])*sigmaN;  #Normalized eccentrily

	DN=1/epsN
	tauN=DN/sum(DN)
	tauD=tauN
	output={"tauD":tauD}
	return (output,rho,sigmaN,epsN)
