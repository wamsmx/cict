import numpy as np 
from scipy.stats import pearsonr
import math
from .step1 import step1

#step4

def calculaVencedor(ponto,colecao):
    diferenca=colecao-ponto
    saida=np.sum(pow(diferenca,2),axis=1)
    return np.argmin(saida)

def calculaDesvioPadrao(dados):
	media=np.mean(dados,axis=0)
	diferenca=dados-media
	normas=np.zeros(diferenca.shape[0])
	for i in range(diferenca.shape[0]):
		normas[i]=np.linalg.norm(diferenca[i][:])
	std=np.mean(normas,axis=0)
	return std

def step4(uStarStar,x):
	idCluster=np.zeros(x.shape[0])
	for i in range(x.shape[0]):
		idCluster[i]=calculaVencedor(x[i],uStarStar)
	output={ 
	"x" : x,
	"Cluster" : idCluster,
	"idCluster" : np.sort(np.unique(idCluster)),
	}
	output["ClusterMean"] =np.zeros([output["idCluster"].shape[0],uStarStar.shape[1]])
	output["stdDev"]= np.zeros(output["idCluster"].shape)
	for i in range(output["idCluster"].shape[0]):
	    cluster=output["idCluster"][i]
	    if (x[idCluster==cluster][:].shape[0]>1):
	    	output["ClusterMean"][i,:]=np.mean(x[idCluster==cluster][:],axis=0)
	    else  :
	    	output["ClusterMean"][i,:]=np.mean(x[idCluster==cluster][:])
	    output["stdDev"][i]=calculaDesvioPadrao(x[idCluster==cluster][:])

	output["nCluster"]=output["idCluster"].shape[0]
	#Calculating tauD for everyCluster
	output["tauMean"]=np.zeros(output["nCluster"])
	output["tauData"]=np.zeros(x.shape[0])
	for i in range(output["idCluster"].shape[0]):
		cluster=output["idCluster"][i]
		clusterDataSet=output["x"][output["Cluster"]==cluster]
		clusterDataSet=np.vstack([clusterDataSet,output["ClusterMean"][i][:]+1e-9])
		saida=step1(clusterDataSet)
		output["tauData"][output["Cluster"]==cluster]=saida["tauD"][0:-1]
		output["tauMean"][i]=saida["tauD"][-1]
	return output

