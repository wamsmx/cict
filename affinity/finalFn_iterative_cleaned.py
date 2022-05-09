import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import math
from .step1_forPwrSys import step1_forPwrSys 
from .step1 import step1
from .step4 import step4
from .step5 import step5
from time import time

from plotly.offline import plot

#finalFn_iteractive_cleaned

def dsearchn(x,xi):
    t=np.zeros([xi.shape[0]][0])
    d=np.zeros([xi.shape[0]][0])
    import numpy.matlib
    for i in range(xi.shape[0]):
        yi=xi[i,:]
        t[i]=np.min(np.sum(pow(x-yi,2),axis=0))
        d[i]=np.sqrt(np.min(np.sum(pow(x-yi,2),axis=0)))
    return d, t

def encontrarPontoMaisProximo(ponto,colecao):
    distancias,_=dsearchn(ponto,colecao)
    idx=np.argmin(distancias)
    return idx

def  procuraMaximosLocais(vector):

    isMaximum=np.zeros(vector.shape)

    vectorMinus1=vector[0:-3][:]
    vectorCompar=vector[1:-2][:]
    vectorPlus1=vector[3:-1][:]

    parcelaMis=vectorCompar>vectorMinus1
    parcelaPlus=vectorCompar>vectorPlus1
    return parcelaPlus



#-------------------------------------------------------------------------------

def finalFn_iterative_cleaned(data,dt,t0,tf,cte):
    print("$$$$$$$$$$",data)
    if tf>(data.shape)[0]:
        x=data[(t0):]
    else:
        x=data[t0:(tf),:]
    
    dima=(np.shape(x)[0])
    dimb=(np.shape(x)[1])
    time=np.transpose((np.linspace(0,(dima*dt),dima)))

    ED=np.zeros([dimb,dimb])
    LID=np.zeros([dimb,dimb])
    RHO=np.zeros([dimb,dimb])
    for k in range(dimb):
        for t in range(dimb):
            ED[k][t]=np.sqrt(sum(pow((x[:,k]-x[:,t]),2)))
            LID[k][t]=sum(x[:,k]-x[:,t])
            RHO[k][t],_=pearsonr(x[:,k],x[:,t])
    
    import time
    tStart=time.time()

    aux=ED
    x=aux;
    datasetX=aux;

    (resultado_step1, rho, sigmaN, epsN)=step1_forPwrSys(datasetX)

    tauD=resultado_step1["tauD"]
    tauD=abs(tauD)
    idxMax=np.argmax(tauD) #comprobar
    u=x;#unique values

    uStar=np.zeros([2,dimb])
    uStar[-1]=u[idxMax]
    tauStar=tauD[idxMax];

    u=np.delete(u,idxMax,0)
    tauD=np.delete(tauD,idxMax)


    #Step V from algoritm 1.

    while not((np.prod(u.shape)==0)):
        idxMax=encontrarPontoMaisProximo(uStar[-1],u)
        uStar=np.vstack([uStar,u[idxMax][:]])
        tauStar=np.vstack([tauStar,tauD[idxMax]])
        u=np.delete(u,idxMax,0)
        tauD=np.delete(tauD,idxMax)
    uStar=np.delete(uStar, 0, axis=0)


    #finding peaks
    from scipy.signal import find_peaks
    maxIdx,maxima = find_peaks(tauStar[:,0], height=0)
    maxima=(np.array(list(maxima.items()), dtype=type))[0,1]

    if (np.prod(maxima.shape)==0):
        maxima=tauStar[0]
        auxy=np.argwhere(tauStar<tauStar[0])
        auxy=auxy[-1]
        maxima=np.append(maxima,[tauStar(auxy)])
        maxIdx=np.append(1,[maxIdx])

    maxima=np.append(max(tauStar),[maxima])
    maxIdx=np.append(0, [maxIdx])


    uStarStar=uStar[maxIdx]
    tauStarStar=tauStar[maxIdx]

    #Important step

    for q in range(5):
        if (q>0):
            uStarStar=resultado_step4["ClusterMean"][np.logical_not(resultado_step5["muR"])] 
            resultado_step4=step4(uStarStar,uStar)
        else :
            resultado_step4=step4(uStarStar,uStar)
        #resultado_step5=step5(np.round(resultado_step4["ClusterMean"],4),np.round(resultado_step4["stdDev"],4),np.round(resultado_step4["tauMean"],4),cte)
        resultado_step5=step5(resultado_step4["ClusterMean"],resultado_step4["stdDev"],resultado_step4["tauMean"],cte)
    
    #end of for 10

    tauData=resultado_step4["tauData"]
    tauMean=resultado_step4["tauMean"]


    import time
    test_time=time.time()-tStart
    bus=np.zeros([datasetX.shape[0],2])


    #Fiding the bus corresponding to the datapoint
    for w in range(datasetX.shape[0]):
        for p in range(uStar.shape[0]):
            if (sum(datasetX[w]==uStar[p])>0):
                if (bus[p,0]==0):
                    bus[p,0]=w;
                    bus[p,1]=resultado_step4["Cluster"][p]
                pass

    _,a=np.unique(bus[:,1],return_counts=True)
    groups=np.zeros([np.max(a),int(np.amax(bus[:,1])+1)])
    aa=np.zeros(groups.shape[1])

    for i in range(int(np.amax(bus[:,1])+1)):
        a2=bus[bus[:,1]==i,0]+1
        aa[i]=a2.shape[0]
        groups[:a2.shape[0],i]=a2

    i=0
    while i<10000:
        if (i>=groups.shape[1]):
            break
        else:
            if(np.any(groups[:,i])==0):
                groups=np.delete(groups,0,i)
            pass
        i=i+1

    busIC=np.zeros(groups.shape[1])

    for j in range(groups.shape[1]):
        busIC[j]=encontrarPontoMaisProximo((resultado_step4['ClusterMean'][j]),datasetX[:int(aa[j]):])


    return test_time, groups, busIC, tauData, tauMean, rho

