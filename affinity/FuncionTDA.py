import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio 
from .finalFn_iterative_cleaned import finalFn_iterative_cleaned 
from statistics import mean

import plotly.graph_objects as go
import plotly.express as px

#Funtions------------------------------------
def movmedian(a,ki):
    import math
    from statistics import median
    t=math.trunc(ki/2)
    q=round(ki/2+0.00000001)
    k=np.array([])
    resultado=np.zeros(a.shape[0])

    for i in range(a.shape[0]):
        if(q<ki):
            k=a[:q]
            resultado[i]=median(k)
            q=q+1 
        elif(i==a.shape[0]-t):
            k=a[i-t:]
            resultado[i]=median(k)
        else:
            k=a[i-t:i+ki-t]
            resultado[i]=median(k)
    return resultado

def create_stem_data(x,y, baseline=0.):
    '''makes y data passing 0 before inbetween actual value to create data for a stem plot
    x,y are 3 times the original length
    '''
    x=np.repeat(x,3)
    y=np.repeat(y,3)
    y[::3]=y[2::3]=baseline
    return x,y


def plot_map(labels, data):
    colors_html={0:'#1f77b4', 1:'#ff7f0e', 2:'#2ca02c', 3:'#d62728', 4:'#9467bd', 5:'#8c564b', 
             6:'#e377c2', 7:'#7f7f7f', 8:'#bcbd22', 9:'#17becf',10:'#c65758',11:'#777fff'} 
    fig = go.Figure() # Creamos una figura
    k=0
    for l in set(labels): # iteramos para cada area
        df=data[labels==l] # seleccionamos los datos correspondientes a la etiqueta l
        if(k>11):
            k=0
        fig.add_trace(go.Scattermapbox( mode='markers', lon = df.lon, lat = df.lat, text=df.Node_name,
                                       hoverinfo='text', name=f"Area {int(l+1)}", 
                                       marker = { 'size': 8, 'color':colors_html[k]}))
        k=k+1
    fig.update_layout( mapbox = { 'style': "open-street-map", 'zoom':2.5,
                                 'center': {'lon': df.lon.mean(), 'lat': df.lat.mean()}},
                      showlegend = True, height=600, width=750)        

    return fig

def scalematrix(m, scale=10):
  # Crear matriz con ceros del tamaño apropiado
  r = np.zeros(((m.shape[0]-1)*scale+1, (m.shape[1]-1)*scale+1))

  # Rellenar filas multiplo de scale (interpolando entre valores de los elementos de la fila)
  for fil in range(m.shape[0]):
    for col in range(m.shape[1]-1):
      r[fil*scale, col*scale:(col+1)*scale+1] = np.linspace(m[fil,col], m[fil,col+1], scale+1)

  # Rellenar resto de ceros, interpolando entre elementos de las columnas
  for fil in range(m.shape[0]-1):
    for col in range(r.shape[1]):
      r[fil*scale:(fil+1)*scale + 1, col] = np.linspace(r[fil*scale,col], r[(fil+1)*scale, col], scale+1)

  
  xr=np.arange(0,r.shape[0]/10,0.1)
  yr=xr
  return r,xr,yr

#Main----------------------------------

def FuntionTDA(freq, dt,t0,tf,cte):
    print("zzzzzzz",freq)
    freq=freq.to_numpy()
    print("yyyyyyy",freq)
    #freq=freq[0:,1:-2]
    
    freq=np.transpose(freq)

    #freq=data
    origin=freq;

    x=freq-60;
    x=x[:,0:]
    x=np.transpose(x)
    for i in range(x.shape[0]):
        x[i]=movmedian(x[i],5)
    x=np.transpose(x)

    #gen_buses=np.vstack([np.arange(1,ngen+1,1),np.zeros(ngen)])

    (test_time,groups,centroid,typicality,clusterMean,rho)=finalFn_iterative_cleaned(x,dt,t0,tf,cte)

    centroid=np.diag(groups[centroid.astype(int),:])


    x1,y1 = create_stem_data(np.arange(typicality.shape[0] ),np.transpose(typicality))
    fig1 = px.line(x=x1,y=y1)

    fig1.add_trace( go.Scatter(x=np.arange(typicality.shape[0]), y=np.transpose(typicality),  mode='markers', 
                                marker=dict(color='blue'), showlegend=False ))
    
    fig1.add_trace(go.Scatter(x=centroid,y=typicality[centroid.astype(int)], mode='markers', marker=dict(color='red'), showlegend=False))


    for i in range(centroid.shape[0]):
        txt='Area '+str(i+1)
        fig1.add_trace(go.Scatter(x=[centroid[i]+0.5,centroid[i]+0.5],y=[0.0,0.04], mode='lines',name=txt))

    fig1.update_layout(title='Tau(k) after ranking ', xaxis_title='ranked_τ (k)  ', yaxis_title="Ranked Typicality τ_k (ν_k) ")


    #-histogram---------------------------------------
    rho,xr,yr=scalematrix(rho)
    fig2=px.imshow(np.rot90(rho), 
                labels=dict(x="Bus k", y="Bus j", title='Correlation (ρ) between buses k and j', color='Correlation ρ between buses k and j'), 
                color_continuous_scale='jet',  origin='lower',x=xr,y=yr)


    #Generating figures

    plot_time=np.arange(0,(x.shape[0]*dt),dt)
    TFR_Groups=[]

    for k in range(groups.shape[1]):
        txt=np.array([])
        fig=go.Figure()
        for t in range(groups.shape[0]):
        #for t in range(2):
            if (groups[t,k]!=0):
                txt='Bus '+str(int(groups[t,k]))
                fig=fig.add_trace( go.Scatter(x=plot_time, y=x[:,int(groups[t,k])-1],  mode='lines', name=txt))
    
        txt='TFR of Group '+str(k+1)

        fig.update_layout(title=txt, xaxis_title='Time (s)', yaxis_title='Δf (Hz)')
        TFR_Groups.append(fig)
    
    
    compa=np.arange(1,int(np.max(groups)+1),1)
    labels=np.zeros(int(np.max(groups)))

    for j in range(groups.shape[1]):
        for i in range(groups.shape[0]):
            if groups[i,j]>0:
                labels[np.where(groups[i,j]==compa)[0][0]]=j

    
    #figMap=None
    #datapd[['lat','lon']].head()
    #figMap=plot_map(labels,datapd)

    return groups, fig1, fig2,TFR_Groups #, figMap
