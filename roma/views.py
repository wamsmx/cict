from django.shortcuts import render
from .metodos import Prony as prony
import os
from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
import pandas as pd
import plotly.express as px
from pathlib import Path
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
# Create your views here.


def pronyView(request):
    #return HttpResponse("Hola Mundo")
    algorithms={0:'Prony', 1:'Matrix Pencil', 3:'ERA'}
    data=pd.DataFrame()
    context={'algorithms':algorithms,
             'texto':"Hola mundo",
             'alg':0,}
    if request.method=='POST':
        try:
            fn=request.FILES['filename']
            alg=int(request.POST.getlist('algorithm')[0])
            data=pd.read_csv(fn, sep=',')
            data=data.drop(columns=['lat','lon'])
            data=data.set_index('Node_name')
            y,t=data.loc['601'],data.loc['time'].values
            fig=px.line(y=y,x=t)
            context['fig']=fig.to_html()
            dt=t[1]-t[0]
            N=len(y)
            if alg==0:
                modes,mag,ang,damp,freq,damprat,enrgy,roots=prony(y,N,dt,5)
            elif alg==1:
                #Pencil
                modes,mag,ang,damp,freq,damprat,enrgy,roots=np.zeros(8)
            else:
                #ERA
                modes,mag,ang,damp,freq,damprat,enrgy,roots=np.zeros(8)
            context['alg']=alg    
            context['parameters']={"Modes":modes,"Amplitude":mag,
                                 "Phase":ang,"Damping":damp,
                                 "Frequency":freq,"Damping Ratio":damprat,
                                 "Energy":enrgy,"Poles":roots}
        except Exception as e:
            print(">>>>>>>>>>>>", e)
            pass
    return render(request, 'roma/prony.html',context)
    #fn=f"data/kundur.csv"
    #fn=os.path.join(BASE_DIR,fn)
    #if request.method=='POST':
    #    try: 
    #        fn=request.FILES['filename']
    #        event=fn.name
    #    except:
    #        pass
  
    #return render(request, 'af/clusters.html',{'tabs':tabs})
