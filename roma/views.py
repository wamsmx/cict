from django.shortcuts import render
from .ringdown import Prony, ERA, Matrix_Pencil as MP
import os
from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
import pandas as pd
import plotly.express as px
from pathlib import Path
import numpy as np
import plotly.graph_objects as go

BASE_DIR = Path(__file__).resolve().parent.parent
# Create your views here.

def pronyView(request):
    #return HttpResponse("Hola Mundo")
    algorithms={0:'Prony', 1:'Matrix Pencil', 3:'ERA'}
    data=pd.DataFrame()
    context={'algorithms':algorithms,
             'alg':1,}
    fn=f"data/Sample_Multi_methods.csv"
    fn=os.path.join(BASE_DIR,fn)
    alg=1
    if request.method=='POST':
            try:
                fn=request.FILES['filename']
            except:
                fn=f"data/Sample_Multi_methods.csv"
                fn=os.path.join(BASE_DIR,fn)
            alg=int(request.POST.getlist('algorithm')[0])
    ts=float(request.POST.get('ti',1))
    tf=float(request.POST.get('tf',5))
    data=pd.read_csv(fn, sep=',')
    if 'lat' in data.columns and  'lon' in data.columns :
        data=data.drop(columns=['lat','lon'])
    if alg==0:
        hp=int(request.POST.get('order',10))
        res,(Ti,Tf)=Prony(data, hp, ts, tf)
        context['order']=hp
    else:
        hp=float(request.POST.get('threshold',0.7))
        if alg==1:
            res,(Ti,Tf)=MP(data, hp, ts, tf)
        else:
            res,(Ti,Tf)=ERA(data, hp, ts, tf)
        context['threshold']=hp
    context['alg']=alg
    context['ts']=ts
    context['tf']=tf
    
    df=data.set_index('Node_name')
    fig = go.Figure()
    for idx in data.index[1:]:
        fig.add_trace(go.Scatter(x=data.loc[data.index[0]][Ti:Tf], y=data.loc[idx][Ti:Tf],mode='lines',name=idx))
    fig.update_layout(xaxis_title='Time (s)')
    context['fig']=fig.to_html()
    table=pd.DataFrame(res)
    table.sort_values(by='Frequency', inplace=True)
    context['tab']={'table':table.to_html(index=False),'latex':table.to_latex(index=False),'csv':table.to_csv(index=False)}
        #except Exception as e:
        #    print(">>>>>>>>>>>>", e)
            #pass
    return render(request, 'roma/prony.html',context)

