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
             'texto':"Hola mundo",
             'alg':0,}
    if request.method=='POST':
        #try:
            fn=request.FILES['filename']
            alg=int(request.POST.getlist('algorithm')[0])
            ts=int(request.POST.get('ti',0))
            tf=int(request.POST.get('tf'))
            data=pd.read_csv(fn, sep=',')
            if 'lat' in data.columns and  'lon' in data.columns :
                data=data.drop(columns=['lat','lon'])
            if alg==0:
                hp=int(request.POST.get('order',3))
                res=Prony(data, ts, tf, hp)
            else:
                hp=float(request.POST.get('threshold',0.5))
                print(ts,tf,hp)
                if alg==1:
                    res=MP(data, ts, tf, hp)
                else:
                    res=ERA(data, ts, tf, hp)
            context['alg']=alg
            df=data.set_index('Node_name')
            fig = go.Figure()
            for idx in data.index[1:]:
                fig.add_trace(go.Scatter(x=data.loc[data.index[0]][ts:tf], y=data.loc[idx][ts:tf],mode='lines',name=idx))
            context['fig']=fig.to_html()
            table=pd.DataFrame(res)
            context['tab']={'table':table.to_html(index=False),'latex':table.to_latex(index=False),'csv':table.to_csv(index=False)}
        #except Exception as e:
        #    print(">>>>>>>>>>>>", e)
            #pass
    return render(request, 'roma/prony.html',context)

