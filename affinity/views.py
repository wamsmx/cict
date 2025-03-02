from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from .clustering import plotMap
import os
from pathlib import Path
from django.template import loader
import json
import pandas as pd

from bokeh.server.server import Server
from bokeh.application import Application
from bokeh.application.handlers.function import FunctionHandler
from bokeh.plotting import figure, ColumnDataSource
import requests, json
from requests_ntlm import HttpNtlmAuth as wauth
import pandas as pd
import plotly.express as px
from django.core import serializers


BASE_DIR = Path(__file__).resolve().parent.parent
# Create your views here.

import glob

def clusters(request, d='signals.csv'):
    #fn=f"data/Gentrip/data_{d}.csv"
    #gd='data/units_gps.csv'
    #gd=os.path.join(BASE_DIR,gd)
    tabs=[]
    dist='euclidean'
    distances={
        'euclidean': 'Euclidean', 'correlation':'Correlation',
        'jensenshannon':'Jensen Shannon','braycurtis':'Bray-Curtis',
        'chebyshev':'Chebyshev', 'cityblock':'Cityblock' #, 'dcor':'Distance Correlation' 
    }
    algorithms={'tda':'TDA','kmeans':'K-Means','hierarchical':'Hierarchical',
               'affinity':'Affinity Propagation',"vat":"VAT"}
    linkages={'centroid':'Centroid','single':'Single','complete':'Complete','average':'Average',
              'weighted':'Weighted','median':'Median','ward':'Ward'}
    n_clusters={str(i):i for i in range(2,21)}
    #events=['sample','20200716_160955','20200919_020549',
            #'20200909_133110','20200921_122833','20200912_062305',
            #'20200921_153352','20200908_132512','20200918_134511',
            #'20200909_115746']
    #event=events[0]
    event='Kundur sample'
    fn=f"data/kundur.csv"
    fn=os.path.join(BASE_DIR,fn)
    alg='tda'
    #if request.method=='GET':
    #for i,event in enumerate(events):
    nc,lk=2,'centroid'
    if request.method=='POST':
        try: 
            fn=request.FILES['filename']
            event=fn.name
        except:
            pass
        print(request.POST)
        dist=request.POST.getlist('distance')[0]
        alg=request.POST.getlist('algorithm')[0]
        nc=int(request.POST.getlist('nclusters')[0])
        lk=request.POST.getlist('linkage')[0]
    fig,table,groups,mapa,hmap,labels_,tda=plotMap(fn,alg=alg,dist=dist,n_clusters=nc,link=lk)
    i=0
    with pd.option_context("max_colwidth", 1000):
        latex_table= table.to_latex(index=False)
        csv_table= table.to_csv(index=False)
    tabs.append({'id':event, 'name':event, 'fig':fig.to_html(),'curves':groups,'csv':csv_table,
                     'table':table.to_html(index=None),'active': i==0, 'latex':latex_table})
        #return HttpResponse(pd.read_csv(fn, sep='').to_html())
    return render(request, 'af/clusters.html',{'tabs':tabs,'dist':dist,'alg':alg,'lk':lk,'nc':nc,
                                               'algorithms':algorithms,'nclusters':n_clusters, 'mapa':mapa,
                                               'distances':distances,'linkages':linkages, 'hmap':hmap,'tda':tda})


def mapmx(request):
    return render(request, 'index.html')

def signals(request):
    return render(request, 'signals.html')

def Home (request):
    return render(request, 'mapdex.html')

def Zonemap(request):
    return render(request, 'zonedex.html')



def frequency(request):
    session=requests.Session()
    session.auth=wauth('potencia2021','W@ms_project202x')
    url="http://148.216.38.78:6152/historian/timeseriesdata/read/current/646,7,59,201,81,54,49,631,806,796,286,786/json"
    r=session.get(url)
    source=(r.content)
    source=json.loads(r.content)
    data2=json.dumps(source)
    return JsonResponse({'data':data2})

def angle(request):
    session=requests.Session()
    session.auth=wauth('potencia2021','W@ms_project202x')
    url="http://localhost:6152/historian/timeseriesdata/read/current/650,11,63,205,85,58,53,635,260,250,290/json"
    r=session.get(url)
    source=(r.content)
    source=json.loads(r.content)
    data=json.dumps(source)
    return JsonResponse({'data':data})
   


