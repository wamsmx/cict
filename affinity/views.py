from django.shortcuts import render
from django.http import HttpResponse
from .clustering import plotMap
import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent
# Create your views here.

import glob
def clusters(request, d='20200716_035506'):
    #fn=f"data/Gentrip/data_{d}.csv"
    events=['20200716_035506','20200716_160955','20200919_020549',
            '20200909_133110','20200921_122833','20200912_062305',
            '20200921_153352','20200908_132512','20200918_134511',
            '20200909_115746']
    tabs=[]
    gd='data/units_gps.csv'
    gd=os.path.join(BASE_DIR,gd)
    for i,event in enumerate(events):
        fn=f"data/Gentrip/data_{event}.csv"
        fn=os.path.join(BASE_DIR,fn)
        fig,table=plotMap(fn,gd)
        tabs.append({'id':event, 'name':event, 'fig':fig.to_html(),
                     'table':table.to_html(index=None),'active': i==0})
    return render(request, 'af/clusters.html',{'tabs':tabs})
