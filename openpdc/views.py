import json
from datetime import datetime, timedelta
from pathlib import Path

import pandas as pd
import requests
from dateutil import parser
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render
from django.template import loader

from .conf import *

# Create your views here

def get_data(request):
    #final = datetime.utcnow()
    #e = final.strftime("%m-%d-%Y %H:%M:%S")
    #if start_date is None:
    #    start_date=end_date-timedelta(days=3)
    #    start_date.strftime("%m-%d-%Y %H:%M:%S")
    #idx = ",".join([str(k) for k in pmu_list])
    #r = session.get(hist_url.format(idx, start_date, end_date))
    #data = r.content
    #print("ZZZZZZZZZZZZZZZZ", hist_url.format(idx, start_date, end_date))
    #data = json.loads(r.content).get("TimeSeriesDataPoints", [])
    #res = {}
    #for record in data:
    #    t = record["Time"]
    #    if t not in res:
    #        res[t] = {f"HID_{k}": np.NaN for k in UN.keys()}
    #        res[t]["Time"] = t
    #    res[t]["HID_{}".format(record["HistorianID"])] = record["Value"]
    #df = pd.DataFrame(res.values())
    #print(df.head())
    #df["Time"] = pd.to_datetime(df["Time"])
    context={"pmus":UN.items()}
    if request.method=='POST':
        context['selected_pmus'] = request.POST.getlist('pmus')
        date=parser.parse(request.POST.get('from'))
        context['start_date'] = date.strftime("%m-%d-%Y %H:%M:%S")
        context['end_date']=(date+timedelta(days=1)).strftime("%m-%d-%Y %H:%M:%S")
        idx = ",".join([str(k) for k in  context['selected_pmus']])
        r = session.get(hist_url.format(idx, context['start_date'],  context['end_date']))
        data = json.loads(r.content).get("TimeSeriesDataPoints", [])
        res = {}
        print("Qurey end")
        for record in data:
            t = record["Time"]
            if t not in res:
                res[t] = {UN[k]: np.NaN for k in UN.keys()}
                res[t]["Time"] = t
            res[t][UN[(record["HistorianID"])]] = record["Value"]
        df = pd.DataFrame(res.values())
        context["df"]=df
        df_t = df.T.reset_index()
        df_t.columns = ['Node_name'] + [i for i in range(len(df_t.columns)-1)]
        
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = 'attachment; filename="datos.csv"'
        df_t.to_csv(path_or_buf=response, index=False)
        return response 
    return render(request, 'openpdc/download.html',context)
