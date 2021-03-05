import pandas as pd
import numpy as np
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import pdist
#from affinity import Affinity
import plotly.graph_objects as go
import plotly.express as px

def clusterGenerators(data,id_machines):
    gen_data=data[id_machines,:]
    non_gen_data=data[np.max(id_machines)+1:,:]
    pref=np.min(-pdist(data,'correlation'))
    S=-1*pairwise_distances(gen_data,metric='correlation')
    for ii in range(len(gen_data)): S[ii,ii]=pref
    aclf=AffinityPropagation(affinity='precomputed',random_state=None).fit(S)
    nearest=[]
    if len(non_gen_data):
        nearest=np.argmin(pairwise_distances(non_gen_data,gen_data,metric='correlation'), 
                          axis=1)
    labels=[]
    nl=[aclf.cluster_centers_indices_[l] for l in aclf.labels_]
    #print(nl)
    #print(nearest)
    for nn in nearest:
        labels.append(nl[nn])
    return aclf, nl+labels


def subRanges(l):
        if  len(l)==0:
            return ""
        i,k,f,labels=l[0],l[0],l[0],[]
        l=list(l)
        l.sort
        if len(l)==1:
            return str(l[0])
        for j in l[1:]:
            if f==j-1:
                f=j
            else:
                #print(i,f,j)
                labels.append("%(i)01d" %locals())
                #labels.append("N%(i)02d" %locals())
                if i!=f:
                    labels[-1]=("%(i)01d-%(f)01d" %locals())
                i=j
                f=j
        labels.append("%(i)01d" %locals())
        if i!=f:
            labels[-1]=("%(i)01d-%(f)01d" %locals())
        #return "[%s]" %",".join(labels)
        return ",".join(labels)

def labels_to_table(labels,fnames):#data,file):
    #labels=data[file]
    #print(len(labels))
    L=set(labels)
    df=[]
    for l in L:
        s=np.where(labels==l)[0]
        #print(s)
        g,ng=[],[]
        for i in s:
            g.append(fnames[i])
        #df.append([','.join(g),','.join(ng),fnames[l]])
        df.append([','.join(g),fnames[l]])
    #df=pd.DataFrame(df,columns=['Coherent Generators','Associated Non-Generator Buses','Centroid'])
    df=pd.DataFrame(df,columns=['Sources (buses)','Centroid'])
    df['Area id']=[i+1 for i in range(len(df))]
    #df=df.set_index('Area')
    return df[['Area id','Sources (buses)','Centroid']] 

def plotMap(fn,gd):
    fnames=[x.replace('source','') for x in pd.read_csv(fn, sep=',').columns[1:-1]]
    fnet01=pd.read_csv(fn, sep=',', skiprows=[0], header=None)
    gps=pd.read_csv(gd)
    data=fnet01.values[:,1:-4]
    clf,labels=clusterGenerators(data.transpose(),[i for i in range(len(fnames))])
    ddf=[]
    for i,l in enumerate(labels):
        idnode=int(fnames[i])
        #print(idnode)
        #print(gps.FDRid)
        if idnode not in gps.FDRid.values:
            print(idnode)
            continue
        lon, lat=gps[gps.FDRid==idnode][['Longitude','Latitude']].values[0]
        ddf.append({'node':fnames[i], 'cluster_id':fnames[l], 'lat':lat, 'lon':lon})
    ddf=pd.DataFrame(ddf)
    fig = px.scatter_geo(ddf,lat=ddf.lat, lon=ddf.lon,color=ddf.cluster_id,hover_name="node")
    fig.update_layout(geo=dict(lataxis={'range':(23,55)},lonaxis={'range':(-110,-62)}))
    table=labels_to_table(labels,fnames)
    return fig,table
    
import glob
if __name__=="__main__":
    #fn="data/Gentrip/data_20200716_160955.csv"
    gd='data/units_gps.csv'
    for fn in glob.glob("data/Gentrip/*.csv"):
        #fn="data/Gentrip/data_20200908_132512.csv"
        f,table=plotMap(fn,gd)
    #print(table.to_html(index=None))
