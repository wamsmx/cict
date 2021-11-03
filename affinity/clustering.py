import pandas as pd
import numpy as np
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import pairwise_distances
from sklearn.metrics import adjusted_rand_score as ari
from sklearn.metrics import silhouette_score as sil
from sklearn.metrics import mean_squared_error as mse
from scipy.spatial.distance import pdist,squareform
from scipy.signal import resample
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree,fcluster
import plotly.graph_objects as go
#from affinity import Affinity
import plotly.graph_objects as go
import plotly.express as px
from dcor import distance_correlation as dc
from plotly.subplots import make_subplots
from matplotlib.colors import rgb2hex
from matplotlib import pyplot as plt
from sklearn.cluster import KMeans

cmap=plt.get_cmap('tab20')


def dcorr(X):
    n=len(X)
    d=[dc(X[i],X[j]) for i in range(n) for j in range(i+1,n)]
    return np.array(d)

#def groupsPlot(data, table, ddf):
#    n=len(set(ddf.cluster_id))+1
#    rows = n//3
#    re=n%3
#    specs=[]
#    fnames=data.index
#    if re==1:
#        specs=[[{"colspan":3},None,None]]
#        rows=rows+1
#    elif re==2:
#        specs=[[{"colspan":2},None,{}]]
#        rows=rows+1
#    specs+=[[{},{},{}] for i in range(len(specs),rows)]
#    fig=make_subplots(rows=rows, cols=3, specs=specs,shared_yaxes=True,
#                      shared_xaxes=True)
#    k=0
#    for i in range(0,rows):
#        for j in range(0,3):
#                if k==0:
#                    print("XXXXXXX",i,j,k)
#                    k=k+1
#                    for y in data.values:
#                        fig.add_trace(go.Scatter(y=y,mode='lines',showlegend = False),
#                                            #  line_color=rgb2hex(cm(i%20))),
#                                  row=i+1, col=j+1
#                        )        
#                elif specs[i][j]!=None:
#                    print("XXXXXXX",i,j,k)
#                    idx=ddf.cluster_id==table.Centroid.values[k-1]
#                    df=data[idx].values
#                    fn=fnames[idx]
#                    k=k+1
#                    for y,nm in zip(df,fn):
#                        fig.add_trace(go.Scatter(y=y,mode='lines', legendgroup =f'{k}', name=f'{nm}'),
#                                  row=i+1, col=j+1
#                        )
#                    print("XXXXXXX",i,j)
#    fig.update_layout(height=800, width=800, title_text="Stacked Subplots")
#    return fig

def groupsPlot(data,table,ddf):
    #fig=make_subplots(rows=2, cols=1,specs=[[{}],[{}]],shared_xaxes=True)
    fig=go.Figure()
    options=[0 for i in range(len(table.Center)+1)]
    all=[True for i in range(len(table.Center))]
    for y,nm in zip(data.values,data.index):
        fig.add_trace(go.Scatter(y=y,mode='lines', name=f'{nm}',visible=True,showlegend=False))
    options[0]=dict(label="All",method="update",args=[{"visible":all}])
    for i,cid in enumerate(table.Center.values):
        idx= ddf.cluster_id==cid
        options[i+1]=dict(label=f"Area {i+1}",method="update",args=[{"visible":idx.values}])
        #df=data[ddf.cluster_id==cid].values
        #fn=data.index[ddf.cluster_id==cid].values
    fig.update_layout(updatemenus=[dict(y=1.1,x=0.1,active=0,buttons=options)])
    fig.update_layout(height=600, width=800)
    return fig  
    

def clusterGenerators(data,id_machines=[],distance='correlation'):
    if not id_machines:
        id_machines=[i for i in range(len(data))]
    gen_data=data[id_machines,:]
    non_gen_data=data[np.max(id_machines)+1:,:]
    if distance!='dcor':
        S=-pdist(data,distance)
    else:
        S=-(1-dcorr(data))
    pref=np.median(S)
    S=squareform(S)
    #print("##########",S.shape)
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

def _centers(data,labels):
    clabels=np.array([l for l in labels])
    for l in set(labels):
        centroid=np.mean(data[labels==l], axis=0)
        dists=pairwise_distances(data.values,[centroid])
        clabels[labels==l]=np.argmin(dists)
    return clabels

def kmeans(data, n_clusters):
    kclf   = KMeans(n_clusters=n_clusters, random_state=33).fit(data.values)
    labels = kclf.predict(data.values)
    return kclf,_centers(data,labels)

def hierarchical(data, link='average'):
        Z = linkage(data.values, link)
        last = Z[:, 2]
        last_rev = last[::-1]
        #last=last[:-4]
        #acceleration = np.diff(last, 2)
        #acceleration_rev = acceleration[::-1]
        #print(acceleration_rev)
        #k = acceleration_rev.argmax() + 2
        #print("kkkkkkk ",k)
        n=len(last_rev)-1
        a=(last_rev[n]-last_rev[0])/(n-0)
        b=n*last_rev[0]/(n-0)
        D=np.zeros(n+1)
        for i in range(len(D)):
            D[i]=np.abs(-last_rev[i]+a*i+b)/np.sqrt(a**2+1)
        k=np.argmax(D)+1
        labels=fcluster(Z, k, criterion='maxclust')
        clabels=_centers(data,labels)
        return last_rev,clabels

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
        print(l,type(l))
        df.append([', '.join(g),fnames[l]])
    #df=pd.DataFrame(df,columns=['Coherent Generators','Associated Non-Generator Buses','Centroid'])
    df=pd.DataFrame(df,columns=['Sources (buses)','Center'])
    df['Area id']=[i+1 for i in range(len(df))]
    #df=df.set_index('Area')
    return df[['Area id','Sources (buses)','Center']] 

def plotMap(fn,alg='affinity', dist='correlation', n_clusters=2, link='ward'):
    #fnames=[x.replace('source','') for x in pd.read_csv(fn, sep=',').columns[1:-1]]
    #if type(fn)!=str:
    #    fn.seek(0)
    data=pd.read_csv(fn, sep=',',index_col=0)
    #gps=pd.read_csv(gd)
    ddf=[]
    if 'lat' in data.columns and 'lon' in data.columns:
       ddf=data[['lat','lon']]
       data=data.iloc[:,0:-2]
    if 'time' in data.index:
        data=data[1:]
        if len(ddf): 
            ddf=ddf[1:]
    if alg=='affinity':
        clf,labels=clusterGenerators(data.values,distance=dist)
    elif alg=='hierarchical':
        clf,labels=hierarchical(data,link)
    else:
        clf,labels=kmeans(data,n_clusters)
    fnames=data.index
    fig = go.Figure()
    table=labels_to_table(labels,fnames)
    if len(ddf):
        ddf['node']=data.index
        ddf['cluster_id']=[fnames[l] for l in labels]
        for i,idx in enumerate(table.Center):
            df=ddf[ddf.cluster_id==idx]
            fig.add_trace(go.Scattermapbox(mode='markers',lon = df.lon, lat = df.lat,
                          text=df.node, hoverinfo='text', name=f"Area {i+1}", 
                          marker = { 'size': 8}))
        
        fig.update_layout( mapbox = { 'style': "stamen-terrain", 'zoom':2.7,
        'center': {'lon': ddf.lon.mean(), 'lat': ddf.lat.mean() }},showlegend = True, 
        height=600, width=750,)
    #fig = px.scatter_geo(ddf,lat=ddf.lat, lon=ddf.lon,color=ddf.cluster_id,hover_name="node")
    #fig.update_layout(geo=dict(lataxis={'range':(23,55)},lonaxis={'range':(-110,-62)}))
    curves=groupsPlot(data,table,ddf)
    return fig,table,curves
    
import glob
if __name__=="__main__":
    data=pd.read_csv('data/sample.csv', sep=',',index_col=0)
    if 'lat' in data.columns and 'lon' in data.columns:
       ddf=data[['lat','lon']]
       data=data.iloc[:,0:-2]
    if 'time' in data.index:
        data=data[1:]
        if len(ddf): 
            ddf=ddf[1:]
    rev,labels=hierarchical(data,'single')
    print(labels)
    plt.plot(rev)
    plt.plot([0,len(rev)-1],[rev[0],rev[-1]])
    D=np.zeros(len(rev))
    n=len(rev)-1
    a=(rev[n]-rev[0])/(n-0)
    b=n*rev[0]/(n-0)
    print("eeeee",b,-rev[0],rev[-1])
    for i in range(len(D)):
        D[i]=np.abs(-rev[i]+a*i+b)/np.sqrt(a**2+1)
    print(np.argmax(D))
        


    plt.show()
    #fn="data/Gentrip/data_20200716_160955.csv"
    #gd='data/units_gps.csv'
    #for fn in glob.glob("data/Gentrip/*.csv"):
        #fn="data/Gentrip/data_20200908_132512.csv"
    #    f,table=plotMap(fn,gd)
    #print(table.to_html(index=None))
