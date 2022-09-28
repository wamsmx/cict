import numpy as np
import pandas as pd

def sortrows(a):
	new=np.zeros(a.shape)
	b=np.sort(a,axis=0)
	for i in range(b.shape[0]):
		c=np.where(b[i,0]==a[:,0])
		new[i,:]=a[c[0][0],:]
	return new

def get_bounds(data, ti,tf=np.inf):
	df=data.set_index('Node_name')
	Tr=df.loc['t'].values
	print(np.max(Tr))
	T=[(i,t) for i,t in enumerate(Tr) if ti<=t and t<=tf]
	Fs=1 /(Tr[1]-Tr[0])
	print(ti,tf,T[0],T[-1])
	return Fs,T   

def Prony(data,Ka,ti,tf=np.inf):
	Fs,T=get_bounds(data,ti,tf)
	Ti,Tfi=T[0]
	Tf,Tff=T[-1]
	data=data.to_numpy()
	print(data)
	t=data[0,1:]
	stna=data[1:,1:]

	TS=1/Fs

	ns=stna.shape
	st1=np.zeros(ns)
	for k in range(ns[0]):
		st1[k,:]=stna[k,:]-np.mean(stna[k,:])
		pass
	st=st1[:,Ti:Tf]
	ns,Nm=st.shape #orden Prony

	#--------Prony 
	longg=np.arange(Ka+1,Nm+1)
	longgs=longg.shape[0]
	cte=0
	print("CCCCCCC",Ka, type(Ka))
	Zigma=np.zeros([ns*longgs,Ka])
	sn=np.zeros([ns*longgs,1])

	for l in range(ns):
		for m in range(longgs):
			for k in range(Ka):
				Zigma[m+cte,k]=st[l,Ka-k+m-1]
				pass
			sn[m+cte,0]=st[l,Ka+m]
			pass
		cte=cte+longgs
		pass

	a=np.dot(np.linalg.pinv(Zigma),sn) 
	z=np.ones(a.shape[0]+1)
	z[1:]=-1*a.T

	zi=np.roots(z)
	Z=np.zeros([Nm,Ka],dtype=complex)

	for m in range(Nm):
		for k in range(Ka):
			Z[m,k]=zi[k]**(m)
			pass
		pass
	B=np.dot(np.linalg.pinv(Z),np.transpose(st))
	sk=np.dot(Z,B)

	#iDENTIFICACION
	lamdas=np.log(zi)/TS
	Amp=2*abs(B)
	tetas=np.angle(B)
	damp=-1*np.real(lamdas)
	omega=np.imag(lamdas)
	freq=omega/(2*np.pi)
	damp_ratio=100*damp/omega
	idx=[i for i,v in enumerate(freq) if v>0 and v<10]
	resultado={
		'Frequency':freq[idx],
		'Damping':damp[idx],
		'Damp.ratio':damp_ratio[idx],
	}
	return resultado, (Ti,Tf)

def ERA(data,threshold,ti,tf=np.inf):
	Fs,T=get_bounds(data,ti,tf)
	T0,Tfi=T[0]
	T_end,Tff=T[-1]
	data=data.to_numpy()
	t=data[0,1:]
	stna=data[1:,1:]
	
	print("XXXXXXXXX",stna.shape)
	
	stax=stna
	Ts=1/Fs
	ns=stax.shape

	st1=np.zeros(ns)

	for k in range(ns[0]):
		st1[k,:]=stna[k,:]-np.mean(stna[k,:])
		pass

	st=np.transpose(st1[:,T0:T_end+1])
	Nm=st.shape
	tao=Ts

	r=round(Nm[0]/2)
	c1=0
	c2=1

	H0=np.zeros([ns[0]*r,r])
	H1=np.zeros([ns[0]*r,r])
	for m in range(r):
		for l in range(r):
			H0[c1*ns[0]:c2*ns[0],l]=np.transpose(st[l+c2,:])
			H1[c1*ns[0]:c2*ns[0],l]=np.transpose(st[l+c2+1,:])
			pass
		c1=c1+1
		c2=c2+1
		pass
	#SVD


	U,S,V=np.linalg.svd(H0)

	#threshold
	s_11=S
	S_11=np.sum(s_11)
	En1=0.1
	flag1=0
	Ene1=0
	energia=threshold #entrada
	while (En1<energia):
		flag1=flag1+1
		Ene1=Ene1+s_11[flag1-1]
		En1=Ene1/S_11
		pass
	nx=flag1 ##Numero de modos (pares conjugados)
	
	Sn=np.zeros([nx,nx])
	np.fill_diagonal(Sn,S)

	Un=U[:,:nx]
	Vt=V
	Vn=Vt[:nx,:]

	b=Sn**(-1/2)
	b[np.isinf(b)]=0
	A=np.dot(np.dot(np.dot(np.dot(b,np.transpose(Un)),H1),np.transpose(Vn)),b)

	z=np.linalg.eig(A)

	#identification
	lambda1=np.log(z[0])*Fs
	lambdaR=np.real(lambda1)
	lambdaI=np.imag(lambda1)

	freq=lambdaI/(2*np.pi)
	damp=-lambdaR
	damp_ratio=100*damp/lambdaI
    
	idx=[i for i,v in enumerate(freq) if v>0 and v<10]
	resultado={
		'Frequency':freq[idx],
		'Damping':damp[idx],
		'Damp.ratio':damp_ratio[idx],
	}
	return resultado,(T0,T_end)


def Matrix_Pencil(data,threshold,ti,tf=np.inf):
	Fs,T=get_bounds(data,ti,tf)
	T0,Tfi=T[0]
	T_end,Tii=T[-1]
	Fs=60
	data=data.to_numpy()
	print(data.shape)
	t=data[0,1:]
	stna=data[1:,1:]
	
	stax=stna
	Ts=1/Fs
	ns=stax.shape

	st1=np.zeros(ns)

	for k in range(ns[0]):
		st1[k,:]=stna[k,:]-np.mean(stna[k,:])
		pass

	st=np.transpose(st1[:,T0:T_end+1])
	Nm=st.shape
	tao=Ts

	r=round(Nm[0]/2)
	c1=0
	c2=1

	H0=np.zeros([ns[0]*r,r])
	H1=np.zeros([ns[0]*r,r])
	for m in range(r):
		for l in range(r):
			H0[c1*ns[0]:c2*ns[0],l]=np.transpose(st[l+c2,:])
			H1[c1*ns[0]:c2*ns[0],l]=np.transpose(st[l+c2+1,:])
			pass
		c1=c1+1
		c2=c2+1
		pass

	U,S,V=np.linalg.svd(H0)

	#threshold
	s_11=S
	S_11=np.sum(s_11)
	En1=0.1
	flag1=0
	Ene1=0
	energia=threshold #entrada
	while (En1<energia):
		flag1=flag1+1
		Ene1=Ene1+s_11[flag1-1]
		En1=Ene1/S_11
		pass
	nx=flag1

	V=np.transpose(V)
	Vnt=V[:,:nx]

	Vs1=Vnt[0:(r-1),:]
	Vs2=Vnt[1:r,:]

	Y1=np.dot(np.transpose(Vs1),Vs1)
	Y2=np.dot(np.transpose(Vs2),Vs1)
	Y1=np.transpose(Y1)
	z=np.linalg.eig(np.dot(np.linalg.inv(Y1),Y2))


	lambda1=np.log(z[0])*Fs
	lambdaR=np.real(lambda1)
	lambdaI=np.imag(lambda1)

	freq=lambdaI/(2*np.pi)
	damp=-lambdaR
	damp_ratio=100*damp/lambdaI

	#parameters1=np.vstack([freq,damp])
	#parameters1=np.vstack([parameters1,damp_ratio])

	#freq_pos=np.argwhere((parameters1[0]>0.2) & (parameters1[0]<1.0))
	#parameters2=parameters1[:,freq_pos]
	#parameters=sortrows(np.transpose(parameters2[:,:,0]))


	#resultado={
	#	'Frecuency':parameters[:,0],
	#	'Damping':parameters[:,1],
	#	'Damp.ratio':parameters[:,2]
	#}
	idx=[i for i,v in enumerate(freq) if v>0 and v<10]
	resultado={
		'Frequency':freq[idx],
		'Damping':damp[idx],
		'Damp.ratio':damp_ratio[idx],
	}
	return resultado,(T0,T_end)

if __name__=="__main__":
    filename = 'Sample_Multi_methods.csv'
    data = pd.read_csv(filename)    
    T0=64 #Starting time
    T_end=2000  #Ending time
    Threshold=0.5 # Threshold
    Data_Matrix_Pencil=Matrix_Pencil(data,T0,T_end,Threshold)
    #Data_Matrix_Pencil=Prony(data,T0,T_end,5)
    print(pd.DataFrame(Data_Matrix_Pencil))
