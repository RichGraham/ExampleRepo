# -*- coding: utf-8 -*-
#I added a new comment
import math
import scipy
from scipy import optimize
import pandas as pa
import numpy as np
import re
import sys
import os

print ('Modules loaded')

#Extract the command line arguments and 
barrierFile= sys.argv[1]
FelFile= sys.argv[2]
pGOFile= sys.argv[3]

#construct the file names
base=os.path.basename(barrierFile)
barrierName=os.path.splitext(base)[0]
base=os.path.basename(FelFile)
FelName=os.path.splitext(base)[0]

#outFileName="Daniel-Thry.dat"
outFileName=barrierName+"-"+FelName+"-Strd.dat"


#Read quiescent barrier input file
with open(barrierFile, 'r') as infile:
    for line in infile:
        if 'ebulk' in line: 
            E0=[float(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print ('Found ebulk ='+str(E0))
        if 'esurface' in line: 
            mus=[float(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print ('Found esurface ='+str(mus))
        if 'Kappa' in line: 
            Kappa0=[float(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print ('Found Kappa0 ='+str(Kappa0))
        if 'MAX_SEGS' in line: 
            maxNT=[int(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print ('Found MAX_SEGS ='+str(maxNT))

ratio=2.0
with open(pGOFile, 'r') as infile:
    for line in infile:
        if 'Qs' in line: 
            Qs0=[float(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print ('Found Qs0 ='+str(Qs0))
        if 'ratio' in line: 
            ratio=[float(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]

print ('Found ratio ='+str(ratio))

            
#Load in datafile specifying Phi DFel
count = len(open(FelFile).readlines(  ))
grid = pa.read_csv(FelFile, delim_whitespace=True, header=None,skiprows=[count-1,0],error_bad_lines=False)  
grid= grid.values
df = grid[:,-1] 
phi = grid[:,-2]

numc=phi.size
thetaMin=1e-300

#a_r
arsq=9.0/16.0*math.pi
ar=math.sqrt(arsq)

def wfun( phi, Df, P, B, NS):
    return Qs*phi/( NS + (Qs-NS)*B*np.exp(-LL*Df - Df**2/2.0/kappa  +  Df/kappa*P))

def afun(x, NS):
    P=x[0]
    B=x[1]
    sum1=0.0
    sum2=0.0
    for i in range(numc):
        wi=wfun( phi[i], df[i],P,B ,NS)
        sum1+=wi
        sum2+=wi*df[i]
    return np.array([sum1-1.0 , sum2-P])

def Free2(Ns,NT):
    global LL,Qs,kappa, Pprevious, Bprevious
    NS=Ns[0]
    LL=NT/NS
    Qs=Qs0*NS
    kappa=Kappa0+1.0/LL**2
    
    sol = scipy.optimize.root(afun, np.array([ max(0.00001,Pprevious), Bprevious]),  args=(NS),method='hybr', jac = False)
    P = sol.x[0]
    B = sol.x[1]
    Pprevious=P
    Bprevious=B
    
    w=np.zeros(numc)
    theta=np.zeros(numc)
    sum1=0
    sum2=0
    sum3=0
    sum4=0
    sum5=0
    for i in range(numc):
        w[i]=wfun( phi[i], df[i], P ,B, NS)
        theta[i] = max(thetaMin,Qs/(Qs - NS)*phi[i] - NS/ (Qs - NS)*w[i])
        sum1 += theta[i]*np.log(theta[i])
        sum2 += phi[i]*np.log(phi[i])
        sum3 += w[i]*np.log(w[i])
        sum4 += w[i]*df[i]
        sum5 += w[i]*df[i]**2
    
    #print w
    #print theta
    
    FF = (Qs-NS)*sum1 - Qs*sum2 + NS*sum3 - 0.5*(NS-1)*np.log(2.0*math.pi/kappa)+0.5*np.log(NS) \
        - NT*sum4 - 0.5*NS/kappa*(sum5 - sum4**2) - E0*NT
    #surface terms
    aspect=NS**3/NT**2/arsq
    if aspect<1:
        ep=math.sqrt(1.0-aspect)
        Stil=2.0*NS+2.0*ar*NT*math.asin(ep)/ep/math.sqrt(NS)
    elif aspect>1:
        eps0=math.sqrt(1.0-1.0/aspect)
        Stil=2.0*NS+arsq*NT**2*math.log((1.0+eps0)/(1.0-eps0))/eps0/NS**2
    else:
        Stil=2.0*NS+2.0*ar*NT/math.sqrt(NS)
    FF+=mus*Stil
    return FF

def Free1(NT):
    x0 = NSprevious
    res=scipy.optimize.minimize(Free2, x0,method='Nelder-Mead' ,args=(NT))
    #res=scipy.optimize.minimize_scalar(Free2, bounds=(1,0.999999*NT), args=(NT), method='bounded')
    print res
    return res.fun
    
def Freefluc(NT):
    global NSprevious
    x0 = NSprevious
    res=scipy.optimize.minimize(Free2, x0,method='Nelder-Mead' ,args=(NT))
    #print res
    #res=scipy.optimize.minimize_scalar(Free2, bounds=(1,0.999999*NT), args=(NT), method='bounded')
    nsmid=res.x[0]
    NSprevious=nsmid
    ##second derivative
    d2fdn2=(Free2([nsmid+0.1],NT)+Free2([nsmid-0.1],NT)-2*Free2([nsmid],NT))/0.01
    fNT=res.fun+math.log(d2fdn2/2/math.pi)
    return fNT

    
NTlist=[]
Flist=[]
Fluclist=[]
Barrierlist=[]

BestBarrier=-10000.0
NTlist=[]
Flist=[]
Fluclist=[]
Barrierlist=[]
NSprevious=1.1
Pprevious=0.5*(max(df)-min(df))
Bprevious=1.0


#Wipe the output file
f= open(outFileName,"w+")
f.close

for i in range(int((maxNT-1)/ratio)):
    NT=2.0+i*ratio
    NTlist.append(NT)
    ans=Freefluc(1.0*NT)
    Fluclist.append(ans)

    
    if(ans>BestBarrier): 
        BestBarrier=ans

    f= open(outFileName,"a+")
    f.write(str(NT)+' '+str( Fluclist[i] )+'\n')
    f.close()

    if( BestBarrier-ans > 1.0):
        break


    if( i % 100 == 0):
        print (str(NT+Ns)+ " / "+  str(maxNT) )

