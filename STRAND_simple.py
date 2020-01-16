# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 13:15:57 2018

@author: phydjr
"""

import math
import scipy
#import matplotlib.pyplot as plt
from scipy import optimize
import pandas as pa
import numpy as np
import re
import sys
import os

print 'Modules loaded'

#Extract the command line arguments and 
barrierFile= sys.argv[1]
FelFile= sys.argv[2]

#construct the file names
base=os.path.basename(barrierFile)
barrierName=os.path.splitext(base)[0]
base=os.path.basename(FelFile)
FelName=os.path.splitext(base)[0]

outFileName="Daniel-Thry.dat"
quiescentname=barrierName+"-"+FelName+"-calc.dat"
outFileName=barrierName+"-"+FelName+"-Thry.dat"


#Read quiescent barrier input file
with open(barrierFile, 'r') as infile:
    for line in infile:
        if 'ebulk' in line: 
            E0=[float(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print 'Found ebulk ='+str(E0)
        if 'esurface' in line: 
            mus=[float(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print 'Found esurface ='+str(mus)
        if 'MAX_SEGS' in line: 
            maxNT=[int(s) for s in re.findall("[-+]?[.]?[\d]+(?:,\d\d\d)*[\.]?\d*(?:[eE][-+]?\d+)?", line)][0]
            print 'Found MAX_SEGS ='+str(maxNT)


#Load in accurate quiescent landscape
grid = pa.read_csv(quiescentname, delim_whitespace=True, header=None)
grid= grid.values
trueQuiescent = grid[:,-1]
index = grid[:,-2]


#Load in datafile specifying Phi DFel

count = len(open(FelFile).readlines(  ))
grid = pa.read_csv(FelFile, delim_whitespace=True, header=None,skiprows=[count-1,0],error_bad_lines=False)  
grid= grid.values
df = grid[:,-1] 
phi = grid[:,-2]



#numc=10
#array for phi_i
#phi=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
#array for deltaf_i
#df=[0.3,0.1,0.1,0,0,0,0,0,0,0]
#numc=2
#phi= np.array([0.2,0.8])
#df= np.array([0.1,0])


numc=phi.size

#bulk and surface energy
#E0=0.3
#mus=1.4
#maxNT=300


#a_r
arsq=9.0/16.0*math.pi
ar=math.sqrt(arsq)

#setting up some parameters
#edf=[0,0,0,0,0,0,0,0,0,0]
edf = [0] * numc
edfmax=0
logphi=[]

for i in range(numc):
    edf[i]=math.exp(df[i])
    logphi.append(math.log(phi[i]))
    edfmax=max(edfmax,edf[i])
iedfmax=0.999999999999999/edfmax

#def afun(A,LL):
#    sum1=0.0
#    sum2=0.0
#    for i in range(numc):
#        tem=1.0-A*edf[i]
#        sum1+=phi[i]*edf[i]/tem
#        sum2+=phi[i]*edf[i]/tem**2
#    tem=LL*sum1/sum2-1.0
#    return tem

def Free2(NS,NT):
    LL=NT/NS

#    A=scipy.optimize.brenth(afun,0,iedfmax,args=(LL))

    sum1=0.0
    for i in range(numc):
        #tem=1.0-A*edf[i]
	#sum1+=phi[i]*edf[i]/tem

        tem=math.exp(LL) + edf[i]
        sum1+=phi[i]*tem

#    AB=1.0/sum1
#    w=[]
#    v=[]

    denom = sum1
    w=[]

    for i in range(numc):
#        tem=1.0-A*edf[i]
#        w.append(AB*phi[i]*edf[i]/tem)
#        v.append(w[i]/tem/LL)

         w.append(phi[i]*(math.exp(LL)+edf[i])/denom)


    # now, put together free energy terms
    sum1=0.0
    for i in range(numc):
#        logw=math.log(w[i])
#        logv=math.log(v[i])
#        logc=math.log(v[i]-w[i]/LL)
#        sum1+=w[i]*(2*logw-logphi[i])/LL-v[i]*logv+(v[i]-w[i]/LL)*logc-v[i]*df[i]    

        logw=math.log(w[i])
        sum1+=w[i]*(logw-logphi[i]-LL*df[i])

    #FF=NT*sum1-NS*math.log(LL)-NT*E0

    #LEADING ORDER
    logNT = math.log(NT)
    logNS = math.log(NS)
    logc = math.log(NT-NS)
    FN=NT*logNT-(NT-NS)*logc-NS*logNS

    #STIRLING APPROX
    #logNT = math.log(NT-1.0)
    #logNS = math.log(NS-1.0)
    #logc = math.log(NT-NS)
    #FN=(NT-1.0)*logNT-(NT-NS)*logc-(NS-1.0)*logNS

 
    #FULL FACTORIAL FORMULA
    #NTfact = math.factorial(NT-1)
    #NSfact = math.factorial(int(NS)-1)
    #NTSfact = math.factorial(NT-int(NS))  
    #logNT = math.log(NTfact)
    #logNS = math.log(NSfact)
    #logNTS = math.log(NTSfact)
    #FN = logNT - logNTS - logNS


    FF=NS*sum1-NT*E0-FN

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
        
def Freequi(NS,NT):
    LL=NT/NS    
    FF=-NS*math.log(LL)-NT*E0+(NT-NS)*math.log(1.0-1.0/LL)
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

#def Free1(NT):
#    res=scipy.optimize.minimize_scalar(Free2, bounds=(1,0.999999*NT), args=(NT), method='bounded')
#    return res.fun
    
def Freefluc(NT):
    res=scipy.optimize.minimize_scalar(Free2, bounds=(1,0.999999*NT), args=(NT), method='bounded')
    nsmid=res.x
    ##second derivative
    d2fdn2=(Free2(nsmid+0.1,NT)+Free2(nsmid-0.1,NT)-2*Free2(nsmid,NT))/0.01
    fNT=res.fun+math.log(d2fdn2/2/math.pi)
    return fNT


#def Free1qui(NT):
#    res=scipy.optimize.minimize_scalar(Freequi, bounds=(1,0.999999*NT), args=(NT), method='bounded')
#    return res.fun
    
def Freeflucqui(NT):
    res=scipy.optimize.minimize_scalar(Freequi, bounds=(1,0.999999*NT), args=(NT), method='bounded')
    nsmid=res.x
   ##second derivative
    d2fdn2=(Freequi(nsmid+0.1,NT)+Freequi(nsmid-0.1,NT)-2*Freequi(nsmid,NT))/0.01
    fNT=res.fun+math.log(d2fdn2/2/math.pi)
    return fNT

#def Freesum(NT):
#    sum=0.0
#    NS=1
#    while NS<NT:
#        sum=sum+math.exp(-Free2(NS,NT))
#        NS=NS+1
#    fren=-math.log(sum)
#    return fren

#def Freesumqui(NT):
#    sum=0.0
#    NS=1
#    while NS<NT:
#        sum=sum+math.exp(-Freequi(NS,NT))
#        NS=NS+1
#    fren=-math.log(sum)
#    return fren

    
#NT=100.0
#res=scipy.optimize.minimize_scalar(Free2, bounds=(1,0.999999*NT), args=(NT), method='bounded')
#print(res)
#print(res.x)
#print(res.fun)

#NSlist=[]
#Flist=[]
#Fqlist=[]

#for i in range(99):
#    NS=(i+1)*1.0
#    NSlist.append(NS)
#    Flist.append(Free2(NS,NT))
#    Fqlist.append(Freequi(NS,NT))
    

#plt.plot(NSlist,Flist)
#plt.plot(NSlist,Fqlist,color='red')
#plt.plot([res.x],[res.fun],'x')
#plt.show()


NTlist=[]
Flist=[]
Fluclist=[]
Barrierlist=[]


#Wipe the output file
f= open(outFileName,"w+")
f.close

for i in range((maxNT-1)/2):
    NT=(i+1)*2.0

    NTlist.append(NT)
#    Flist.append(Free1(NT))
#    Flist.append(Freesum(NT))
#    Fluclist.append(Freefluc(NT))

#   SUM OVER FOR LIST
#    Flist.append(Freesum(NT)-Freesumqui(NT))

#   OPTIMIZE FREE2 FOR LIST
    #Fluclist.append(Free1(NT)-Free1qui(NT))
#   OPTIMIZE FREE2 WITH FLUCS LIST?
    Fluclist.append(Freefluc(NT)-Freeflucqui(NT))

    f= open(outFileName,"a+")
    #print str(NT)+' '+str(index[NT-1] )
    f.write(str(NT)+' '+str( trueQuiescent[NT-1]+Fluclist[i] )+'\n')
    f.close()
#    print str(NT)+' '+str(Free1(NT))+' '+str(Freequi(NT/2,NT))+' '+str(Flist[i])+' '+str(Freesumqui(NT))#
     

#plt.plot(NTlist,Flist,color='red')
#plt.plot(NTlist,Fluclist,color='blue')
#plt.show()

