import pickle

def readpickle(fname):
    dataOut=pickle.load(open(fname,'rb'))

    [clutFL,reliabFlagNSL,reliabFlagMSL,pathAttenNSL,\
     pathAttenMSL,binSfcNSL,binSTopL,binZeroDegL,\
     zKuL,zKaL,cmbSfcRainL,gvSfcRainL,cmbPWCL,piaL]=dataOut
    #print(zKuL[0].shape)
    return clutFL,reliabFlagNSL,reliabFlagMSL,pathAttenNSL,\
     pathAttenMSL,binSfcNSL,binSTopL,binZeroDegL,\
     zKuL,zKaL,cmbSfcRainL,gvSfcRainL,cmbPWCL,piaL


clutFL,reliabFlagNSL,reliabFlagMSL,pathAttenNSL,\
    pathAttenMSL,binSfcNSL,binSTopL,binZeroDegL,\
    zKuL,zKaL,cmbSfcRainL,gvSfcRainL,cmbPWCL,piaL=readpickle("dataVN_3D.pklz")

ic=0

from sklearn.neighbors import NearestNeighbors
neigh = NearestNeighbors(n_neighbors=10)
from netCDF4 import Dataset
fh=Dataset("simZDataset_3.nc")
zT=fh["zSim"][:,:]
neigh.fit(zT)
rmsL=[]
from numpy import *
import matplotlib.pyplot as plt
h1=(176-arange(99,165))*0.125
for i,zKu in enumerate(zKuL[2]):
    if clutFL[2][i,1,1]>165:
        zObs=zKu[1,1,99:165]
        ic+=1
        zObs[zObs<0]=0
        rms,ind=neigh.kneighbors(zObs.reshape(1,-1))
        rmsL.append(rms[0][0])
        if rms[0][0]>40:
            plt.plot(zT[ind[0],:].mean(axis=0),h1,color='pink')
            #plt.plot(zT[ind[0][1],:],h1,color='pink')
            #plt.plot(zT[ind[0][2],:],h1,color='pink')
            plt.plot(zObs,h1,color='blue')
            plt.show()


from scipy.ndimage import gaussian_filter
noise=random.randn(66)*1.5
noise=gaussian_filter(noise,3)
plt.plot(noise,h1)
