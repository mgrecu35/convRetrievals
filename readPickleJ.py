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
