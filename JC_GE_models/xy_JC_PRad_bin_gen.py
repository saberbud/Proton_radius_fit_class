import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import interpolate

JC_para_0=np.loadtxt("TableGE0Parametrizaton.dat",skiprows=0)
JC_para_mu=np.loadtxt("TableGEmuHParametrizaton.dat",skiprows=0)
JC_para_co=np.loadtxt("TableGECODATAParametrizaton.dat",skiprows=0)

TableGEpParameterization=JC_para_0

class JC_class:
    def __init__(self):
        self.initD=0

    def init_JC(self):
        print("Initializing.")
        x=[]
        y=[]

        for i in range(len(TableGEpParameterization)):
            x.append(TableGEpParameterization[i][0]*25.7)
            y.append(TableGEpParameterization[i][1])

        self.Jose=interpolate.interp1d(x,y,kind='cubic')
        self.initD=1

    def calc_val(self,xin):
        if self.initD<1:
            self.init_JC()

        return self.Jose(xin)

    def __str__(self):
        if self.initD<1:
            return "JC not initialized."
        return "JC initialized"


# JC init
JC=JC_class()
print(JC.initD)
print(JC)
#TableGEpParameterization=JC_para_co   #CODATA para
#TableGEpParameterization=JC_para_mu   #Mu para
JC.init_JC()

data= np.loadtxt("bin_set_1.txt",skiprows=1)
Q2=data[:,0]
GE_0=data[:,1]
dGE=data[:,2]

x_in=Q2
y=JC.calc_val(x_in)
print(y)

dout=np.array([Q2, y, dGE])
dout=np.transpose(dout)
print(dout.shape)
print(dout)

np.savetxt("out.txt", dout, fmt='%.9e', delimiter=' ', newline='\n')






















