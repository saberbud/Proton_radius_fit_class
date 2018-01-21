import numpy as np
from lmfit import Model
import matplotlib.pyplot as plt
import statsmodels.api as sm

#
# First Just A Normal Least Sqaures Fit
#

#filename="xy_gen_Q2_GE_41_dipole.txt"
#filename="xy_gen_Q2_GE_41_dipole_fluc1.txt"
filename="xy_gen_Q2_GE_41_dipole_fluc2.txt"

data=np.loadtxt(filename,skiprows=1)
q2=data[:,0]
ge=data[:,1]
dge=data[:,2]

#def quad_ge(q2,a,b):
#    return (a+b*q2+0.01*q2**2)

def rational_ge(q2,a,b,c):
    return (a+b*q2)/(1+c*q2)

print("single fit")
model=Model(rational_ge)
result=model.fit(ge,q2=q2,a=1,b=-0.01,c=0.1,weights=1/dge)
print(result.fit_report())
print("a= " + str(result.best_values['a']))
print("b= " + str(result.best_values['b']))
print("c= " + str(result.best_values['c']))
print("R= " + str(np.sqrt(6*(-1*result.best_values['b']/result.best_values['a']+result.best_values['c']))))

#
# Now Errors From Statistical Bootstrap
#
print("\nBootstrap fit")

loop=1000

a_boot=[]
b_boot=[]
c_boot=[]
radius_boot=[]

slen=len(q2)  #full resampling
#slen=100  #choose for partial resampling

for l in range(loop):
    if(np.mod(l,100) == 0):
        print("Running: " + str(l) + " | " + str(loop))

    idx=np.random.randint(0,len(q2),size=slen)
    idx=np.sort(idx)
    subq2=q2[idx]; subge=ge[idx]; subdge=dge[idx]
    if(np.mod(l,100) == 0):
        print("  sizes q2, ge, dge: " + str(subq2.size) + " , " + str(subge.size) + " , " + str(subdge.size))

    model=Model(rational_ge)
    result=model.fit(subge,q2=subq2,a=1,b=-0.01,c=0.1,weights=1/subdge)

    a_boot.append(result.best_values['a'])
    b_boot.append(result.best_values['b'])
    c_boot.append(result.best_values['c'])
    radius_boot.append(np.sqrt(6*(-1*result.best_values['b']/result.best_values['a']+result.best_values['c'])))
    


# Bootstrap result:

print('a: ',np.mean(a_boot),' +/- ',(np.std(a_boot)), ' | size= ', len(a_boot))
print('b: ',np.mean(b_boot),' +/- ',(np.std(b_boot)), ' | size= ', len(b_boot))
print('c: ',np.mean(c_boot),' +/- ',(np.std(c_boot)), ' | size= ', len(c_boot))
print('R: ',np.mean(radius_boot),' +/- ',(np.std(radius_boot)), ' | size= ', len(radius_boot))

'''
# qq plots:
normqq=np.sort(rational_ge(subq2,result.best_values['a'],result.best_values['b'],result.best_values['c'])-subge)
sm.qqplot(normqq,line='r')
plt.show()
'''

