import numpy as np
from lmfit import Model

def rational_R_ge(q2,a,b,R):
    return a*(1.-R*R*q2/6.+b*q2)/(1+b*q2)

class xy_bootstrap(object):
    def read_in(self,filename,seed=7):
        data=np.loadtxt(filename,skiprows=1)
        self._q2=data[:,0]
        self._ge=data[:,1]
        self._dge=data[:,2]
        np.random.seed(seed)

    def single_fit(self):
        print("single R fit")
        ge=self._ge; q2=self._q2; dge=self._dge
        #print("ge q2 dge shape" + str(ge.shape)+ " , " + str(q2.shape) + " , " + str(dge.shape))
        modelr=Model(rational_R_ge)
        result=modelr.fit(ge,q2=q2,a=1,b=-0.01,R=1.0,weights=1/dge)
        print(result.fit_report())
        print("a= " + str(result.best_values['a']))
        print("b= " + str(result.best_values['b']))
        print("R= " + str(result.best_values['R']))

    def bs_replace(self, loop, slen):
        ge=self._ge; q2=self._q2; dge=self._dge
        print("Bootstrap with replacement: nbin tot= " + str(len(q2)) + " , sample nbin= " + str(slen) + " , loop= " + str(loop))
        assert(slen <= len(q2))

        a_boot=[]
        b_boot=[]
        radius_boot=[]
        r_err_boot=[]

        for l in range(loop):
            if(np.mod(l,100) == 0):
                print("Running: " + str(l) + " | " + str(loop))

            idx=np.random.randint(0,len(q2),size=slen)
            idx=np.sort(idx)
            subq2=q2[idx]; subge=ge[idx]; subdge=dge[idx]
            if(np.mod(l,100) == 0):
                print("  sizes q2, ge, dge: " + str(subq2.size) + " , " + str(subge.size) + " , " + str(subdge.size))
            modelr=Model(rational_R_ge)
            result=modelr.fit(subge,q2=subq2,a=1,b=-0.01,R=1.0,weights=1/subdge)
            a_boot.append(result.best_values['a'])
            b_boot.append(result.best_values['b'])
            radius_boot.append(result.best_values['R'])
            r_err_boot.append(result.params['R'].stderr)

        print('a: ',np.mean(a_boot),' +/- ',(np.std(a_boot)), ' | size= ', len(a_boot))
        print('b: ',np.mean(b_boot),' +/- ',(np.std(b_boot)), ' | size= ', len(b_boot))
        print('R: ',np.mean(radius_boot),' +/- ',(np.std(radius_boot)), ' | size= ', len(radius_boot))
        print('r err: ',np.mean(r_err_boot),' +/- ',(np.std(r_err_boot)), ' | size= ', len(r_err_boot))

        return [np.mean(radius_boot), (np.std(radius_boot)), len(radius_boot)]


    def bs_noreplace(self, loop, slen):
        ge=self._ge; q2=self._q2; dge=self._dge
        print("Bootstrap without replacement: nbin tot= " + str(len(q2)) + " , sample nbin= " + str(slen) + " , loop= " + str(loop))
        assert(slen <= len(q2))
        
        a_boot=[]
        b_boot=[]
        radius_boot=[]
        r_err_boot=[]
        D={}

        for l in range(loop):
            if(np.mod(l,100) == 0):
                print("Running: " + str(l) + " | " + str(loop))
    
            perm=np.arange(len(q2))
            np.random.shuffle(perm)
            idx=perm[0:slen]
            idx=np.sort(idx)

            ts=str(idx)
            if (ts in D): continue
            D[ts]=1

            subq2=q2[idx]; subge=ge[idx]; subdge=dge[idx]
            if(np.mod(l,100) == 0):
                print("  sizes q2, ge, dge: " + str(subq2.size) + " , " + str(subge.size) + " , " + str(subdge.size))
            modelr=Model(rational_R_ge)
            result=modelr.fit(subge,q2=subq2,a=1,b=-0.01,R=1.0,weights=1/subdge)
            a_boot.append(result.best_values['a'])
            b_boot.append(result.best_values['b'])
            radius_boot.append(result.best_values['R'])
            r_err_boot.append(result.params['R'].stderr)

        print('a: ',np.mean(a_boot),' +/- ',(np.std(a_boot)), ' | size= ', len(a_boot))
        print('b: ',np.mean(b_boot),' +/- ',(np.std(b_boot)), ' | size= ', len(b_boot))
        print('R: ',np.mean(radius_boot),' +/- ',(np.std(radius_boot)), ' | size= ', len(radius_boot))
        print('r err: ',np.mean(r_err_boot),' +/- ',(np.std(r_err_boot)), ' | size= ', len(r_err_boot))

        return [np.mean(radius_boot), (np.std(radius_boot)), len(radius_boot)]



# Test
xyb=xy_bootstrap()

#filename="xy_gen_Q2_GE_41_dipole.txt"
filename="xy_gen_Q2_GE_41_dipole_fluc1.txt"
#filename="xy_gen_Q2_GE_41_dipole_fluc2.txt"
xyb.read_in(filename)

xyb.single_fit()
#[Rmean,Rrms,Rsize]=xyb.bs_replace(1000, 25)
[Rmean,Rrms,Rsize]=xyb.bs_noreplace(100, 40)

print("Rmean,Rrms,Rsize= " + str([Rmean,Rrms,Rsize]))





