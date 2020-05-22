from numba import njit

@njit
def gridlookup2(x0,xgrid):

    nx = np.shape(xgrid)[0]

    ix = 0
    for jx in range(nx):
        if x0<=xgrid[jx]:
            break
        ix = ix+1

    ix = min(max(1,ix),nx-1)

    return ix-1

@njit
def driver(r0,vmat0,mu0,Ge,Pe,knotsb,β,α,δ,znow,lnow,critin,critmu):

    vmat1 = np.zeros((nb,ne))
    gmat0 = np.zeros((nb,ne))

    # Factor Price
    m0 = lnow*(r0/α/znow)**(1/(α-1))
    w0 = (1-α)*znow*m0**(α)*lnow**(-α)

    # instanteneous utility
    umat0 = np.zeros((nb,nb,ne))
    for ie in range(ne):
        for ib in range(nb):
            for jb in range(nb):

                enow = Ge[ie]
                know = knotsb[ib]
                kp   = knotsb[jb]
                cnow = w0*enow + (1.0+r0-δ)*know - kp;

                if cnow>0:
                    umat0[jb,ib,ie] = np.log(cnow)
                else:
                    umat0[jb,ib,ie] = -1e+10

    # value function iteration
    diff   = 1e+4
    iterin = 0

    while diff>critin:

        # NOTE: np.max with axis argument not supported by numba
#         for ie in range(ne):

#             utemp = umat0[:,:,ie]
#             vcond = Pe[ie,0]*vmat0[:,0] + Pe[ie,1]*vmat0[:,1]
#             vtemp = utemp + β*vcond.reshape(nb,1) # BUG: vcond should have been (nbx1)
#             vmat1[:,ie] = np.max(vtemp,axis=0)

        for ie in range(ne):

            for ib in range(nb):

                utemp = umat0[:,ib,ie]
                vcond = Pe[ie,0]*vmat0[:,0] + Pe[ie,1]*vmat0[:,1]
                vtemp = utemp + β*vcond
                vmat1[ib,ie] = np.max(vtemp)

        diff = np.max(np.abs(vmat1-vmat0))
        iterin = iterin+1
        vmat0 = np.copy(vmat1)
#         print([iterin,diff])

    # policy function
    for ie in range(ne):

        for ib in range(nb):

            utemp = umat0[:,ib,ie]
            vcond = Pe[ie,0]*vmat0[:,0] + Pe[ie,1]*vmat0[:,1]
            vtemp = utemp + β*vcond
            jb = np.argmax(vtemp)
            gmat0[ib,ie] = knotsb[jb]

    # transition matrix
    AA = np.zeros((nb*ne,nb*ne))
    wb = np.zeros((nb,ne))

    for ie in range(ne):

        for ib in range(nb):

            know = knotsb[ib]
            kp = gmat0[ib,ie]
            kb = gridlookup2(kp,knotsb)
            wb[ib,ie] = (knotsb[kb+1]-kp)/(knotsb[kb+1]-knotsb[kb])

            for je in range(ne):

                ia = nb*ie+ib
                ja = nb*je+kb
                AA[ia,ja]   = wb[ib,ie]*Pe[ie,je]
                AA[ia,ja+1] = (1.0-wb[ib,ie])*Pe[ie,je]

    # distribution
    diffmu = 1e+4
    dist   = np.zeros(2*nb)
    mu1    = np.zeros((nb,ne))

    while diffmu>critmu:

        dist[0:nb]    = mu0[:,0]
        dist[nb:2*nb] = mu0[:,1]
#         dist = np.reshape(mu0,(2*nb,1),order='F') # NOTE: cannot use with numba
        dist = AA.T@dist
        mu1[:,0] = dist[0:nb]
        mu1[:,1] = dist[nb:2*nb]
#         mu1  = np.reshape(dist,(nb,2),order='F')

        diffmu = np.max(np.abs(mu1-mu0))
        mu0  = mu1/np.sum(mu1)

    # Calculate K
    m1 = 0.0
    for ie in range(ne):

#         m1 = m1 + mu0[:,ie].T@knotsb
        m1 = m1 + np.dot(mu0[:,ie],knotsb)

    r1 = (α)*znow*m1**(α-1)*lnow**(1-α)

    return r1, vmat0, mu0


import numpy as np
import matplotlib.pyplot as plt
import time

β  = 0.96
α  = 0.36
δ  = 1-0.92
σ  = 1.5
μ  = 0.05
pee = 0.925
puu = 0.5

# grid points
ne = 2
Ge = np.array([1.0,μ])
Pe = np.array([[pee,1-pee],
    [1-puu,puu]])

nb = 1001
kmin = 0
kmax = 20.0
#kbounds = np.array([kmin,kmax])
knotsb = np.linspace(kmin,kmax,nb)

mue = np.identity(ne)
for i in range(10000):
    mue = Pe@mue

mue = mue[0,:]
znow = 1.0
# efficiency unit of labor
lnow = Ge.T@mue

# initial distribution
vmat0 = np.zeros((nb,ne))
mu0   = np.ones((nb,ne))/(nb*ne)

# initial value of r
mnow = lnow*(α*β/(1.0-β*(1.0-δ)))**(1.0/(1.0-α))
mnow = 5.2074 # from My_Aiyagari.m
m0 = mnow
r0 = (α)*znow*mnow**(α-1)*lnow**(1-α);


start = time.time()
critin  = 1e-4
critmu  = 1e-8
critout = 1e-3
diffout = 1e+3
damp = 0.01
iter = 0

# for test
# r1,vmat0,mu0 = driver(r0,vmat0,mu0,Ge,Pe,knotsb,β,α,δ,znow,lnow,critin,critmu)

while diffout>critout:

    r1,vmat0,mu0 = driver(r0,vmat0,mu0,Ge,Pe,knotsb,β,α,δ,znow,lnow,critin,critmu)

    diffout = np.abs(np.log(r1)-np.log(r0))
    iter = iter+1
    print("iter = %4d, diff = %5.6f, oldr = %5.6f, newr = %5.6f" % (iter,diffout,r0,r1))

    # Update K
    r0 = damp*r1 + (1.0-damp)*r0

eptime = time.time()-start
print("Elapsed time is %5.10f." % eptime)
