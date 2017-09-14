from __future__ import division
import numpy as np
import scipy as sc
from scipy.optimize import minimize_scalar
import functools

def XYZ2OE(XX,constants):
    XX=XX.reshape((-1,1))

    X=X0[0:3]
    V=X0[3:6]

    MU=constants['MU']

    rm = np.norm( X); 
    vm = np.norm(V);

    a = 1/(2/rm - vm^2/mu)

    hbar = np.cross(X,V,axisa=0, axisb=0, axisc=0)
    cbar = np.cross(V,hbar,axisa=0, axisb=0, axisc=0)-MU*X/rm

    e = np.norm(cbar)/MU

    hm = np.norm(hbar)
    ih = hbar/hm;

    ie = cbar/(MU*e)

    ip = np.cross(ih,ie,axisa=0, axisb=0, axisc=0)

    i = np.arccos(ih[2,0])

    w = np.arctan2(ie[2,0],ip[2,0])

    Om = np.arctan2(ih[0,0],-ih[1,0])

    sig = np.dot(X.transpose(),V)/np.sqrt(MU)

    E = np.arctan2(sig/np.sqrt(a),1-rm/a);

    return np.array([a,e,E,w,i,Om]).reshape((-1,1))

def keplerSolver(M,e,constants):
    tol=constants['tol']

    # % Compute Eccentric Anomayl
    if (M > -np.pi and M < 0 or M > np.pi):
        E1 = M - e;
    elif (M < -np.pi or M > 0 and M < np.pi):
        E1 = M + e
    elif (M == 0):
        E1 = M + e

    E = E1 + (M - E1 + e*np.sin(E1))/(1 - e*np.cos(E1))

    while (abs(E - E1) > tol):
        E1 = E
        E = E1 + (M - E1 + e*sin(E1))/(1 - e*cos(E1))

    E = (E/(2*np.pi) - np.floor(E/(2*np.pi)))*2*np.pi

    f = 2*np.arctan2(np.sqrt(1+e)*np.tan(E/2),np.sqrt(1-e))

    return (E,f)


def OE2XYZ(Xorb,constants):
    tol=constants['tol']
    a,e,E,w,inc,Om=Xorb.reshape(1,-1)[0]
    M=E-e*np.sin(E)

    tol     = 1e-10

    p       = a*(1 - e^2)

    E,f = keplerSolver(M,e,constants)

    

    c_f     = np.cos(f);
    s_f     = np.sin(f);
    r_pqw   = np.divide(np.array([ p*c_f,p*s_f,0]) , (1 + e*c_f) ).reshape(-1,1)
    v_pqw   = np.multiply(np.array([ -s_f,(e + c_f),0]),np.sqrt(mu/p) ).reshape(-1,1)


    s_Om    = sin(Om);
    c_Om    = cos(Om);
    s_w     = sin(w);
    c_w     = cos(w);
    s_i     = sin(inc);
    c_i     = cos(inc);

    ROT  = [[c_Om*c_w-s_Om*s_w*c_i,-c_Om*s_w-s_Om*c_w*c_i,s_Om*s_i],
        [s_Om*c_w+c_Om*s_w*c_i,-s_Om*s_w+c_Om*c_w*c_i,-c_Om*s_i],
        [s_w*s_i,c_w*s_i,c_i]]
    
    ROT=np.array(ROT)

    r       = np.dot(ROT,r_pqw)
    v       = np.dot(ROT,v_pqw)
    X=np.vstack((r,v))

    return X





def orbit_propagate_kepler(t,X0,constants):
    X0=X0.reshape((-1,1))
    MU=constants['MU']

    OE0_mu = XYZ2OE(X0,MU);
    a0,ecc0,E0,w0,inc0,Omega0=OE0_mu.reshape(1,-1)[0]

    M0=E0-ecc0*np.sin(E0);
    # X= OE2XYZ(OE0_mu,constants);
    # Xorb0_mu=[r;v];

    n=np.sqrt(MU/a0^3);

    M=n*t + M0;
    E=keplerSolver(M,ecc0,constants);
    Xorb=np.array([a0,ecc0,E,w0,inc0,Omega0]).reshape(-1,1)
    X= OE2XYZ(Xorb,constants);
    return (X,Xorb)

def transcendentalfunc(y):
    if y>0:
        sqrty=np.sqrt(y)
        C=(1-np.cos( sqrty ))/y
        S=(np.sqrt(y)-np.sin(np.sqrt(y)))/np.sqrt(y^3)
        
    elif y<0:
        sqrtmy=np.sqrt(-y)
        C=(np.cosh(sqrtmy)-1)/(-y)
        S=(np.sinh(sqrtmy)-sqrtmy)/np.sqrt(-y^3)
    else:
       C=0.5
       S=0.1667
    return (C,S)


def univkepler_err(MU,alpha,sig0,normr0,t,t0,x):

    [C,S]=transcendentalfunc(alpha*x^2);
    err=sig0*x**2*C+(1-normr0*alpha)*x**3*S+normr0*x-np.sqrt(MU)*(t-t0);
    return err

def universalkepler(a,t0,X0,t,options):
    x=1
    F=functools.partial(univkepler_err,)
    x=minimize_scalar(F)

    return x

def forwardSTM(t,X0,computePHI,computeProb,probFunc0,options,constants):
    # tol=constants['tol']

    X0=X0.reshape((-1,1))
    MU=constants['MU']
    # function [Xt,Xtorb,X0,X0orb,ft,f0,PHI,p0,pt]=sim1_forwardSTM(t,X0,MU,computePHI,computeProb,probFunc0)

    r0=np.norm(X0[0:3])
    v0=np.norm(X0[3:6])

    r0vec=X0[0:3]
    v0vec=X0[3:6]

    X0orb = XYZ2OE(X0,constants)
    a0,ecc0,E0,w0,inc0,Omega0=X0orb.reshape(1,-1)[0]

    # M0=E0-ecc0*np.sin(E0)
    # X0orb=[a0,ecc0,E0,w0,inc0,Omega0]
    f0 = 2*np.arctan2(np.sqrt(1+ecc0)*np.tan(E0/2),np.sqrt(1-ecc0))

    # n=np.sqrt(MU/a0^3)

    Xt,Xtorb=orbit_propagate_kepler(t,X0,constants)
    a,ecc,E,w,inc,Omega=Xtorb.reshape(1,-1)[0]
    # M=E-ecc*np.sin(E)
    ft = 2*np.arctan2(np.sqrt(1+ecc)*np.tan(E/2),np.sqrt(1-ecc))

    if options['computePHI']==False:
        PHI=np.nan
    else:
        rtvec=Xt[0:3]
        vtvec=Xt[3:6]
        rt=np.norm(rtvec)
        vt=np.norm(vtvec)
        
        
        th=ft-f0
        sig0=np.dot(r0vec.transpose(),v0vec)/np.sqrt(MU)
        p=a*(1-ecc^2)
        r=p*r0/( r0+(p-r0)*np.cos(th)-np.sqrt(p)*sig0*np.sin(th) )
        
        F=1-(r/p)*(1-np.cos(th))
        G=r*r0*np.sin(th)/np.sqrt(MU*p)
        Ft=np.sqrt(MU)/(r0*p)*(sig0*(1-np.cos(th))-np.sqrt(p)*np.sin(th))
        Gt=1-f0/p*(1-np.cos(th))
        
        Xfg=np.vstack((F*X0[0:3]+G*X0[3:6],Ft*X0[0:3]+Gt*X0[3:6]))
        
        # % PHI=[F,G;Ft,Gt]
        t0=0
        x=universalkepler(a,t0,X0,t,options)
        alpha=1/a
        C,S=transcendentalfunc(alpha*x**2)
        
        Q=x^3/alpha*( C-3*S)-np.sqrt(MU)*(t-t0)*x**2*C
        Q=Q/np.sqrt(MU)
        
        
        R=r0/MU*(1-F)*( (rtvec-r0vec)*v0vec.transpose()-(vtvec-v0vec)*r0vec.transpose()  )+Q/MU*vtvec*v0vec.transpose()+G*np.identity(3)
        V=r0/MU*(vtvec-v0vec)*(vtvec-v0vec).transpose()+1/rt^3*( r0*(1-F)*rtvec*r0vec.transpose()-Q*rtvec*v0vec.transpose() )+Gt*np.identity(3)
        Rt=rt/MU*(vtvec-v0vec)*(vtvec-v0vec).transpose()+1/r0^3*(r0*(1-F)*rtvec*r0vec.transpose()+Q*vtvec*r0vec.transpose() )+F*np.identity(3)
        Vt=(-1/r0^2*(vtvec-v0vec)*r0vec.transpose()-1/rt^2*rtvec*(vtvec-v0vec).transpose()+Ft*( np.identity(3) - 1/rt^2*rtvec*rtvec.transpose()+
            1/(MU*rt)*(rtvec*vtvec.transpose()-vtvec*rtvec.transpose())*rtvec*(vtvec-v0vec).transpose())-MU*Q*rtvec*r0vec.transpose()/(rt^3*r0^3) )
        
        
        
        
        PHI=np.block([[Rt,R],[Vt,V]])

    if options['computeProb']==False:
        p0=np.nan
        pt=np.nan
    else:
       p0=probFunc0(X0)
       pt=p0/abs(np.linalg.det(PHI)) 

    return (Xt,Xtorb,X0,X0orb,ft,f0,PHI,p0,pt)