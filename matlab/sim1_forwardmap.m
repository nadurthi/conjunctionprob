function [X,Xorb]=sim1_forwardmap(t,X0,MU)

OE0_mu = XYZ2OE_m(X0,MU);
% [p,a,ecc,inc,Omega,argp,nu,m,arglat,truelon,lonper ] = rv2coe (X0_mu(1:3),X0_mu(4:6), MU)

% OE=[a,e,E,w,i,Om];
a0=OE0_mu(1);
ecc0=OE0_mu(2);
E0=OE0_mu(3);
w0=OE0_mu(4);
inc0=OE0_mu(5);
Omega0=OE0_mu(6);
M0=E0-ecc0*sin(E0);
[r,v]= elm2rv(a0,ecc0,inc0,Omega0,w0,M0,0,MU);
Xorb0_mu=[r;v];

n=sqrt(MU/a0^3);

M=n*t + M0;
E=kepler_rob(M,ecc0,1e-12);
Xorb=[a0,ecc0,E,w0,inc0,Omega0];
[r,v]= elm2rv(a0,ecc0,inc0,Omega0,w0,M,0,MU);
X=[r;v];
