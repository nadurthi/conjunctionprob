function [Xt,Xtorb,X0,X0orb,ft,f0,PHI]=sim1_forwardSTM(t,X0,MU)
r0=norm(X0(1:3));
v0=norm(X0(4:6));

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
X0orb=[a0,ecc0,E0,w0,inc0,Omega0];
f0 = 2*atan2(sqrt(1+ecc0)*tan(E0/2),sqrt(1-ecc0));



[Xt,Xtorb]=sim1_forwardmap(t,X0,MU);
a=Xtorb(1);
ecc=Xtorb(2);
E=Xtorb(3);
w=Xtorb(4);
inc=Xtorb(5);
Omega=Xtorb(6);
M=E-ecc*sin(E);
ft = 2*atan2(sqrt(1+ecc)*tan(E/2),sqrt(1-ecc));

th=ft-f0;
sig0=X0(1:3)'*X0(4:6)/sqrt(MU);
p=a*(1-ecc^2);
r=p*r0/( r0+(p-r0)*cos(th)-sqrt(p)*sig0*sin(th) );

F=1-(r/p)*(1-cos(th));
G=r*r0*sin(th)/sqrt(MU*p);
Ft=sqrt(MU)/(r0*p)*(sig0*(1-cos(th))-sqrt(p)*sin(th));
Gt=1-f0/p*(1-cos(th));


PHI=[F,G;Ft,Gt];


% 
% [E0,f0] = kepler_rob(M0,ecc0,1e-12);


end

