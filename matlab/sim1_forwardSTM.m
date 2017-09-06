function [Xt,Xtorb,X0,X0orb,ft,f0,PHI,p0,pt]=sim1_forwardSTM(t,X0,MU,computePHI,computeProb,probFunc0)
r0=norm(X0(1:3));
v0=norm(X0(4:6));

r0vec=X0(1:3);
v0vec=X0(4:6);

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

n=sqrt(MU/a0^3);



[Xt,Xtorb]=sim1_forwardmap(t,X0,MU);
a=Xtorb(1);
ecc=Xtorb(2);
E=Xtorb(3);
w=Xtorb(4);
inc=Xtorb(5);
Omega=Xtorb(6);
M=E-ecc*sin(E);
ft = 2*atan2(sqrt(1+ecc)*tan(E/2),sqrt(1-ecc));

if computePHI==0
    PHI=NaN;
else
    rtvec=Xt(1:3);
    vtvec=Xt(4:6);
    rt=norm(rtvec);
    vt=norm(vtvec);
    
    
    th=ft-f0;
    sig0=X0(1:3)'*X0(4:6)/sqrt(MU);
    p=a*(1-ecc^2);
    r=p*r0/( r0+(p-r0)*cos(th)-sqrt(p)*sig0*sin(th) );
    
    F=1-(r/p)*(1-cos(th));
    G=r*r0*sin(th)/sqrt(MU*p);
    Ft=sqrt(MU)/(r0*p)*(sig0*(1-cos(th))-sqrt(p)*sin(th));
    Gt=1-f0/p*(1-cos(th));
    
    % keyboard
    Xfg=[F*X0(1:3)+G*X0(4:6); Ft*X0(1:3)+Gt*X0(4:6)];
    
    % PHI=[F,G;Ft,Gt];
    t0=0;
    x=universalkepler(a,t0,X0,t,MU,1e-12);
    alpha=1/a;
    [C,S]=transcendentalfunc(alpha*x^2);
    
    % t0=M0/n;
    Q=x^3/alpha*( C-3*S)-sqrt(MU)*(t-t0)*x^2*C;
    Q=Q/sqrt(MU);
    
    
    R=r0/MU*(1-F)*( (rtvec-r0vec)*v0vec'-(vtvec-v0vec)*r0vec'  )+Q/MU*vtvec*v0vec'+G*eye(3);
    V=r0/MU*(vtvec-v0vec)*(vtvec-v0vec)'+1/rt^3*( r0*(1-F)*rtvec*r0vec'-Q*rtvec*v0vec' )+Gt*eye(3);
    Rt=rt/MU*(vtvec-v0vec)*(vtvec-v0vec)'+1/r0^3*(r0*(1-F)*rtvec*r0vec'+Q*vtvec*r0vec' )+F*eye(3);
    Vt=-1/r0^2*(vtvec-v0vec)*r0vec'-1/rt^2*rtvec*(vtvec-v0vec)'+Ft*( eye(3) - 1/rt^2*rtvec*rtvec'+...
        1/(MU*rt)*(rtvec*vtvec'-vtvec*rtvec')*rtvec*(vtvec-v0vec)')-MU*Q*rtvec*r0vec'/(rt^3*r0^3);
    
    
    
    %         keyboard
    % rqx=[(1-x^2/r0*C)*X0(1:3)+(t-t0-x^3/sqrt(MU)*S)*X0(4:6); x*sqrt(MU)/(rt*r0)*(alpha*x^2*S-1)*X0(1:3)+(1-x^2/rt*C)*X0(4:6)]
    
    PHI=[Rt,R;Vt,V];
    % PHI=[V',-R';-Vt',Rt'];
    
    % [Xt,Xfg,rqx,inv(PHI)*X0]
    %
    % [E0,f0] = kepler_rob(M0,ecc0,1e-12);
end

if computeProb==false
    p0=NaN;
    pt=NaN;
else
   p0=probFunc0(X0);
   pt=p0/abs(det(PHI)); 
    
end

end

