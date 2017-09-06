function x=universalkepler(a,t0,X0,t,MU,tol)
X0=X0(:);
r0=X0(1:3);
v0=X0(4:6);
normr0=norm(r0);

nn=5;
sig0=r0'*v0/sqrt(MU);
alpha=1/a;
% err=1000;
% while err>tol
%     [C,S]=transcendentalfunc(alpha*x^2);
%     F=sig0*x^2*C+(1-normr0*alpha)*x^3*S+normr0*x-sqrt(MU)*(t-t0);
%     
%     dF=sig0*x*(1-alpha*x^2*S)+(1-normr0*alpha)*x^2*C+normr0;
%     
%     n*F/( dF+sign(dF)*abs( (n-1)^2*dF^2-n*(n-1)*F   )  )
%     
% end
% rp=a-a*e;
% x0=MU*(t-t0)^2/( rp*  )
% keyboard
opt=optimset('TolX',1e-12,'TolFun',1e-12);
x=fsolve(@(x)F(x,MU,alpha,sig0,normr0,t,t0),0,opt);



end


function err=F(x,MU,alpha,sig0,normr0,t,t0)

 [C,S]=transcendentalfunc(alpha*x^2);
 err=sig0*x^2*C+(1-normr0*alpha)*x^3*S+normr0*x-sqrt(MU)*(t-t0);


end