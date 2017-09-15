function  [T,Y] = FGSTTs(x0,t)
% F and G series solution - analytic solution to the 2 Body Problem
% Analytic computation of the state transition tensors

mue = 398601.2;  % Gravitational parameter of earth
tol = 1e-15; N = 6;
R0 = x0(1:3,:); V0 = x0(4:6,:);
r0 = sqrt(R0'*R0); v0 = sqrt(V0'*V0);
a = 1/(2/r0 - v0^2/mue);
sig0 = R0'*V0/sqrt(mue);
I6 = eye(6); Y(1,1:N+N^2+N^3)=0;
Y(1,1:N)=x0'; T = t;
Y(1,N+1:N+N^2) = [I6([1:N^2])]';
Ehat = 0; I3 = eye(3);

for k = 2:length(t)
    dt = t(k)-t(1);
    Ehat0 = Ehat;
    [Ehat]=keplermod(dt, r0, a, sig0, Ehat0, tol);
    Ehatt(k,1)=Ehat;
    
    r = a + (r0 - a)*cos(Ehat) + sqrt(a)*sig0*sin(Ehat);
    % F and G functions
    F(k,1) = 1-a*(1-cos(Ehat))/r0;
    G(k,1) = dt + sqrt(a^3/mue)*(sin(Ehat)-Ehat);
    % Ft and Gt functions
    F(k,2) = -sqrt(mue*a)*sin(Ehat)/(r*r0);
    G(k,2) = 1+a*(cos(Ehat)-1)/r;
    
    Rt = F(k,1)*R0 + G(k,1)*V0;
    Vt = F(k,2)*R0 + G(k,2)*V0;
    Y(k,1:3) = Rt'; 
    Y(k,4:6) = Vt';
    
    %%% State Transition Tensors
    %% First Order
    dr0dR0 = R0/r0; dr0dV0 = 0*R0; 
    dv0dR0 = 0*V0; dv0dV0 = V0/v0;
    dadR0 = 2*a^2*R0/(r0^3); dadV0 = 2*a^2*V0/(mue);
    dsig0dR0 = V0/sqrt(mue); dsig0dV0 = R0/sqrt(mue);
    Ea = (sig0*(1-cos(Ehat))/(2*r*sqrt(a))-1.5*sqrt(mue/a^3)*dt/r +...
        r0*sin(Ehat)/(a*r));
    Er0 = -sin(Ehat)/r; 
    Esig0 = -sqrt(a)*(1-cos(Ehat))/r;
    dEhatdR0 = Ea*dadR0 + Er0*dr0dR0 + Esig0*dsig0dR0;
    dEhatdV0 = Ea*dadV0 + Er0*dr0dV0 + Esig0*dsig0dV0;
    rhoa = 1 - cos(Ehat) + sig0*sin(Ehat)/(2*sqrt(a));
    rhoE = sig0*sqrt(a)*cos(Ehat) - (r0-a)*sin(Ehat);
    rhosig0 = sqrt(a)*sin(Ehat);
    rhor0 = cos(Ehat);
    drdR0 = rhoa*dadR0+rhoE*dEhatdR0+rhor0*dr0dR0+rhosig0*dsig0dR0;
    drdV0 = rhoa*dadV0+rhoE*dEhatdV0+rhor0*dr0dV0+rhosig0*dsig0dV0;
    %%%%% F and G general sensitivity coefficients
    sinE = sin(Ehat); cosE=cos(Ehat);
    Fa=(cosE-1)/r0; Fr0=a*(1-cosE)/(r0^2); FE=-a*sinE/r0;
    dFdR0 = Fa*dadR0 + Fr0*dr0dR0 + FE*dEhatdR0;
    dFdV0 = Fa*dadV0 + Fr0*dr0dV0 + FE*dEhatdV0;
    
    Ga=1.5*sqrt(a/mue)*(sinE-Ehat);GE=sqrt(a^3/mue)*(cosE-1);
    dGdR0 = Ga*dadR0 +  GE*dEhatdR0;
    dGdV0 = Ga*dadV0 +  GE*dEhatdV0;
    
    Fta=-sqrt(mue/a)*sinE/(2*r*r0);Ftr=sqrt(mue*a)*sinE/(r0*r^2);
    Ftr0=sqrt(mue*a)*sinE/(r*r0^2);FtE=-sqrt(mue*a)*cosE/(r0*r);
    dFtdR0 = Fta*dadR0 + Ftr*drdR0 + Ftr0*dr0dR0 + FtE*dEhatdR0;
    dFtdV0 = Fta*dadV0 + Ftr*drdV0 + Ftr0*dr0dV0 + FtE*dEhatdV0;
    
    Gta = -(1-cosE)/r; Gtr=a*(1-cosE)/r^2; GtE=-a*sinE/r;
    dGtdR0 = Gta*dadR0 + Gtr*drdR0 + GtE*dEhatdR0;
    dGtdV0 = Gta*dadV0 + Gtr*drdV0 + GtE*dEhatdV0;
    
    Phi11 = F(k,1)*I3 + R0*dFdR0' + V0*dGdR0';
    Phi21 = F(k,2)*I3 + R0*dFtdR0' + V0*dGtdR0';
    Phi12 = G(k,1)*I3 + R0*dFdV0' + V0*dGdV0';
    Phi22 = G(k,2)*I3 + R0*dFtdV0' + V0*dGtdV0';
    Phit = [Phi11, Phi12; Phi21, Phi22];
    size(Phit)
    Y(k,N+1:N+N^2) = reshape(Phit,1,N^2);
    
    %%% Second Order State Transition Tensor
%     d2r0d2R0=I3/r0-R0*R0'/r0^3; d2v0d2V0=I3/v0-V0*V0'/v0^3;
%     d2r0dR0dV0 = 0*I3; d2r0dV0dR0 = 0*I3;
%     d2sig0dR0dV0=I3/sqrt(mue); d2sig0dV0dR0=I3/sqrt(mue);
%     d2sig0d2R0=0*I3; d2sig0d2V0=0*I3;
%     d2ad2R0 = 2*a^2*(I3 + (4*a/r0^3 - 3/r0^2)*R0*R0')/r0^3;
%     d2adR0dV0 = 8*a^3*R0*V0'/(mue*r0^3); d2adV0R0 = d2adR0V0';
%     d2ad2V0 = 2*(a^2)*I3/mue + 8*(a^3)*V0*V0'/(mue^2);
%     
%     Ea_r = -Ea/r; Ea_r0 = sinE/(r*a); 
%     Ea_E = (r0*cosE/(r*a)+sig0*sinE/(2*sqrt(a)*r));
%     Ea_a = -(r0*sinE/(r*a^2)+sig0*(1-cosE)/(4*r*sqrt(a^3))- ...
%         9*sqrt(mue/a^5)*dt/(4*r));
%     Ea_sig0 = (1-cosE)/(2*sqrt(a)*r);
%     
%     dEadR0 = Ea_r*drdR0 + Ea_E*dEhatdR0 + Ea_a*dadR0 ...
%                 + Ea_sig0*dsig0dR0 + Ea_r0*dr0dR0;
%     dEadV0 = Ea_r*drdV0 + Ea_E*dEhatdV0 + Ea_a*dadV0 ...
%                 + Ea_sig0*dsig0dV0 + Ea_r0*dr0dV0;
%             
%     Esig0_a = Esig0/(2*a); Esig0_r=-Esig0/r; Esig0_E=-sqrt(a)*sinE/r;
%     
%     dEsig0dR0 = Esig0_a*dadR0+Esig0_r*drdR0+Esig0_E*dEhatdR0;
%     dEsig0dV0 = Esig0_a*dadV0+Esig0_r*drdV0+Esig0_E*dEhatdV0;
%     
%     Er0_r= -Er0/r; Er0_E=-cosE/r; 
%     dEr0dR0 = Er0_r*drdR0 + Er0_E*dEhatdR0;
%     dEr0dV0 = Er0_r*drdV0 + Er0_E*dEhatdV0;
%     
%     rhoa_a = -(rhoa + sig0*sinE/(4*sqrt(a)))/a; rhoa_r = 1/a; 
%     rhoa_r0 = -cosE/a; rhoa_sig0=-sinE/(2*sqrt(a)); 
%     rhoa_E = (r0*sinE/a - sig0*cosE/(2*sqrt(a)));
%     drhoadR0 = rhoa_a*dadR0+rhoa_r*drdR0+rhoa_r0*dr0dR0+...
%                rhoa_sig0*dsig0dR0+rhoa_E*dEhatdR0;
%     drhoadV0 = rhoa_a*dadV0+rhoa_r*drdV0+rhoa_r0*dr0dV0+...
%                rhoa_sig0*dsig0dV0+rhoa_E*dEhatdV0;
%     
%     rhoE_E = -(sig0*sqrt(a)*sinE + (r0-a)*cosE); rhoE_sig0 = sqrt(a)*cosE;
%     rhoE_r0 = -sinE; rhoE_a=(sig0*cosE/sqrt(a)+sinE);
%     drhoEdR0 = rhoE_E*dEhatdR0+rhoE_sig0*dsig0dR0;
end

