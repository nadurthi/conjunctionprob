% conjunctuion probability by GMM

MU      = 3.986004418e5;     % Gravitational Const
Re      = 6378.137;          % Earth radius (km)
opt = odeset('reltol',1e-12,'abstol',1e-12);

% satellite 1:
X0_mu = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
P0 =blkdiag(0.01,0.01,0.01,1e-6,1e-6,1e-6);
T=[linspace(0,160690,10),[160691:160790]];

% [X,Xorb]=sim1_forwardmap(T(85),X0_mu,MU);
XX0=X0_mu; %Xmc(i,:)'
[Xt,Xtorb,X0,X0orb,ft,f0,PHI]=sim1_forwardSTM(T(85),XX0,MU);

Xm=Xmc(10,:)';
[t,X1]= ode45(@twoBody,T,Xm,opt);
 [X1(85,:)',Xt,Xt+PHI*(Xm-XX0)]
 
 
%%
% [Xcut8,wcut8]=conjugate_dir_gausspts_till_8moment(mu,P);
Nmc=10000;
Xmc=mvnrnd(X0_mu,P0,Nmc);
wmc=ones(Nmc,1)/Nmc;
XT=zeros(Nmc,6);
for i=1:Nmc
[X,Xorb]=sim1_forwardmap(T(85),Xmc(i,:)',MU);
XT(i,:)=X;
end

plot3(XT(:,1),XT(:,2),XT(:,3),'ro')

figure
plot3(XT(:,4),XT(:,5),XT(:,6),'ro')


% OE0_mu = XYZ2OE_m(X0_mu,MU);
% 
% % OE=[a,e,E,w,i,Om];
% a0=OE0_mu(1);
% ecc0=OE0_mu(2);
% E0=OE0_mu(3);
% w0=OE0_mu(4);
% inc0=OE0_mu(5);
% Omega0=OE0_mu(6);
% M0=E0-ecc*sin(E0);
% [r,v]= elm2rv(a0,ecc0,inc0,Omega0,w0,M0,0,MU);
% Xorb0_mu=[r;v];
% 
% 
% 
% [X0_mu,Xorb]
% 
% n=sqrt(MU/a^3);
% 
% 
% 
% M=n*T(65) + M0;
% E=kepler_rob(M,ecc,1e-12);
% [r,v]= elm2rv(a,ecc,inc,Omega,w,M,0,MU);
% XorbT=[r;v];
% 
% % [X,w]=conjugate_dir_gausspts_till_8moment(mu,P);
% [t,X1]= ode45(@twoBody,T,X0_mu,opt);
% 
% [XorbT,X1(65,:)']
