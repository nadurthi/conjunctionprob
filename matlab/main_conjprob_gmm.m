% conjunctuion probability by GMM

MU      = 3.986004418e5;     % Gravitational Const
Re      = 6378.137;          % Earth radius (km)
opt = odeset('reltol',1e-12,'abstol',1e-12);

% satellite 1:
X0_mu = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
P0 =blkdiag(0.01,0.01,0.01,1e-6,1e-6,1e-6);

probFunc0=@(X)mvnpdf(X,X0_mu,P0);

T=[linspace(0,160690,10),[160691:160790]];

% [X,Xorb]=sim1_forwardmap(T(85),X0_mu,MU);
XX0=X0_mu; %Xmc(i,:)'
[Xt,Xtorb,X0,X0orb,ft,f0,PHI]=sim1_forwardSTM(T(85),XX0,MU,true);

Xm=Xmc(10,:)';
[t,X1]= ode45(@twoBody,T,Xm,opt);
 [X1(85,:)',Xt,Xt+PHI*(Xm-XX0)]
 
 
%%
% [Xcut8,wcut8]=conjugate_dir_gausspts_till_8moment(mu,P);
Nmc=10000;
Xmc=mvnrnd(X0_mu,P0,Nmc);
wmc=ones(Nmc,1)/Nmc;
XT=zeros(Nmc,6);
pdf=zeros(Nmc,1);
for i=1:Nmc
    [X,Xorb,X0,X0orb,ft,f0,PHI,p0,pt]=sim1_forwardSTM(T(85),Xmc(i,:)',MU,true,true,probFunc0);
% [X,Xorb]=sim1_forwardmap(T(85),Xmc(i,:)',MU);
XT(i,:)=X;
pdf(i)=pt;
i
end

plot3(XT(:,1),XT(:,2),XT(:,3),'ro')

figure
plot3(XT(:,4),XT(:,5),XT(:,6),'ro')

figure
plot3(XT(:,1),XT(:,2),abs(pdf),'ro')

%%
muT=mean(XT);
PT=cov(XT);
pdfGasssT=mvnpdf(XT,muT,PT);

figure
plot3(XT(:,1),XT(:,2),abs(pdf),'ro')
figure
plot3(XT(:,1),XT(:,2),abs(pdfGasssT),'kx')


