% conjunctuion probability by GMM

MU      = 3.986004418e5;     % Gravitational Const
Re      = 6378.137;          % Earth radius (km)
opt = odeset('reltol',1e-12,'abstol',1e-12);

% satellite 1:
X0_mu = [7000 0 0 1.0374090357 -1.0374090357 7.4771288355]';
P0 =blkdiag(0.000001,0.000001,0.000001,1e-10,1e-10,1e-10);

probFunc0=@(X)mvnpdf(X,X0_mu,P0);

T=[linspace(0,160690,10),[160691:160790]];

% [X,Xorb]=sim1_forwardmap(T(85),X0_mu,MU);
XX0=X0_mu; %Xmc(i,:)'
[Xt,Xtorb,X0,X0orb,ft,f0,PHI]=sim1_forwardSTM(T(85),XX0,MU,false,true,probFunc0);
% [T,Y] = FGSTTs(XX0,[0,T(85)]);
% reshape(Y(end,7:7+35),6,6)
% Xm=Xmc(10,:)';
% [t,X1]= ode45(@twoBody,T,Xm,opt);
%  [X1(85,:)',Xt,Xt+PHI*(Xm-XX0)]
 
 
%%
% [Xcut8,wcut8]=conjugate_dir_gausspts_till_8moment(X0_mu,P0);
Nmc=100000;
Xmc=mvnrnd(X0_mu,P0,Nmc);
% Xmc=[mvnrnd(X0_mu,0.5*P0,1500);mvnrnd(X0_mu,1*P0,1500);mvnrnd(X0_mu,2*P0,1500);mvnrnd(X0_mu,3*P0,1500)];
pdf0=mvnpdf(Xmc,X0_mu',P0);

% for i=[3]
%     X=getptsonsphere(X0_mu,i^6*P0,1);
%     Xmc=vertcat(Xmc,X);
% end
% for i=[3]
%     X=getptsonsphere(X0_mu,i^6*P0,6);
%     Xmc=vertcat(Xmc,X);
% end

Nmc=size(Xmc,1);
wmc=ones(Nmc,1)/Nmc;
XT=zeros(Nmc,6);
pdf=zeros(Nmc,1);
for i=1:Nmc
    [X,Xorb,X0,X0orb,ft,f0,PHI,p0,pt]=sim1_forwardSTM(T(85),Xmc(i,:)',MU,false,true,probFunc0);
% [X,Xorb]=sim1_forwardmap(T(85),Xmc(i,:)',MU);
XT(i,:)=X;
pdf(i)=pt;
i
end

% Ncut=length(wcut8);
% XTcut=zeros(Ncut,6);
% pdfcut=zeros(Ncut,1);
% for i=1:Ncut
%     [X,Xorb,X0,X0orb,ft,f0,PHI,p0,pt]=sim1_forwardSTM(T(85),Xcut8(i,:)',MU,false,true,probFunc0);
% % [X,Xorb]=sim1_forwardmap(T(85),Xmc(i,:)',MU);
% XTcut(i,:)=X;
% pdfcut(i)=pt;
% i
% end

%%
IDX = kmeans(XT,20);
C=['r','b','g','k','y','m','c'];
M=['o','+','s','*','d','^'];
colmark=cell(1,length(C)*length(M));
k=1;
for j=1:length(M)
    for i=1:length(C)
        colmark{k}=strcat(C(i),M(j));
        k=k+1;
    end
end

close all
figure
hold on
for id =unique(IDX)'
    plot3(XT(IDX==id,1),XT(IDX==id,2),XT(IDX==id,3),colmark{id})
end

figure
hold on
for id =unique(IDX)'
    plot3(XT(IDX==id,4),XT(IDX==id,5),XT(IDX==id,6),colmark{id})
end

figure
plot3(XT(:,4),XT(:,5),XT(:,6),'ro')

figure
plot3(XT(:,1),XT(:,2),pdf,'ro')


[M,P]=MeanCov(XT,pdf);
Y=zeros(size(Xmc));
A=inv(sqrtm(P));
for i=1:Nmc
    Y(i,:)=A*(XT(i,:)'-M);
end

figure
plot3(Y(:,1),Y(:,2),Y(:,3),'ro')

figure
hold on
for id =unique(IDX)'
    plot3(Y(IDX==id,1),Y(IDX==id,2),Y(IDX==id,3),colmark{id})
end

IDY = kmeans(Y,25);

figure
hold on
for id =unique(IDY)'
    plot3(Y(IDY==id,1),Y(IDY==id,2),Y(IDY==id,3),colmark{id})
end

figure
hold on
for id =unique(IDX)'
    plot3(XT(IDY==id,1),XT(IDY==id,2),XT(IDY==id,3),colmark{id})
end


%%
muT=mean(XT);
PT=cov(XT);
pdfGasssT=mvnpdf(XT,muT,PT);

figure
plot(XT(:,1),pdf,'bo')
figure
plot(XT(:,1),pdfGasssT,'bo')

figure
plot(XT(:,2),pdf,'bo')
figure
plot(XT(:,2),pdfGasssT,'bo')

figure
plot(XT(:,3),pdf,'bo')
figure
plot(XT(:,3),pdfGasssT,'bo')


figure
plot(XT(:,4),pdf,'bo')
figure
plot(XT(:,4),pdfGasssT,'bo')


figure
plot(XT(:,5),pdf,'bo')
figure
plot(XT(:,5),pdfGasssT,'bo')


figure
plot(XT(:,6),pdf,'bo')
figure
plot(XT(:,6),pdfGasssT,'bo')
