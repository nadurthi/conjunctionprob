x0=Xmc(200,:);

[Xcut,wcut]=conjugate_dir_gausspts_till_8moment(x0',P0/10000);
Ycut=zeros(size(Xcut));

for i=1:length(wcut)
    [X,Xorb,X0,X0orb,ft,f0,PHI,p0,pt]=sim1_forwardSTM(T(85),Xcut(i,:)',MU,false,true,probFunc0);
    Ycut(i,:)=X;
end
[mx,Px]=MeanCov(Xcut,wcut);
[my,Py]=MeanCov(Ycut,wcut);

Pxy=CrossCov(Xcut,mx,Ycut,my,wcut);
A=Pxy'*inv(Px);
b=my-A*mx;
det(A)
%%