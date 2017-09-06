function [p0,pt,Xt]=sim1_forwardProb(t,X0,MU,probFunc0)
% [Xt,Xtorb]=sim1_forwardmap(t,X0,MU);
% PHI=sim1_forwardSTM(t,X0,MU);
[Xt,Xtorb,X0,X0orb,ft,f0,PHI]=sim1_forwardSTM(t,X0,MU,true);

p0=probFunc0(X0);
pt=p0/det(PHI);
