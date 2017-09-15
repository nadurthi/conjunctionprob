function dy = STTNint4(t,y)
% State Transition Tensor Integration
% General Code
% Computes upto 3rd order STM about a nominal
% convention - 1st order standard STM.
% For verification of the 2nd order STTs, 3rd order is set to zero.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem Specific Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global mu1 c1;
nt = 2; N = 6;
x(1:N,1)=y(1:N,1);   % Nominal Trajectory
N4 = N^4; N3 = N^3; N2 = N^2; N5 = N*N4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% High order partials of the system model
[F]=funvalo(x,nt+1,N);  %N*(1-N^(nt+1))/(1-N) dimensional array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dy(N*(1-N^(nt+1))/(1-N),1)=0;

VD2F= F(N2+N+1:N3+N2+N,1);
%VD3F = F(N3+N2+N+1:N4+N3+N2+N,1);
%VD4F = F(N4+N3+N2+N+1:N5+N4+N3+N2+N,1);

% Zeroth order   ( Nominal Motion)
dy(1:N,1)=F(1:N,1);

% Matricized state transition tensors
    DF = reshape(F(N+1:N2+N,1),N,[]);
    Phi2 = reshape(y(N+1:N2+N,1),N,[]);
    Phi3 = reshape(y(N2+N+1:N3+N2+N,1),N,[]);
%     Phi4 = reshape(y(N3+N2+N+1:N4+N3+N2+N,1),N,[]);
    
% First Order
IN = eye(N); TJac = kron(IN,DF); 
dy(N+1:N2+N,1) = TJac*y(N+1:N2+N,1);

% Second Order
TJac = kron(IN,TJac); Flas = kron(Phi2',kron(Phi2',IN));

dy(N2+N+1:N3+N2+N,1)=TJac*y(N2+N+1:N3+N2+N,1)+Flas*VD2F;


% % Third Order
% TJac = kron(IN,TJac); Flas = kron(Phi2',Flas);
% in = N3+N2+N;
% 
% vect21 = kron(Phi2',kron(Phi3',IN))*VD2F;
% T2 = reshape(vect21,N,N,N,N); 
% T2 = T2 + permute(T2,[1 2 4 3]);
% vect22 = kron(Phi3',kron(Phi2',IN))*VD2F;
% vect2 = T2([1:N4]') + vect22;
% 
% dy(in+1:in+N^4,1)= TJac*y(in+1:in+N^4,1) + vect2 + Flas*VD3F;
% 

t
