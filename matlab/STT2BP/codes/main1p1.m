clear all; close all; clc;
global mue mu1 c1;

a0 = 7100; e0 = 0.2; i0 = 35*pi/180; Om0 = 20*pi/180; w0 = 15*pi/180; 
M0 = 0.00*pi/3; tol = 1e-15; E0g = 1; % Mean anamoly for the Kepler's equation solution
[E0]=kepler(M0, e0, E0g, tol) % E - Eccentric anamoly

re = 6372.797;   % Mean earth radius (Km)
mue = 398601.2;  % Gravitational parameter of earth (Km^3/Sec^2)
% Orbit initial orientation parameters
R = [cos(Om0), -sin(Om0), 0; sin(Om0), cos(Om0), 0; 0 0 1];
R = R*[1, 0, 0; 0, cos(i0), -sin(i0); 0 sin(i0) cos(i0)];
R = R*[cos(w0), -sin(w0), 0; sin(w0), cos(w0), 0; 0 0 1];
% initial conditions in the orbit frame {ie,ib,ih}
r0 = a0*(1-e0*cos(E0));
xb = [a0*(cos(E0)-e0); a0*sqrt(1-e0^2)*sin(E0); 0]; 
xdotb = [-sqrt(mue*a0)*sin(E0)/r0; sqrt(mue*a0*(1-e0^2))*cos(E0)/r0; 0 ];

% initial position and velocity (inertial frame ECEF)
X = R*xb; Xdot = R*xdotb;   xb0 = [X;Xdot];
% orbital period
tp = 2*pi*sqrt(a0^3/mue);
t = linspace(0,2*tp,1000);
tspan = t;

options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(6,1));
[T,Y] = ode45('orbitrect',tspan,xb0,options);
% F&G series solution
[T1,Y1,F,G,Ehatt] = FGsolve(xb0,t);

figure
% plot3(Y(:,1),Y(:,2),Y(:,3),'-',Y1(:,1),Y1(:,2),Y1(:,3),'-.r');
subplot(311)
plot(T,abs(Y(:,1)-Y1(:,1)));
subplot(312)
plot(T,abs(Y(:,2)-Y1(:,2)));
subplot(313)
plot(T,abs(Y(:,3)-Y1(:,3)));

figure
subplot(311)
plot(T,abs(Y(:,4)-Y1(:,4)));
subplot(312)
plot(T,abs(Y(:,5)-Y1(:,5)));
subplot(313)
plot(T,abs(Y(:,6)-Y1(:,6)));

N=6; mu1 = mue; c1 = 0;

t1 = linspace(0,2*tp,1000); tspan1 = t1;
y0a(N^3+N^2+N,1)=0; y0a(1:N,1)=xb0; I6 = eye(6);
y0a(N+1:N+N^2,1)= [I6([1:N^2])]';
options1 = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(N^3+N^2+N,1));
[Tstma,Ystma] = ode45(@STTNint4,tspan1,y0a,options1);

[Tstmfg1,Ystmfg1] = FGSTTs(xb0,t1);
%[Tstmfg1,Ystmfg1] = FGSTT1Battin(xb0,t1);

err = abs(Ystma(:,N+1:N+N^2)-Ystmfg1(:,N+1:N+N^2));

figure
plot(Tstma,err);
title('Error Between Numerical and Analytical Solution: First Order STM');

 
%%% verify
for j = 1:length(t)
    
    PhiN(:,:, j) = reshape(Ystma(j,N+1:N+N^2), N, N);
    Phifg(:,:,j) = reshape(Ystmfg1(j,N+1:N+N^2),N,N);
    
end



