function [dy]=orbitrect(t,y)
dy(6,1) = 0;

mue = 398601.2;  % Gravitational parameter of earth

mu1 = mue;
x = y(1:3,1); xdot = y(4:6,1);

r = sqrt(x'*x); v = sqrt(xdot'*xdot);

dx = xdot;

dxdot(1,1) = -mu1*x(1)/(r^3);
dxdot(2,1) = -mu1*x(2)/(r^3);
dxdot(3,1) = -mu1*x(3)/(r^3);

dy = [dx;dxdot];
