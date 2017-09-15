clear all; close all; clc;

%% Semi analytical orbit transfer

AU = 1.49597871*1e8; muS = 1.3271244*1e11;
TU = sqrt(AU^3/muS)
% location in Earth's orbit
Re = 1.0; 
Rm = 2.279391*1e8/AU;
R1 = [Re; 0; 0]; r1 = norm(R1);
R2 = [0; Rm; 0]; r2 = norm(R2);

cbar = R2 - R1;
am = norm(cbar)/2;



