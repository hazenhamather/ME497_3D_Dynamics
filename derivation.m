clear variables
close all
clc

syms psi(t) theta(t) phi(t)

syms Jb1B Jb2B Jb3B Jb1T Jb2T Jb3T mb mt L1 L2 tdata eta etadot gamma
assume((Jb1B>0) & (Jb2B>0) & (Jb3B>0)&(Jb1T>0)&(Jb2T>0)&(Jb3T>0)&...
    (mb>0)&(mt>0)&(L1>0)&(L2>0));

%Moments of inertia of the Body
% Jb1B = 7.37e-7; %kg-m^2
% Jb2B = 5.42e-8;
% Jb3B = 7.81e-7;

%Moments of inertia of the Tail
% Jb1T = 2.34e-8; %kg-m^2
% Jb2T = 2.34e-8;
% Jb3T = 8.86e-10;

%Masses and Length
% mb = 3.04e-3; %mass of body, kg
% L1 = 2.53e-2; %Length of body, m
% mt = 2.10e-4; %mass of tail, kg
% L2 = 1.35e-2; %Length of tail, m

%Body and Tail inertia matrix
JB = [Jb1B,0,0;0,Jb2B,0;0,0,Jb3B];
JT = [Jb1T,0,0;0,Jb2T,0;0,0,Jb3T];

%Body Rotation Matrices
RB1 = [cos(psi), sin(psi), 0; -sin(psi), cos(psi), 0;0,0,1];
RB2 = [1,0,0;0,cos(theta), sin(theta); 0, -sin(theta), cos(theta)];
RB3 = [cos(phi), 0, -sin(phi);0,1,0;sin(phi), 0, cos(phi)];

%Omega in the body fixed frame
omegab = (RB3*RB2*RB1)*[0;0;diff(psi)] + (RB3*RB2)*[diff(theta);0;0] + ...
    RB3*[0;diff(phi);0];

%Tail Rotation Matrices, 
RT1 = [cos(eta), 0, -sin(eta);0,1,0;sin(eta), 0, cos(eta)];
RT2 = [1,0,0;0,cos(gamma), sin(gamma);0,-sin(gamma), cos(gamma)];

%Omega in the tail fixed frame
temp = omegab + [0;etadot;0];
omegat = (RT2*RT1)*temp;

%Angular Momentum of Body
LB = JB*omegab; %about the body center
LBT = (RT2*RT1) * LB; %in the tail axes

%Angular Momentum of Tail
LT = JT*omegat;

%Angular Momentum about G
rtb = [0;0;L2] + (RT2*RT1)*[0;L1;0]; %in tail centererd axes
vtb = L1*(RT2*RT1)*(cross(omegab,[0;1;0])) + L2*cross(omegat, [0;0;1]);
LsysG = LBT+LT + ((mb*mt)/(mb+mt))*(cross(rtb,vtb));

%Create M and f files
omegabint = (RT2*RT1)*omegab; %moves omegab to t-frame
ODEsRot = LsysG == 0;
stateEqs = simplify(ODEsRot);
stateVars = [psi;theta;phi];
[Msym, fsym] = massMatrixForm(stateEqs, stateVars);
Msym = simplify(Msym);
fsym = simplify(fsym);
M = odeFunction(Msym, stateVars,Jb1B, Jb2B, Jb3B, Jb1T, Jb2T, Jb3T, mb, mt, L1, L2, tdata, eta, etadot, gamma,'file','M');
f = odeFunction(fsym, stateVars,Jb1B, Jb2B, Jb3B, Jb1T, Jb2T, Jb3T, mb, mt, L1, L2, tdata, eta, etadot, gamma,'file','f');