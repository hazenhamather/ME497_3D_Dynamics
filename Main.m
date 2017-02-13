clear variables
close all
clc

%Moments of inertia of the Body
Jb1B = 7.37e-7; %kg-m^2
Jb2B = 5.42e-8;
Jb3B = 7.81e-7;

%Moments of inertia of the Tail
Jb1T = 2.34e-8; %kg-m^2
Jb2T = 2.34e-8;
Jb3T = 8.86e-10;

%Masses and Length
mb = 3.04e-3; %mass of body, kg
L1 = 2.53e-2; %Length of body, m
mt = 2.10e-4; %mass of tail, kg
L2 = 1.35e-2; %Length of tail, m

%Define tdata for 120ms
tdata = 0:0.0001:0.120;

%Period
T = 0.120;

%Eta, etadot, and gamma
eta = pi*(1-cos(pi*tdata/T));
etadot = (pi^2/T)*sin(pi*tdata/T);
gamma = 0;

%Initial Conditions
ICs = [0;0;-pi];

%Using ODE45
tol = 1e-6;
options = odeset('mass', @M, 'abstol',tol,'reltol',tol);
[t,y] = ode45(@f,tdata,ICs,options,Jb1B, Jb2B, Jb3B, Jb1T, Jb2T, Jb3T, mb, mt, L1, L2, tdata, eta, etadot, gamma);

tms = t*1000;
psi = y(:,1) * 180/pi;