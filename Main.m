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
eta = pi*(1-cos(pi*tdata/T)); %radians
etadot = (pi^2/T)*sin(pi*tdata/T); %radians/s
gamma = -20*180/pi; %gecko tail not perpendicular with body
% gamma = 0; %Gecko tail perpendicular to body

%Initial Conditions
ICs = [0;0;-pi];

%Using ODE45
tol = 1e-6;
options = odeset('mass', @M, 'abstol',tol,'reltol',tol);
[t,y] = ode45(@f,tdata,ICs,options,Jb1B, Jb2B, Jb3B, Jb1T, Jb2T, Jb3T, mb, mt, L1, L2, tdata, eta, etadot, gamma);
tms = t*1000;

%Extracting solutions
psi = y(:,1)*180/pi; %degrees
theta = y(:,2)*180/pi; %degrees
phi = y(:,3)*180/pi; %degrees
% eta = eta*180/pi;

%First subplot showing our Euler angles
figure(1)
subplot(3,1,1);
plot(tms,psi);
axis([0 120 -5 5]);
xlabel('Time (ms)');
ylabel('\psi (deg)');

subplot(3,1,2);
plot(tms,theta);
axis([0 120 0 10]);
xlabel('Time (ms)');
ylabel('\theta (deg)');

subplot(3,1,3);
plot(tms,phi);
axis([0 120 -400 -100]);
xlabel('Time (ms)');
ylabel('\phi (deg)');

%Animation Function (convert back to radians)
psi = psi*pi/180;
theta = theta*pi/180;
phi = phi*pi/180;
animate_gecko(psi,theta,phi,eta,gamma,'default');