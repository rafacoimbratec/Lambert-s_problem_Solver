% main.m
clear; clc;
mu = 398600.4418; % km^3/s^2
r1 = [7000;0;0];
th = deg2rad(60);
Rz = [cos(th) -sin(th) 0; sin(th) cos(th) 0; 0 0 1];
r2 = Rz*r1;
dt = 1800;

[v1,v2,info] = lambert_uv(r1,r2,dt,mu,'verbose',true);

% Optional: validate with your Kepler propagator
% [r2p,~] = kepler_universal(r1,v1,dt,mu);
% norm(r2p - r2)

% Show results
disp('Universal Variable solution v1:'); disp(v1');

% Optional plot
%plot_transfer(r1,r2,v1_uv,mu,dt);

