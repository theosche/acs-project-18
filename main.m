%--------------------------------------------------------------------------
%                             Main
%--------------------------------------------------------------------------
%% Initialisation
clear all; clc
addpath('./script_matlab/');
load TorMod.mat;
d = 1;
nB = length(G1.b)-1-d;
nA = length(G1.f)-1;
Ts = G1.Ts;
% load controller
load('RST/RSTG1.mat');
load('RST/RSTG2.mat');
load('RST/RSTG3.mat');
% Polynomial for controller 4
load('RST/polynomial.mat');
% Hr = 0;
% Hs = 0;
Hr = [1 1];
Hs = [1 -1];
nHs = length(Hs)-1; nHr = length(Hr)-1;
nr = nA + nHs + nHr -1;
ns = nB + nHs + nHr +d -1;
nt = 0;
%% Suppervisor parameters
beta = 50;
lambda = 0.02;
DT = 500;
%% Simulation
% sim('project.slx');