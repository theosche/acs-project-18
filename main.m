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
Hs = [1,-1]; 
Hr = [1, 1]; 
% load controller
load('RST/RSTG1.mat');
load('RST/RSTG2.mat');
load('RST/RSTG3.mat');
% Polynomial for controller 4
load('RST/polynomial.mat');
% Suppervisor parameters
beta = 1;
lambda = 0.5;
DT = 10;
%% Simulation
% sim('project.slx');