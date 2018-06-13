%--------------------------------------------------------------------------
%                             Main
%--------------------------------------------------------------------------
%% Initialisation
addpath('./script_matlab/');
load TorMod.mat;
d = 1;
nB = length(G1.b)-1-d;
nA = length(G1.f)-1;
%% Simulation
% sim('project.slx');