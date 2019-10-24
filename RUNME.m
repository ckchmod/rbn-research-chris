%% Instructions on how to run RBNp with examples of Analysis
clc; clear all;

%% Set up a tissue of networks with same parameters as in Villani Paper
% Each network is a Random Boolean Network (Perturbed). It has noise, so
% the overall network is ergodic. the function ssDist() calculuates the
% steady state probability distribution for all the cells in the network.

topology = 'hexagonal';
numCells = 4^2; % 20^2
numGenes = 10;   % 100
k = 4;          % 3 
p = 0.21;
interaction = 1; % number of interacting cells (bandwidth)

% Tissue of Boolean Networks with no Perturbation
villaniOrig = boolCellGrid(topology, numCells, numGenes, k, p, ...
    interaction, [], [], []);

% Copying villaniOrig initial conditions to create villaniOrigP
initTtable = villaniOrig.initTtable;
initvarF = villaniOrig.initvarF;
initStates = villaniOrig.initStates;
villaniOrig.update_all(60);
goodStates = villaniOrig.allStates;

% Tissue of Boolean Networks with Perturbation
perturbation = .1; 
villaniOrigP = boolCellGrid(topology, numCells, numGenes, k, p, ...
    interaction, initStates, initTtable, initvarF, perturbation);
villaniOrigP.update_all(60);
badStates = villaniOrigP.allStates;

%==========================================================================
%% Example Analysis
% Create a matrix of pairwise KLD symmetric.

% Get the Steady State Distribution of object
badState_ssDist = ssDist(villaniOrigP);

% Get the pairwise Kulle-Leiber Distances (Symetric)
KLDMatrix = KLDPairwise(badState_ssDist);

% Kull-Leiber Distance Distribution Visualization
histogram(KLDMatrix)
xlabel('KLD Distance'); ylabel('count');
title(['KLD Distance Distribution - numCells: ' num2str(numCells) ...
    ' numGenes: ' num2str(numGenes) ' k: ' num2str(k) ' p: ' num2str(p) ...
    ' bandwidth: ' num2str(interaction) ' criticality: ' ...
    villaniOrigP.criticality]);