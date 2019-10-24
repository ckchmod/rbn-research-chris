%% BN Sampled Receiving Nodes
% A single cell dynamic, where the receiving node is sampled from a
% probability distribution.
% --------------------------------------

clc; clear all;

%p = genpath('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris');
%addpath (p);

%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat');
load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpOrdered.mat');

topology = 'single'; 
interaction = 1;
numCells = 1;
perturb = .05; 
steps = 100000;

RBNpnew = boolCellGrid(topology, numCells, RBNp.numGenes, RBNp.k, RBNp.p, ...
     interaction, RBNp.initState, RBNp.initTtable, RBNp.initVar, perturb);
RBNpnew.update_all(steps);
RBNp_ssDist = ssDist(RBNpnew);

[lattice_size, state_size] = size(RBNp_ssDist);
lattice_size = 9; % tbd

 subplot(sqrt(lattice_size),4, 11);
 hold on
 s1 = bar(RBNp_ssDist);
 str=sprintf('Single Cell Sample: [0 1]');
 set(s1,'FaceAlpha',0.7);
 title(str);axis([0 1024 0 0.5]);
 set(gca, 'YScale', 'log')
 hold off
 
%  
%  change interaction rule for one cell
% check finding entropy
% check pairwise kullback
% 
% email by wednesdya for meeting seutp 