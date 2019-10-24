%% To do:
% find grn for cell motiliy
% rbn for cell motility
% machine learning boolean network
% 
% --------------------------------------
% periodic boundary conditions
% 
% Keep Receive constant in alternating pattern

% 5-13-19
% Adaption of the environment 
% transition of sampling probabiliyt on receiving node
% M entropy vs order parameter (p)


clc; clear all;

p = genpath('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris');
addpath (p);

load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat');
%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpChaotic.mat');


%inNodes = RBNp.allStates(:,1,:); %First column is all in Nodes.
%outNodes = RBNp.allStates(:,end,:); %last column is all end Nodes.

topology = 'symmetric'; 
numCells = 10^2; 
numGenes = 10;   
interaction = 1;
perturb = .01; 
steps = 3000;

threshold = 0.02;

RBNp = boolCellGrid(topology, numCells,numGenes, RBNp.k, RBNp.p, ...
                interaction, RBNp.initStates, RBNp.initTtable, RBNp.initvarF, perturb);
RBNp.update_all(steps);

RBNp_ssDist = ssDist(RBNp);
figure(1);
[lattice_size, state_size] = size(RBNp_ssDist);
for i = 1:lattice_size
     %j=i;
     subplot(sqrt(lattice_size),sqrt(lattice_size), i);
     hold on
     set(gca, 'YScale', 'log')
     s1 = stem(1:length(RBNp_ssDist(i,:)), RBNp_ssDist(i,:));
     str=sprintf('Cell Number: %d', i);
     title(str);axis([0 1024 0 0.04]);
     hold off
end