%% Large scale simulations of RBNp's 
clc; clear all;

topology = 'symmetric'; 
numCells = 10^2; 
numGenes = 10;   
interaction = 1;
perturb = .01; 
steps = 5;

k=1; p=.2; 
threshold = 0.02;

RBNp = boolCellGrid(topology, numCells,numGenes, k, p, ...
    interaction, [], [], [], perturb);
RBNp.update_all(steps);
