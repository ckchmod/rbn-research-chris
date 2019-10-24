%clear all; clc;
 topology = 'symmetric'; 
 numCells = 1^2; 
 numGenes = 5;   
 interaction = 0;
 k = 2;
 p = .5;
 perturb = .1;
% RBNp = boolCellGrid(topology, numCells,numGenes, k, p, ...
%     interaction, [], [], [], perturb);
% 
% lambda = log(mean(sum(bnActivity(RBNp.initTtable),2)))
%F = RBNp.initTtable; varF = RBNp.initvarF; nv = [1, 1, 1, 1, 1];
[A, Avec] = bnAsparse(F, varF, nv); attractors = bnAttractor(Avec);
steps = 100;
RBNp = boolCellGrid(topology, numCells, numGenes, k, p, ...
     interaction, [0,0,0,0,0], F, varF, perturb);
RBNp.update_all(steps);
lastRBNpStates = RBNp.allStates(:,:,end);                   
RBNpdummy = boolCellGrid(topology, numCells,numGenes, k, p, ...
    interaction, lastRBNpStates, RBNp.initTtable, RBNp.initvarF, perturb);
delta = size(RBNp.allStates,3)-1;
RBNpdummy.update_all(delta);
lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end);

%Calculate Steady State Distribution Here;
RBNp_ssDist = ssDist(RBNp); %nnz(RBNp_ssDist); 
RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
RBNpstar_ssDist = ssDist(RBNp);%nnz(RBNpstar_ssDist);

% Difference of two distributions:
KLDoutput = .5*(KLD(RBNp_ssDist, RBNpstar_ssDist) + KLD(RBNpstar_ssDist, RBNp_ssDist))

figure;
hold on
scatter(1:32, RBNp_ssDist);
scatter(1:32, RBNpstar_ssDist);
str = sprintf('g = 1.0, p = 0.1, t=100, Difference of Distrib: %f', KLDoutput);
xlabel('states in decimal'); ylabel('distribution'); title(str);
hold off


