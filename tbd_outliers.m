clear all; clc;
% o1 = load('outlier_140.mat'); o1 = o1.RBNp;
% o2 = load('outlier_201.mat'); o2 = o2.RBNp;
% o3 = load('outlier_212.mat'); o3 = o3.RBNp;
% o4 = load('outlier_266.mat'); o4 = o4.RBNp;
 o5 = load('/home/charlestreykang/Desktop/MATLAB/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat'); o5 = o5.RBNp;
% 
% o1_ssd = ssDist(o1);
% o2_ssd = ssDist(o2);
% o3_ssd = ssDist(o3);
% o4_ssd = ssDist(o4);
 o5_ssd = ssDist(o5);
% 
% o1_Tt = o1.initTtable;
% o2_Tt = o2.initTtable;
% o3_Tt = o3.initTtable;
% o4_Tt = o4.initTtable;

RBNp = o5; %test
k=2; p=.55; topology = 'symmetric'; numCells = 9^2; 
threshold = .02;
numGenes = 10; interaction = 1;  perturb = .00; % noise
lastRBNpStates = RBNp.allStates(:,:,end);                
RBNpdummy = boolCellGrid(topology, numCells,numGenes, k, p, ...
    interaction, [], RBNp.initTtable, RBNp.initvarF, perturb);
delta = size(RBNp.allStates,3)-1;
RBNpdummy.update_all(65536);
lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end); 
RBNpstar_ssDist = lastRBNpdummyStates;
RBNpstar_ssDist=ssDist(RBNpdummy);
% RBNp_ssDist = o4_ssd;
% RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
% RBNpstar_ssDist = ssDist(RBNp);   
% for i = 1:numCells
%     P = RBNp_ssDist(i, :);
%     P_delta = RBNpstar_ssDist(i, :);
%     converge = .5*(KLD(P,P_delta) + KLD(P_delta,P));
%     if (converge < threshold)
%         % continue loop to check every KLD(P,P*) converges;
%         noconv = false;
%     else
%         noconv = true;
%         break
%     end
% end  
KLDMatrix = KLDPairwise(RBNpstar_ssDist);
newKLD = mean(KLDMatrix)
newKLDv = var(KLDMatrix)