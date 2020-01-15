% Large scale simulations of RBNp's 

clc; clear all;
load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat')

steps = 2^14;
interaction = 1;
topology = 'ising';
perturbation = 0.0;
numCells = 400;
J = -1;
init_temperature = 2.26918531421;
h_interac = 1;
 
% Need 
% low T and h =0, PRE + Prob canalyzing function
% low T and h != 0, ???
% high T and h != 0, ???    
%
% question, is the hamiltonian only computed inside the monte-carlo
% simmulation?

RBNpnew = boolCellGrid_v2( topology, numCells, RBNp.numGenes, RBNp.k, RBNp.p, ...
     interaction, [], RBNp.initTtable, RBNp.initvarF, perturbation, J, init_temperature, h_interac);
RBNpnew.update_all(steps);
RBNp_ssDist = ssDist(RBNpnew);

% figure(3);
% for i = 1:numCells
%     subplot(sqrt(numCells), sqrt(numCells), i)
%     hold on
%     set(gca, 'YScale', 'log')
%     s1 = bar(1:length(RBNp_ssDist(i,:)), RBNp_ssDist(i,:), 5);
%     str=sprintf('Cell Number: %d', i);
%     title(str);axis([0 1024 0 0.04]);
%     hold off
%     alpha(s1,.5);
% end
 
% figure(1);
% for i = 1:viz_cells
%     if (i <= 3)
%         j = i;     
%         if (mod(j,2) == 1)
%             col = 'g';
%         else
%             col =  'b';
%         end
%      elseif (i <=6)
%         j = i+7;
%         if (mod(j,2) == 1)
%             col = 'b';
%         else
%             col =  'g';
%         end 
%      else
%         j = i+14;
%         if (mod(j,2) == 1)
%             col = 'g';
%         else
%             col = 'b';
%         end
%     end  
%      %j=i;
%      subplot(3,3, i)
%      hold on
%      set(gca, 'YScale', 'log')
%      s1 = stem(1:length(RBNp_ssDist(j,:)), RBNp_ssDist(j,:),col);
% 
%      %s1 = bar(1:length(RBNp_ssDist(j,:)), RBNp_ssDist(j,:), 5)
%      str=sprintf('Cell Number: %d', j);
%      title(str);axis([0 1024 0 0.04]);
%      hold off
%      alpha(s1,.5);
% end
