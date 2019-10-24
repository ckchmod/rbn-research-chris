clear all; clc;
% % o1 = load('outlier_140.mat'); o1 = o1.RBNp;
% % o2 = load('outlier_201.mat'); o2 = o2.RBNp;
% % o3 = load('outlier_212.mat'); o3 = o3.RBNp;
% % o4 = load('outlier_266.mat'); o4 = o4.RBNp;
% o5 = load('/home/charlestreykang/Desktop/MATLAB/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat');% o5 = o5.RBNp;
% % 
% % o1_ssd = ssDist(o1);
% % o2_ssd = ssDist(o2);
% % o3_ssd = ssDist(o3);
% % o4_ssd = ssDist(o4);
% % o5_ssd = ssDist(o5);
% % 
% % o1_Tt = o1.initTtable;
% % o2_Tt = o2.initTtable;
% % o3_Tt = o3.initTtable;
% % o4_Tt = o4.initTtable;
% 
% o5 = load('objecTest.mat'); o5 = o5.RBNp;
% 
% RBNp = o5; %test
load('objTest_initStates.mat');load('objTest_initT.mat');load('objTest_initvarF.mat');
k=2; p=.55; topology = 'symmetric'; numCells = 10^2; 
threshold = .02;
numGenes = 10; interaction = 1;  perturb = .01; % 0.01 is default noise
%lastRBNpStates = RBNp.allStates(:,:,end);                
RBNpdummy = boolCellGridPerturb(topology, numCells,numGenes, k, p, ...
    interaction, initStates, initT, initvarF, perturb);
%delta = size(RBNp.allStates,3)-1;
RBNpdummy.update_all(65536); %65536
%lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end); 
%RBNpstar_ssDist = lastRBNpdummyStates;
RBNpstar_ssDist = ssDist(RBNpdummy);
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
% KLDMatrix = KLDPairwise(RBNpstar_ssDist);
% newKLD = mean(KLDMatrix)
% newKLDv = var(KLDMatrix)

RBNp_ssDist = RBNpstar_ssDist;
figure(1);
for i = 1:16
%      if (i <= 3)
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
    if (i <= 4)
        j=i;
        if (mod(i,2) == 1)
            col = 'g';
        else
            col = 'b';
        end
    elseif (i <= 8)
        j=i+6;
        if (mod(i,2) == 1)
            col = 'b';
        else
            col = 'g';
        end
    elseif (i <= 12)
        j=i+12;
        if (mod(i,2) == 1)
            col = 'g';
        else
            col = 'b';
        end
    else
        j = i+18;
        if (mod(i,2) == 1)
           col = 'b';
        else
            col = 'g';
        end
    end
     %j=i;
     subplot(4,4, i)
     hold on
     s1 = scatter(1:length(RBNp_ssDist(j,:)), RBNp_ssDist(j,:), 'filled',col)
     str=sprintf('Cell Number: %d', j);
     title(str);axis([0 1024 0 0.04]);
     hold off
     alpha(s1,.5);
end