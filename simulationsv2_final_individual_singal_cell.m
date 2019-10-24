%% Large scale simulations of RBNp's 
clc; clear all;
load('/home/charlestreykang/Desktop/MATLAB/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat')

topology = 'single'; 
numCells = 1; 
numGenes = 10;   
interaction = 1;
perturb = .05; %.01; 
steps = 2^14;

%o5 = load('outlier_481.mat'); o5 = o5.RBNp; RBNp = o5; %test
%outlier example symmetric diversifcation
 k=2; p=.55; 
initState = RBNp.initState;
initTable = RBNp.initTtable; 
initVar = RBNp.initvarF;
% 
RBNpnew = boolCellGrid(topology, numCells,numGenes, k, p, ...
     interaction, initState, initTable, initVar, perturb);
RBNpnew.update_all(steps);
RBNp_ssDist = ssDist(RBNpnew);

% %line example symmetric diversifcation test
% k=2; p=.45; 
% 
% RBNpnew = boolCellGrid(topology, numCells,numGenes, k, p, ...
%      interaction, [], [], [], perturb);
% RBNpnew.update_all(steps);
% RBNp_ssDist = ssDist(RBNpnew);
% 

% hold on
% j =1;  col =  'r';
% set(gca, 'YScale', 'log')
% s1 = stem(1:length(RBNp_ssDist(j,:)), RBNp_ssDist(j,:),col)
% str=sprintf('Cell Number: %d', j);
% title(str);axis([0 1024 0 0.04]);
% hold off
% alpha(s1,.5);

RBNp_ssDist = ssDist(RBNpnew);
figure(2);
% for i = 1:9
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
%      j=i;
%      subplot(3,3, i)
%      hold on
%      set(gca, 'YScale', 'log')
%      s1 = stem(1:length(RBNp_ssDist(j,:)), RBNp_ssDist(j,:),col)
%      str=sprintf('Cell Number: %d', j);
%      title(str);axis([0 1024 0 0.04]);
%      hold off
%      alpha(s1,.5);
% end

for i = 1:1
    col = 'r'
     j=i;
     subplot(1,1, i)
     hold on
     set(gca, 'YScale', 'log')
     s1 = stem(1:length(RBNp_ssDist(j,:)), RBNp_ssDist(j,:),col)
     str=sprintf('Cell Number: %d', j);
     title(str);axis([0 1024 0 0.04]);
     hold off
     alpha(s1,.5);
end


threshold = 0.02;
noconv = false;
dataRow = zeros(1, 6);
convergeRow = [];
dataFile = [];
convergeFile = [];
checkpointFilename = 'perturb.mat';
restart = false;

if (exist(checkpointFilename, 'file') == 2)
    s = load(checkpointFilename);    
    dataFile = s.dataFile;
    convergeFile = s.convergeFile;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulating RBNp
% Create initial RBNp and RBNp*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
tic
RBNp = boolCellGrid(topology, numCells,numGenes, k, p, ...
    interaction, [], [], [], perturb);
lyapunov(RBNp.initTtable,numCells,interaction)
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

for i = 1:numCells
   P = RBNp_ssDist(i, :);
   P_delta = RBNpstar_ssDist(i, :);
   converge = .5*(KLD(P,P_delta) + KLD(P_delta,P));
   if (converge < threshold)
       % continue loop to check every KLD(P,P*) converges;
       noconv = false;
   else
       noconv = true;
       break
   end
end      

while (noconv == true)
    lastRBNpStates = RBNp.allStates(:,:,end);                
    RBNpdummy = boolCellGrid(topology, numCells,numGenes, k, p, ...
        interaction, lastRBNpStates, RBNp.initTtable, RBNp.initvarF, perturb);
    delta = size(RBNp.allStates,3)-1;
    RBNpdummy.update_all(delta);
    lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end);                      
    RBNp_ssDist = RBNpstar_ssDist;
    RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
    RBNpstar_ssDist = ssDist(RBNp);   
    for i = 1:numCells
        P = RBNp_ssDist(i, :);
        P_delta = RBNpstar_ssDist(i, :);
        converge = .5*(KLD(P,P_delta) + KLD(P_delta,P));
        if (converge < threshold)
            % continue loop to check every KLD(P,P*) converges;
            noconv = false;
        else
            noconv = true;
            break
        end
    end  
    convergeRow = [convergeRow, converge];
    converge
end
KLDMatrix = KLDPairwise(RBNpstar_ssDist);
t_final = size(RBNp.allStates,3)-1;
dataRow(1, 1) = k;
dataRow(1, 2) = p;
dataRow(1, 3) = lyapunov(RBNp.initTtable,numCells,interaction);
dataRow(1, 4) = mean(KLDMatrix);
dataRow(1, 5) = var(KLDMatrix);
dataRow(1, 6) = t_final;
toc    
figure
hold on
str = sprintf('g = .95, p = 0.0, t= %f, MeanKLD: %f', t_final, mean(KLDMatrix));
title(str)
for i = 1:9
    subplot(3,3, i)
    hold on
    s1 = scatter(1:length(RBNp_ssDist(i,:)), RBNp_ssDist(i,:),'filled')
    %s2 = scatter(1:length(RBNp_ssDist(i,:)), RBNpstar_ssDist(i,:), 'filled')
    hold off
    alpha(s1,.5);%alpha(s2,.5);
end
hold off


    
    % dataFile = [dataFile; dataRow];
% convergeFile = [convergeFile; convergeRow];
% save('perturb.mat', 'dataFile', 'convergeFile');