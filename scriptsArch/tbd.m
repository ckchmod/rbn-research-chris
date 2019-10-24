%clc; clear all;

topology = 'symmetric';
numCells = 10^2; % 20^2
numGenes = 10;   % 100
interaction = 1; % number of interacting cells (bandwidth)
perturb = .3; % noise
steps = 1;
t_final = 0;
k = 2;
p = 0.1;

numSelect = 10;
sizeG = 2;
threshold = 0.02;
% dataFrame = zeros(17100,4);
% count =1;
%  for k = 3:3
%      for p = .05:.05:.95
%          for iter = 1:100
%              RBNp = boolCellGrid(topology, numCells, numGenes, k, p, ...
%              interaction, [], [], [], perturb);
%              dataFrame(count,1) = k;
%              dataFrame(count,2) = p;
%              dataFrame(count,3) = iter;
%              lam = lyapunov(RBNp.initTtable,numCells,interaction)
%              dataFrame(count,4) = lam;
%              count = count +1
%          end
%      end
%  end

k = 1;
p = .05;
simul_iter_low = 1 ;
simul_iter_high = 1;
for ss = simul_iter_low:simul_iter_high
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Simulating RBNp
            % Create initial RBNp and RBNp*
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    tic
    RBNp = boolCellGrid(topology, numCells, numGenes, k, p, ...
        interaction, [], [], [], perturb);
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check for convergence, 
    % ad hoc with certain threshold\
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
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
    cc=1;
    while (noconv == true)
        lastRBNpStates = RBNp.allStates(:,:,end);                
        RBNpdummy = boolCellGrid(topology, numCells,numGenes, k, p, ...
            interaction, lastRBNpStates, RBNp.initTtable, RBNp.initvarF, perturb);
        delta = size(RBNp.allStates,3)-1
        RBNpdummy.update_all(delta);
        lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end);                      
        RBNp_ssDist = RBNpstar_ssDist;
        RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
        RBNpstar_ssDist = ssDist(RBNp);   
        for i = 1:numCells
            P = RBNp_ssDist(i, :);
            P_delta = RBNpstar_ssDist(i, :);
            converge = .5*(KLD(P,P_delta) + KLD(P_delta,P));
            if (i ==1)
                convergePlot(cc) = converge
                cc = cc+1
            end
            if (converge < threshold)
                % continue loop to check every KLD(P,P*) converges;
                noconv = false;
            else
                noconv = true;
                break
            end
        end  
    end
%     KLDMatrix = KLDPairwise(RBNpstar_ssDist);
%     t_final = size(RBNp.allStates,3)-1;
%     dataRow(1, 1) = indexCount;
%     dataRow(1, 2) = k;
%     dataRow(1, 3) = p;
%     dataRow(1, 4) = ss;
%     dataRow(1, 5) = lyapunov(RBNp.initTtable,numCells,interaction);
%     dataRow(1, 6) = mean(KLDMatrix);
%     dataRow(1, 7) = var(KLDMatrix);
%     dataRow(1, 8) = t_final;
%     dataFile = [dataFile; dataRow];
%     save('checkpoint.mat', 'dataFile'); 
%     indexCount = indexCount +1;
%     simul_iter_low=1;
%     ss
    disp('el fin');
    toc    
end
         
         