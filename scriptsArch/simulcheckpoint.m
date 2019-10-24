%% Large scale simulations of RBNp's 
clc; clear all;

%%%%%%%%%%%%% Checkpoint Setup %%%%%%%%%%%%%%%%%%%%%%%%%
% numIter = 10;
% startIter = 1;
% checkpointFilename = 'checkpoint.mat';
% 
%  if exist(checkpointFilename, 'file')
%     s = load(checkpointFilename);
%     startIter = s.i;
%     fprintf('Restarting from iteration %d\n', startIter);
%  end
% 
%   for i = startIter:numIter
%     fprintf('Starting iteration %d\n', i);
% %     expensiveComputation();
% %     save(checkpointFilename, 'i');
% %   end
%   % We succefully finished. Let's delete our checkpoint file
%   delete(checkpointFilename);
% 
%  % function expensiveComputation()
%     % Pretend to do lots of work!
%     pause(1);
%  % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
topology = 'symmetric';
numCells = 10^2; % 20^2
numGenes = 10;   % 100
interaction = 1; % number of interacting cells (bandwidth)
perturb = .1; % noise
steps = 100;
t_final = 0;

mink = 1; maxk = 2;
minp = .05; maxp = .95;
increment = 0.05;
rowindex = length(minp:increment:maxp);
i = 1; j = 1;
simul_iter = 10; % Number of times you are repeating the same Boolean Network
simresults = zeros((maxk-mink+1)*(rowindex),6); % last two columns moment1 and moment2
threshold = 0.001;
noconv = false;

distributionSize = numCells*(2^numGenes);

for k = mink:maxk
    for p = minp:increment:maxp        
        RBNpTotal = zeros(simul_iter*nchoosek(numCells,2),1);
        LyapunovTotal = zeros(simul_iter);
        for ss = 1:simul_iter
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
            RBNpdummy.update_all(size(RBNp.allStates,3)-1);
            lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end);
            
            %Calculate Steady State Distribution Here;
            RBNp_ssDist = ssDist(RBNp); %nnz(RBNp_ssDist); 
            RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
            RBNpstar_ssDist = ssDist(RBNp);%nnz(RBNpstar_ssDist);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Uncomment for Absolute Ergodicitiy Guarantee for KLD
            % Else, uses approx KLD
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             while ((distributionSize-nnz(RBNp_ssDist)==0)==0 && ...
%                 (distributionSize-nnz(RBNpstar_ssDist)==0)==0)
% 
%                 lastRBNpStates = RBNp.allStates(:,:,end);        
%             
%                 RBNpdummy = boolCellGrid(topology, numCells,numGenes, k, p, ...
%                     interaction, lastRBNpStates, RBNp.initTtable, RBNp.initvarF, perturb);
%                 RBNpdummy.update_all(1);
%                 lastRBNpdummyStates = RBNpdummy.allStates(:,:,end);
%                             
%                 RBNp_ssDist = ssDist(RBNp);
%                 RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
%                 RBNpstar_ssDist = ssDist(RBNp);
%                 % size(RBNp.allStates,3)
%                 % nnz(RBNpstar_ssDist)                
%             end
%             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check for convergence, 
            % ad hoc with certain threshold\
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            for i = 1:numCells
               P = RBNp_ssDist(i, :);
               P_delta = RBNpstar_ssDist(i, :);
               % KLD(P, PD)
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
                RBNpdummy.update_all(size(RBNp.allStates,3)-1);
                lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end);
                            
                RBNp_ssDist = RBNpstar_ssDist;
                RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
                RBNpstar_ssDist = ssDist(RBNp);
    
                size(RBNp.allStates,3)-1
                
                for i = 1:numCells
                    P = RBNp_ssDist(i, :);
                    P_delta = RBNpstar_ssDist(i, :);
                    % KLD(P, PD)
                    converge = .5*(KLD(P,P_delta) + KLD(P_delta,P))
                    if (converge < threshold)
                        % continue loop to check every KLD(P,P*) converges;
                        noconv = false;
                    else
                        noconv = true;
                        break
                    end
                end  
                   
            end
            KLDMatrix = KLDPairwise(RBNpstar_ssDist);
            RBNpTotal(((ss-1)*nchoosek(numCells,2)+1):ss*nchoosek(numCells,2)) ...
                = KLDMatrix;
            LyapunovTotal(ss) = lyapunov(RBNp.initTtable,numCells,interaction);
            t_final = size(RBNp.allStates,3)-1;
            ss
            toc
        end
        
        m1 = mean(RBNpTotal(:,1));  
        m2 = var(RBNpTotal(:,1));
        lyapmean = mean(LyapunovTotal(:,1));
        simresults(j + (i-1)*rowindex,1) = k;
        simresults(j + (i-1)*rowindex,2) = p;
        simresults(j + (i-1)*rowindex,3) = lyapmean;
        simresults(j + (i-1)*rowindex,4) = m1;
        simresults(j + (i-1)*rowindex,5) = m2;
        simresults(j + (i-1)*rowindex,6) = 2*k*p*(1-p);
        j = j+1;
        p
    end 
    i = i +1;
    j = 1;
    k
end
save('simresults.m');