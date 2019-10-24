%% Large scale simulations of RBNp's 
clc; clear all;

topology = 'symmetric';
numCells = 10^2; % 20^2
numGenes = 10;   % 100
interaction = 1; % number of interacting cells (bandwidth)
perturb = .1; % noise
steps = 1;
t_final = 0;

mink = 1; maxk = 2; kincrement = 1;
minp = .05; maxp = .95; pincrement = 0.05;
rowindex = length(minp:pincrement:maxp);
i = 1; j = 1;
simul_iter_low = 1;
simul_iter_high = 10; % Number of times repeating the same Boolean Network
threshold = 0.001;
noconv = false;

indexCount = 1;
dataRow = zeros(1, 8);
dataFile = [];
checkpointFilename = 'checkpoint.mat';
restart = false;

if (exist(checkpointFilename, 'file') == 2)
    s = load(checkpointFilename);
    dataFile = s.dataFile;
    indexCount = dataFile(end,1)+1;
    mink = dataFile(end, 2);
    minp = dataFile(end, 3);
    simul_iter_low = dataFile(end, 4);
    if (simul_iter_low == simul_iter_high)
        simul_iter_low = 1;
        minp = minp + pincrement;
        if (minp == maxp)
            minp = .5;
            mink = mink + kincrement;
        end
    else
        simul_iter_low = simul_iter_low+1;
    end
    fprintf('Restarting from iteration %d, k=%d, p=%d, sim=%d\n',...
        indexCount, mink, minp, simul_iter_low);
end
if ((mink == maxk) && (minp == maxp) && (simul_iter_low==simul_iter_high))
    disp('No more simulations to run...');
end

for k = mink:kincrement:maxk
    for p = minp:pincrement:maxp        
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
            RBNpdummy.update_all(1);
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
               converge = .5*(KLD(P,P_delta) + KLD(P_delta,P));
               if (converge < threshold)
                   % continue loop to check every KLD(P,P*) converges;
                   noconv = false;
               else
                   noconv = true;
                   break
               end
            end      
            %Test TBD
            cc = 1;
            while (noconv == true)
                lastRBNpStates = RBNp.allStates(:,:,end);                
                RBNpdummy = boolCellGrid(topology, numCells,numGenes, k, p, ...
                    interaction, lastRBNpStates, RBNp.initTtable, RBNp.initvarF, perturb);
                delta = size(RBNp.allStates,3)-1;
                RBNpdummy.update_all(1);
                lastRBNpdummyStates = RBNpdummy.allStates(:,:,2:end);                      
                RBNp_ssDist = RBNpstar_ssDist;
                RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
                RBNpstar_ssDist = ssDist(RBNp);   
                for i = 1:numCells
                    P = RBNp_ssDist(i, :);
                    P_delta = RBNpstar_ssDist(i, :);
                    converge = .5*(KLD(P,P_delta) + KLD(P_delta,P));
                    % test - TBD
                    if (i ==1)
                        convergePlot(cc) = converge;
                        cc = cc+1;
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
            KLDMatrix = KLDPairwise(RBNpstar_ssDist);
            t_final = size(RBNp.allStates,3)-1;
            dataRow(1, 1) = indexCount;
            dataRow(1, 2) = k;
            dataRow(1, 3) = p;
            dataRow(1, 4) = ss;
            dataRow(1, 5) = lyapunov(RBNp.initTtable,numCells,interaction);
            dataRow(1, 6) = mean(KLDMatrix);
            dataRow(1, 7) = var(KLDMatrix);
            dataRow(1, 8) = t_final;
            dataFile = [dataFile; dataRow];
            save('checkpoint.mat', 'dataFile'); 
            indexCount = indexCount +1;
            simul_iter_low=1;
            ss
            toc    
        end
        minp = .05;
        simul_iter_low=1;
        p
    end
    minp = .05;
    simul_iter_low=1;
    k
end
disp('Finished Simulation');