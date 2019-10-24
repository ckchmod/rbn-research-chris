%% Large scale simulations of RBNp's 
clc; clear all;

topology = 'symmetric';
numCells = 10^2; % 20^2
numGenes = 10;   % 100
interaction = 1; % number of interacting cells (bandwidth)
perturb = .1; % noise
steps = 20000;

mink = 1; maxk = 2;
minp = .05; maxp = .95;
increment = 0.05;
rowindex = length(minp:increment:maxp);
i = 1; j = 1;
simul_iter = 10; % Number of times you are repeating the same Boolean Network
simresults = zeros((maxk-mink+1)*(rowindex),6); % last two columns moment1 and moment2
threshold = 0.001;

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
            RBNpdummy.update_all(1);
            lastRBNpdummyStates = RBNpdummy.allStates(:,:,end);
            
            %Calculate Steady State Distribution Here;
            RBNp_ssDist = ssDist(RBNp); %nnz(RBNp_ssDist); 
            RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
            RBNpstar_ssDist = ssDist(RBNp);%nnz(RBNpstar_ssDist);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Check that the distribution is ergodic
            % Guaranteeing Ergodicitiy is hard-work
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
            % ad hoc with certain threshold
            % We may need mathematical proof for this?
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %count = 1;
            KLDMatrix = KLDPairwise(RBNp_ssDist);
            KLDMatrixStar = KLDPairwise(RBNpstar_ssDist);
            normKLD = (norm(KLDMatrix - KLDMatrixStar, 2));
            if (normKLD < threshold )
                % pass convergence
            else
                while (normKLD >= threshold)

                lastRBNpStates = RBNp.allStates(:,:,end);        
            
                RBNpdummy = boolCellGrid(topology, numCells,numGenes, k, p, ...
                    interaction, lastRBNpStates, RBNp.initTtable, RBNp.initvarF, perturb);
                RBNpdummy.update_all(1);
                lastRBNpdummyStates = RBNpdummy.allStates(:,:,end);
                            
                RBNp_ssDist = RBNpstar_ssDist;
                RBNp.allStates = cat(3, RBNp.allStates, lastRBNpdummyStates);
                RBNpstar_ssDist = ssDist(RBNp);

                KLDMatrix = KLDMatrixStar;
                KLDMatrixStar = KLDPairwise(RBNpstar_ssDist);
                normKLD = (norm(KLDMatrix - KLDMatrixStar, 2))
                
                %- Uncomment to count track & plot long term behavior -%
                %KLDifferenceVec(count) = normKLD;
                %count = count + 1; 
                    
                end
            end
            %plot([1:count], KLDifferenceVec);
            RBNpTotal(((ss-1)*nchoosek(numCells,2)+1):ss*nchoosek(numCells,2)) ...
                = KLDMatrix;
            LyapunovTotal(ss) = lyapunov(RBNp.initTtable,numCells,interaction);
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
%save('simresults.m');