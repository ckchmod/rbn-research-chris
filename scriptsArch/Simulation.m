%% Large scale simulations of RBNp's
clc; clear all;

topology = 'symmetric';
numCells = 10^2; % 20^2
numGenes = 11;   % 100
interaction = 1; % number of interacting cells (bandwidth)
perturb = .1; % noise
steps = 200;

mink = 1; maxk = 4;
minp = 0; maxp = 1;
increment = 0.020;
rowindex = length(minp:increment:maxp);
i = 1; j = 1;
simul_iter = 100; % Number of times you are repeating the same Boolean Network
simresults = zeros(maxk*(rowindex),5); % last two columns moment1 and moment2

for k = mink:maxk
    for p = minp:increment:maxp        
        %RBNpMoments = zeros(simul_iter, 2);
        RBNpTotal = zeros(simul_iter*nchoosek(numCells,2),1);
        for ss = 1:simul_iter
            % Simulating RBNp
            RBNp = boolCellGrid(topology, numCells, numGenes, k, p, ...
                interaction, [], [], [], perturb);
            RBNp.update_all(steps);
            RBNp_ssDist = ssDist(RBNp);
            KLDMatrix = KLDPairwise(RBNp_ssDist);
    
            %RBNpMoments(ss, 1) = mean(KLDMatrix);
            %RBNpMoments(ss, 2) = var(KLDMatrix);
            RBNpTotal(((ss-1)*nchoosek(numCells,2)+1):ss*nchoosek(numCells,2)) ...
                = KLDMatrix;
        end
        
        m1 = mean(RBNpTotal(:,1));
        m2 = var(RBNpTotal(:,1));
        simresults(j + (i-1)*rowindex,1) = k;
        simresults(j + (i-1)*rowindex,2) = p;
        simresults(j + (i-1)*rowindex,3) = 2*k*p*(1-p);
        simresults(j + (i-1)*rowindex,4) = m1;
        simresults(j + (i-1)*rowindex,5) = m2;        
        j = j+1;
    end 
    i = i +1;
    j = 1;
end
save('simresults.m');