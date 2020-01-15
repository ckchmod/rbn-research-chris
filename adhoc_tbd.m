clear all; clc;
load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpChaotic.mat')
chaotic_states = RBNp.allStates;
load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpOrdered.mat')
ordered_states = RBNp.allStates;
load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat')
critical_states = RBNp.allStates;
% 
% S_(n,m,k,t), where S is the value ( 0 or 1) of the (n,m,k) node, n=1,…,N, m=1,…,M, at genes k=1,…,K, and t=1,…, Ntot (Ntot is the total time steps)
% 
% S_(1,1,1,1) S_(1,1,1,2) ……………...…. S_(1,1,1,Ntot) 
% S_(1,1,2,1) S_(1,1,2,2)                                   .   
% .                                                                       .
% .                                                                       .
% S_(1,1,K,1) S_(1,1,K,2)                                   .
% S_(1,2,1,1) S_(1,2,1,2)                                    .
% S_(1,2,2,1) S_(1,2,2,2)                                    .
% .                                                                        .
% .
% S(N,M,K,1) S_(N,M,K,2) ……………...…. S_(N,M,K,Ntot)

ordered_write = zeros(1000, size(ordered_states,3));
for i = 1:size(ordered_states,3) 
    A = ordered_states(:,:,i)';
    A = A(:);
    ordered_write(:,i) = A;
end

critical_write = zeros(1000, size(critical_states,3));
for i = 1:size(critical_states,3) 
    A = critical_states(:,:,i)';
    A = A(:);
    critical_write(:,i) = A;
end

chaotic_write = zeros(1000, size(chaotic_states,3));
for i = 1:size(chaotic_states,3) 
    A = chaotic_states(:,:,i)';
    A = A(:);
    chaotic_write(:,i) = A;
end
