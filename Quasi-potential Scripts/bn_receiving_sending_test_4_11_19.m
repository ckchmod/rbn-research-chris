%% Large scale simulations of RBNp's 
clc; clear all;

%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat');
%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpOrdered.mat');
load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpChaotic.mat');


inNodes = RBNp.allStates(:,1,:); %First column is all in Nodes.
outNodes = RBNp.allStates(:,end,:); %last column is all end Nodes.

Cells = 3;
%T = 4000;
T = 3000;
figure(1)

hold on
for cells = 1:Cells
    x = reshape(inNodes(cells,1,1:T), [1,T]);
    y = reshape(outNodes(cells,1,1:T), [1,T]);
    subplot(2,1,1);
    plot(1:T, x); ylim([-1 2]); xlabel('time'); ylabel('states'); title('time evolution of states of receiving nodes')
    subplot(2,1,2);
    plot(1:T, y); ylim([-1 2]); xlabel('time'); ylabel('states'); title('time evolution of states of ouptut nodes')
end
hold off

figure(2) 
title('In Nodes')
for i = 1:100
     subplot(10,10, i)
     hold on
     s1 = reshape(inNodes(i,1,1:T), [1,T]);
     histogram(s1)
end
hold off

figure(3)
title('Out Nodes')
for i = 1:100
     subplot(10,10, i)
     hold on
     s1 = reshape(outNodes(i,1,1:T), [1,T]);
     histogram(s1)
    
end
hold off