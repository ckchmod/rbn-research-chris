%% Large scale simulations of RBNp's 
clc; clear all;

%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat');
%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpOrdered.mat');
%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/RBNpChaotic.mat');
%load('/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier01_35.mat');


myDir = '/home/chrisk/Desktop/MATLAB-scripts/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/'; 
myFiles = dir(fullfile(myDir,'out*'));
file_size = size(myFiles,1);
s.name = [];
s.lyap = [];
s(1) = [];
for k = 1:length(myFiles)
  st =struct;
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  fprintf(1, 'Now reading %s\n', baseFileName);
  load(fullFileName);
  lyap = RBNp.crit_val;
  st.name = fullFileName;
  st.lyap = lyap;
  s(end+1) = st;
end


T = struct2table(s); % convert the struct array to a table
sortedT = sortrows(T, 'lyap');

hold on
for k = 1:length(myFiles)
  %baseFileName = myFiles(k).name;
  FileName = sortedT.name{k};
  fprintf(1, 'Now reading %s\n', FileName);
  load(FileName);

    % inNodes = RBNp.allStates(:,1,:); %First column is all in Nodes.
    % outNodes = RBNp.allStates(:,end,:); %last column is all end Nodes.
    initTtable = RBNp.initTtable;
    initvarF = RBNp.initvarF;
    initnv = zeros(1, size(initvarF,2));
    numGenes = RBNp.numGenes;
    for i = 1:size(initvarF,2)
        count = 0;
        for j = 1:size(initvarF, 1)
            if(initvarF(j,i) >= 0)
                count = count + 1;
            end
        end
        initnv(i) = count;
    end

    stateSpace = [];
    for i = 0:(2^numGenes -1) % This should be numGenes^10 -1, but is fine as of now
        num = [i, bi2de(bnNextState(de2bi(i, numGenes), initTtable, initvarF, initnv))];
        stateSpace = [stateSpace; num];
    end
    s = stateSpace(:,1)' + 1;
    t = stateSpace(:,2)' + 1;
    G = digraph(s, t);
    figure(1)
    subplot(13,13,k)
    plot(G,'Layout','force');
    title(RBNp.crit_val);

end
hold off
