numCells = 6^2;
interaction = 2;
numGenes = 10;
k = 5;
p = .25;
RBNp = boolCellGrid('symmetric', numCells, numGenes, k, p, interaction, [], [], [], .2);
ttable = RBNp.initTtable;
row = size(ttable,1);
col = size(ttable,2);
LTcolvector = 2^4; %linear threshold column vector (N,S,E,W)
LTvector = [0; 0; 0; 1; 0; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1];
 
if (row < LTcolvector) 
    bigmatrix = zeros(LTcolvector, numCells*col)-1; 
else
    bigmatrix = zeros(row, numCells,col);
end

for k = 1:row;
    for i = 1:numCells
        for j = 1:col
            bigmatrix(k, j+(i-1)*col) = ttable(k,j);
        end
    end
end
for inter = 1: interaction
    for repcol = 1:numCells
        for repi = 1:size(LTvector,1)
            bigmatrix(repi, (col*(repcol-1)+1)+(inter-1)) = LTvector(repi);
        end
    end   
end
tempA = bnActivity(bigmatrix);
temp = sum(tempA,2);
lambda = log(mean(temp))

