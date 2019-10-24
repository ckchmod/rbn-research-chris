function lambda = lyapunov( F, Cells, Interaction)
    %Takes in the truth table, applys LinearThreshold(2), on colum

    ttable = F;
    numberC = Cells;
    interaction = Interaction;
    row = size(ttable,1);
    col = size(ttable,2);
    LTcolvector = 2^4; %linear threshold column vector (N,S,E,W)
    LTvector = [0; 0; 0; 1; 0; 1; 1; 1; 0; 1; 1; 1; 1; 1; 1; 1];

    if (row < LTcolvector) 
        bigmatrix = zeros(LTcolvector, numberC*col)-1; 
    else
        bigmatrix = zeros(row, numberC,col);
    end

    for k = 1:row;
        for i = 1:numberC
            for j = 1:col
                bigmatrix(k, j+(i-1)*col) = ttable(k,j);
            end
        end
    end
    for interC = 1: interaction
        for repcol = 1:numberC
            for repi = 1:size(LTvector,1)
                bigmatrix(repi, (col*(repcol-1)+1)+(interC-1)) = LTvector(repi);
            end
        end   
    end
    
    tempA = bnActivity(bigmatrix);
    temp = sum(tempA,2);
    lambda = log(mean(temp));

end

