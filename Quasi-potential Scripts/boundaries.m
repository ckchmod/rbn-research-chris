function [N,S,W,E] = boundaries(i,j, matrix)

    N = i-1; S = i+1; E = j+1; W = j-1;
    if N == 0 
        N = size(matrix,1);
    end 
    if S == size(matrix,1)+1
        S = 1;
    end
    if W == 0 
        W = size(matrix,2);
    end 
    if E == size(matrix,2)+1
        E = 1;
    end 

end

