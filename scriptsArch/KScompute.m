function [ KS ] = KScompute( allStates, numCells, numGenes, numSelect, sizeG)

    randSelect = randi([1 numCells*numGenes],1, numSelect);
    mat_steps = size(allStates,3) -1;
    length_group1 = 2:mat_steps/2;
    length_group2 = mat_steps/2+1:mat_steps;
    group1 = zeros(mat_steps/2, numSelect);
    group2 = zeros(mat_steps/2, numSelect);
    for i = 1:mat_steps/2
        for j = 1:numSelect
            row = floor((randSelect(j)+numGenes-1)/numGenes);
            col = mod(randSelect(j), numGenes); if (col==0) col=numGenes; end;
            group1(i, j) = allStates(row, col, i + 1); 
            group2(i, j) = allStates(row, col, i + 1 +mat_steps/2);
        end
    end
    
    gSpacing = sizeG:sizeG:size(group1, 1);
    group1_vec = zeros(1, size(gSpacing,2));
    group2_vec = group1_vec;

    for i = 1:size(gSpacing,2)
        group1_vec(1,i) = bi2de(group1(gspacing(i),:));
        group2_vec(1,i) = bi2de(group2(gspacing(i),:));
    end

    [h,p,k] = kstest2(group1_vec, group2_vec);
    KS = k;    
    
end

