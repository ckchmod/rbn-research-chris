function [ KLDMatrix ] = KLDPairwise( ssD )
    % Pass in Steady State Distribution, and return pairwise KLDMatrix

    pairwise = nchoosek(1:length(ssD(:,1)), 2);
    KLDMatrix = zeros(length(pairwise(:,1)), 1);

    for i = 1:length(KLDMatrix(:,1))
        P = ssD(pairwise(i,1),:);
        Q = ssD(pairwise(i,2),:);
        % KLD Symmetric
        KLDMatrix(i) = .5*(KLD(P,Q) + KLD(Q,P));
    end
    
end

