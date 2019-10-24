orderedLyap = zeros(1,binsize);

Ni = ones(1,binsize);
for l = 1:length(emptyList(:,2))
    for j = 1:length(N)-1
        %simresults(l,5) >= N(n) && simresults(l,5) < N(n+1))
        if (ordered(l,5) >= N(j) && ordered(l,5) < N(j+1))
            orderedLyap(Ni(j),j)= ordered(l,6);
            Ni(j) = Ni(j)+1;
            break;
        else
            %continue;
        end
    end
end
%orderedLyap(orderedLyap==0) = [];


T = bplot(orderedLyap)