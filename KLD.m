function [ rhs ] = KLD(P, Q)

    % The Kullback–Leibler divergence is defined only if Q(i)=0 implies P(i)=0, for all i (absolute continuity). 
    % Common practice (throwout 0's)
    % http://mathoverflow.net/questions/72668/how-to-compute-kl-divergence-when-pmf-contains-0s

    temp =  P.*log(P./Q);
    temp(isnan(temp))=0; % resolving the case when P(i)==0
    temp(isinf(temp))=0; % resolving the case when Q(i)==0 ??
    rhs = sum(temp,2);

end
