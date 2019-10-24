%% All Visualizations Used for Paper
clear all; clc;

%% Scatter Plot Lyapunov Exponent vs. Average KLD
% %load('batchsum.mat');
a = load('gbatchsum10by10t1.mat');
b = load('gbatchsum10by10t2.mat');
c = load('gbatchsum10by10t3.mat');
binsize = 10;
simul_iter=25;

figure
subplot(3,2,1)
scatter(a.dataFile(:,5),a.dataFile(:,6), 'filled','MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Average of D_{KLS}');
%title(['Sensitivity vs mean(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(a)');
axis([-3 3 0 4])
%title(['Sensitivity vs mean(Diversity) 4x4 l_{tv}=3']);

subplot(3,2,2)
scatter(a.dataFile(:,5),a.dataFile(:,7), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Variance of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(b)');
axis([-3 3 0 4])
%title(['Sensitivity vs var(Diversity) 4x4 l_{tv}=3']);
% 
subplot(3,2,3)
scatter(b.dataFile(:,5),b.dataFile(:,6), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Average of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(c)');
axis([-3 3 0 4])
subplot(3,2,4)
scatter(b.dataFile(:,5),b.dataFile(:,7), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Variance of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(d)');
axis([-3 3 0 4])
subplot(3,2,5)
scatter(c.dataFile(:,5),c.dataFile(:,6), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Average of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(e)');
axis([-3 3 0 4])
subplot(3,2,6)
scatter(c.dataFile(:,5),c.dataFile(:,7), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Variance of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(f)');
axis([-3 3 0 4])

% colored 

asked_p = 2;
figure
subplot(3,2,1)
scatter(a.dataFile(:,5),a.dataFile(:,6), 25 ,a.dataFile(:,asked_p), 'filled', 'MarkerFaceAlpha',.5 );xlabel('Lyapunov Exponent'); ylabel('Average of D_{KLS}');
%title(['Sensitivity vs mean(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(a)');
axis([-3 3 0 4])
%title(['Sensitivity vs mean(Diversity) 4x4 l_{tv}=3']);

subplot(3,2,2)
scatter(a.dataFile(:,5),a.dataFile(:,7), 25 ,a.dataFile(:,asked_p), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Variance of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(b)');
axis([-3 3 0 4])
%title(['Sensitivity vs var(Diversity) 4x4 l_{tv}=3']);
% 
subplot(3,2,3)
scatter(b.dataFile(:,5),b.dataFile(:,6), 25 ,b.dataFile(:,asked_p), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Average of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(c)');
axis([-3 3 0 4])
subplot(3,2,4)
scatter(b.dataFile(:,5),b.dataFile(:,7), 25 ,b.dataFile(:,asked_p), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Variance of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(d)');
axis([-3 3 0 4])
%
subplot(3,2,5)
scatter(c.dataFile(:,5),c.dataFile(:,6), 25 , c.dataFile(:,asked_p), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Average of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(e)');
axis([-3 3 0 4])
%
subplot(3,2,6)
scatter(c.dataFile(:,5),c.dataFile(:,7), 25 ,c.dataFile(:,asked_p), 'filled', 'MarkerFaceAlpha',.5);xlabel('Lyapunov Exponent'); ylabel('Variance of D_{KLS}');
%title(['Sensitivity vs var(Diversity) All Experiments (each RBNp repeated ' num2str(simul_iter) ...
%    ' times)']);
title('(f)');
axis([-3 3 0 4])
hp4 = get(subplot(3,2,6),'Position')
h = colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.01  hp4(2)+hp4(3)*2.1])
xlabel(h, 'k')

%% Distribution of Lyapunov exponents
%lyapunov
dataFile = a.dataFile;
[Y, N] = histcounts(dataFile(:,5), binsize, 'Normalization', 'Probability');
iidistrib = 20./Y';
sumdistrib = iidistrib./sum(iidistrib);
ordered = sortrows(dataFile,5);
emptyList = zeros(1,3);
count = 1;
for l = 1:length(dataFile(:,5))
    %n=1;
    for j = 1:length(N)-1
        %simresults(l,5) >= N(n) && simresults(l,5) < N(n+1))
        if (ordered(l,5) >= N(j) && ordered(l,5) < N(j+1))
            break;
        else
            %continue;
        end
    end
    j;
    rando = rand();
    if rando <= sumdistrib(j);
       emptyList(count, 1) = ordered(l, 5);
       emptyList(count, 2) = ordered(l, 6);
       emptyList(count, 3) = ordered(l, 7);
       emptyList(count, 4) = ordered(l, 8);
       count = count+1;
    end
        
end

numberOfBins = 20;
[counts, binValues] = hist(ordered(:,5), numberOfBins);
%bar(binValues, counts, 'barwidth', 1);
normalizedCounts = counts / sum(counts);
bar(binValues, normalizedCounts, 'barwidth', 1);
xlabel('Lyapunov Exponent');
ylabel('% count');
title('Distribution of LEs: 10x10 l_{tv}=1');
axis([-2.5 2.5 0 .13]);
% subplot(2,2,3);
% scatter(emptyList(:,1),emptyList(:,2));xlabel('Lyapunov Exponent'); ylabel('Average of KLD Distances');
% %title(['Sensitivity vs mean(Diversity) IID Selection (each RBNp repeated ' num2str(simul_iter) ...
% %    ' times)']);ylim=[0 0.2];
% title(['Sensitivity vs mean(Diversity) Uniformly Distributed Experiments']);
% 
% subplot(2,2,4);
% scatter(emptyList(:,1),emptyList(:,3));xlabel('Lyapunov Exponent'); ylabel('Variance of KLD Distances');
% %title(['Sensitivity vs var(Diversity) IID Selection (each RBNp repeated ' num2str(simul_iter) ...
% %    ' times)']);ylim=[0 0.2];
% title(['Sensitivity vs var(Diversity) Uniformly Distributed Experiments']);

% %% Convergence Time
% figure
% subplot(2,1,1);
% scatter(dataFile(:,5), dataFile(:,8));xlabel('Lyapunov Exponent'); ylabel('convergence time 2^t');
% title(['Lyapunov Exponent vs convergence time']);
% 
% subplot(2,1,2);
% scatter(emptyList(:,1), emptyList(:,4));xlabel('Lyapunov Exponent'); ylabel('convergence time 2^t');
% title(['Lyapunov Exponent vs convergence time (uniform distribution)']);

%% LOWESS Plots

% datain1 = [dataFile(:,5)';dataFile(:,6)']';
% datain2 = [dataFile(:,5)';dataFile(:,7)']';
% %datain3 = [emptyList(:,1)';emptyList(:,2)']';
% %datain4 = [emptyList(:,1)';emptyList(:,3)']';
% f = 0.2;
% wantplot = 1;
% [dataout1 lowerLimit upperLimit xy] = lowess(datain1,f,0);
% [dataout2 lowerLimit upperLimit xy] = lowess(datain2,f,0);
% %[dataout3 lowerLimit upperLimit xy] = lowess(datain3,f,0);
% %[dataout4 lowerLimit upperLimit xy] = lowess(datain4,f,0);
% subplot(2,2,3)
% plot(dataout1(:,1), dataout1(:,2), 'oblue',dataout1(:,1),dataout1(:,3),...
%        '-red');xlabel('Lyapunov Exponent'); ylabel('Mean KLD Distances');
% title('(c)');
% axis([-3 3 0 .05]);
% %title(['Sensitivity vs mean(Diversity) 4x4 l_{tv}=3 (Near 0)']);
% subplot(2,2,4);
% plot(dataout2(:,1), dataout2(:,2), 'oblue',dataout2(:,1),dataout2(:,3),...
%        '-red');xlabel('Lyapunov Exponent'); ylabel('Variance of KLD Distances');
% title('(d)');
% axis([-3 3 0 .01]);
%title(['Sensitivity vs var(Diversity) 4x4 l_{tv}=3 (Near 0)']);
% subplot(2,2,3);
% plot(dataout3(:,1), dataout3(:,2), 'oblue',dataout3(:,1),dataout3(:,3),...
%        '-red');xlabel('Lyapunov Exponent'); ylabel('Mean KLD Distances');
%    
% title(['Sensitivity vs mean(Diversity) Uniformly Distributed']);
% subplot(2,2,4);
% plot(dataout4(:,1), dataout4(:,2), 'oblue',dataout4(:,1),dataout4(:,3),...
%        '-red');xlabel('Lyapunov Exponent'); ylabel('Variance of KLD Distances');
%    
% title(['Sensitivity vs var(Diversity) Uniformly Distributed']);

% %% Boxplots
% 
% figure
% subplot(1, 2, 1);
% hold on
% A = bplot(emptyList(:,2), 2, 'outliers'); 
% B = bplot(simresults(:,6), 1, 'outliers'); 
% legend(A, 'what');
% xlabel('1 = All Experiments, 2 = Uniformly Distributed'); ylabel('Mean KLD Distances');
% title(['Boxplot of Mean KLD Distances']);
% hold off;
% 
% subplot(1, 2, 2);
% hold on
% A = bplot(emptyList(:,3), 2, 'outliers'); 
% B = bplot(simresults(:,7), 1, 'outliers'); 
% xlabel('1 = All Experiments, 2 = Uniformly Distributed'); ylabel('Variance of KLD Distances');
% legend(A, 'what'); 
% title(['Boxplot of Variance of KLD Distances']);
% hold off;

%% Visual Steady State Distributions of Cells
%
load('/home/charlestreykang/Desktop/MATLAB/rbn-research-chris/results/gbatch_final/10by10output/10by10threshold2/outlier02_288.mat')
RBNp_ssDist = ssDist(RBNp);
figure(1);
for i = 1:9
    if (i <= 3)
        j = i;     
        if (mod(j,2) == 1)
            col = 'g';
        else
            col =  'b';
        end
     elseif (i <=6)
        j = i+7;
        if (mod(j,2) == 1)
            col = 'b';
        else
            col =  'g';
        end 
     else
        j = i+14;
        if (mod(j,2) == 1)
            col = 'g';
        else
            col = 'b';
        end
    end  
     %j=i;
     subplot(3,3, i)
     hold on
     set(gca, 'YScale', 'log')
     s1 = stem(1:length(RBNp_ssDist(j,:)), RBNp_ssDist(j,:),col)
     str=sprintf('Cell Number: %d', j);
     title(str);axis([0 1024 0 0.04]);
     hold off
     alpha(s1,.5);
end