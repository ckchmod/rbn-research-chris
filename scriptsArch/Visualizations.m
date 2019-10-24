%clear all; clc;
a = load('checkpoint12.mat'); b = load('checkpoint34.mat'); c = load('checkpoint56.mat'); d = load('checkpoint56test.mat')
e = load('checkpoint78.mat'); f = load('checkpoint78test.mat'); g = load('checkpoint90.mat');
simresults = cat(1, a.dataFile, b.dataFile);
simresults = cat(1, simresults, c.dataFile);
simresults = cat(1, simresults, d.dataFile);
simresults = cat(1, simresults, e.dataFile);
simresults = cat(1, simresults, f.dataFile);
simresults = cat(1, simresults, g.dataFile);


a = load('batch1'); b= load('batch2'); simresults = cat(1, a.simresults, b.simresults);
simul_iter=max(simresults(:,4));
figure(1)
% 2D Visualization: p
scatter(simresults(:,5),simresults(:,6), 'filled');
xlabel('Lyapunov Exponent'); ylabel('diversity - average KLD Distance');
title(['Sensitivity vs Diversity (each RBNp repeated ' num2str(simul_iter) ...
    ' times)']);

figure(2)
plot(convergePlot); xlabel('t'); ylabel('D_{KLsymmetric}<\epsilon'); 
title('convergence of a distribution k=1, p=0.05')

% Estimate SSD Converge Rate
plot(estimates'); xlabel('2^n Timesteps'); ylabel('KLD_{Symmetric}(P(2^n),P(2^{n+1}))'); title('Convergence of Steady-State Distributions'); legend('Maximum Chaotic: \lambda=1.4420', 'Critical: \lambda=0', 'Maximum Ordered: \lambda=-1.8981')
%

% hold on
% p = polyfit(simresults(:,5), simresults(:,6), 2);
% x1 = linspace(0, 2.5);
% f1 = polyval(p, x1);
% plot(x1, f1);
% legend('KLD Distances','2nd order-polyfit')
% hold off
% 
% 
% figure(2)
% % 3D Visualization: p vs k vs mean(KLD)
% scatter3(simresults(:,2), simresults(:,1), simresults(:,4), 'filled');
% xlabel('p'); ylabel('k'), zlabel('mean KLD');
% hold on
% d = linspace(0,1);
% fs = .5.*(1./(d.*(1-d)));
% z = zeros(1,length(d));
% zmin = min(simresults(:,4)); zmax = max(simresults(:,4));
% plot3(d, fs, z); axis([0 1 0 5 zmin zmax]);
% 
% a = linspace(0,1);
% b = linspace(0,.3);
% [X, Y] = meshgrid(a,b);
% Z = .5.*(1./(X.*(1-X)));
% 
% mesh(X,  Z, Y) 
% xlabel('p'); ylabel('k'); title('3D view of mean KLDs in terms of p and k')
% hold off
