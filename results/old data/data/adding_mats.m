clear all; clc;
a = load('simres1215.mat');
b = load('simres1235.mat');
c = load('simres1255.mat');
d = load('simres3415.mat');
e = load('simres5615.mat');
f = load('simres7815.mat');
g = load('simres9015.mat');

dat = [a.simresults; b.simresults; c.simresults; d.simresults; e.simresults; f.simresults; g.simresults]

scatter(dat(:,3), dat(:,4), 40, 'MarkerEdgeColor',[0 .5 .5], 'MarkerFaceColor', [0 .7 .7]); title('Lyapunov Exponent vs KLD Diversity'); xlabel('{\lambda}'); ylabel('mean KLD')