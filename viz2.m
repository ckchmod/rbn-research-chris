clear all; clc;
a = openfig('ConvergenceSSD600-v2.fig')
xlabel('Timestep $(t_n = 2^n)$', 'Interpreter', 'latex', 'FontSize', 32)
ylabel('$D_{KLS}(P_i(t_n)||P_i(t_{n+1}))$','Interpreter','latex', 'FontSize',32)

%b = openfig('DiversityLyapunov600.fig')
