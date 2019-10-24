% figure(1)
% openfig('single_cell_0_02_288.fig')
% figure(2)
% openfig('single_cell_1_02_288.fig')
% figure(3)
% openfig('single_cell_loop_02_288.fig')

a=hgload('single_cell_0_02_288.fig');
b=hgload('single_cell_1_02_288.fig');
c=hgload('single_cell_loop_02_288.fig');
% Prepare subplots
figure
h(1)=subplot(1,3,1);
h(2)=subplot(1,3,2);
h(3)=subplot(1,3,3);
% Paste figures on the subplots
copyobj(allchild(get(a,'CurrentAxes')),h(1));
copyobj(allchild(get(b,'CurrentAxes')),h(2));
copyobj(allchild(get(c,'CurrentAxes')),h(3));