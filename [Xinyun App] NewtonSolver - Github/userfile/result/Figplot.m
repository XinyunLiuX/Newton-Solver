clc; clear all
file = dir(fullfile('*.mat'));
fname = {file.name};
fnameNum = length(fname);
w_data = zeros(fnameNum,1);
g_data = zeros(fnameNum,1);
for i = 1:fnameNum
    load(fname{i});
    w_data(i) = w;
    g_data(i) = g;
end
plot(g_data,w_data,'.','Markersize',15)
hold on;
plot(-8:0.1:8,m + (-8:0.1:8)/2, 'r','LineWidth',1.2)
plot(-8:0.1:8,-m +(-8:0.1:8)/2, 'r','LineWidth',1.2)
axis([-3,3,-3,3])
xticks([-3 0 3])
yticks([-3 0 3])
plot([-10,10],[1.5,1.5],'k--')
plot([-10,10],[-1.5,-1.5],'k--')
