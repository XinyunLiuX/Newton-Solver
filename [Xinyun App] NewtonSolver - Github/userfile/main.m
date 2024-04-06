clc; clear all
addpath('../solver')
global m e g
m = 1;
for g = -0.6:-0.1:-3
x0_generate
% <<<Xinyun Notes>>> 
% for type-II breather, when g approaching 0, the step size for e should be set small.
% for example, 
% de = 0.001 for g = (+-)0.1;
% de = 0.002 for g = (+-)0.2;
% de = 0.005 for g = (+-)0.3;...
% Also, there might be a bifurcation at g = (+-)4, for type II(-+), as we can see from Fig 3(d)
for e = 0:0.01:0.15
    [x0, fvec, f, check] = newt(x0,@fvecfunc);
    disp(['\omega = ',num2str(x0(end))])
    if abs(e-0.15)<1E-8
        dir = 'result'; mkdir(dir);
        type = 'IIminus';
        filename = [dir,'/m=',num2str(m),'_e=',num2str(e),'_g=',num2str(g),type,'.mat'];
        w = x0(end);
        x_profile = x0(1:end-1);
        save(filename,'x_profile','g','m','e','w')
    end
end
bar(x0(1:end-1))
pause(2)
end