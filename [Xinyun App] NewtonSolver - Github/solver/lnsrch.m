function [x,f,fvec,check] = lnsrch(xold, f_old, g, p, stpmax, func, fvecfunc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Line Search and Backtracking: 
% given xold, fold, J, p, func = fmin = 1/2|F(x)|^2， fvecfunc = F(x) 
% return a new point x and fvec along the direciton p from xold where the function f has decrease sufficiently.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xold: initial n-dimentional point; fold: initial scaler f(xold); p: newton step at xold
% f(x) = 1/2F(x)*F(x) := func(x); G(lambda) = f(xold + lambda*p); fold = G(0); f = G(lambda_1); f2 = G(lambda_2)
% g: vector gradient of f(x) at xold; g = J*F(xold) obtained elsewhere
% alam: lambda_1; alam2: lambda_2; tmplam: latest lambda
% stpmax: maximum step of p
% slope = G'(lambda),satisfy G'(0) < 0 for newton direction is always a direction of descent for infinitesimal step
% ALF: alpha; eps = epsilon
ALF = 1E-4; TOLX = eps;
check = false;
summ = norm(p);
if summ > stpmax          
    p = p*stpmax/summ;  % Scale if attempted step is too big
end
slope = g*p;    % g is a row vector; initialization =  G'(0) 
if  slope >= 0 
    error('Roundoff problem in lnsrch'); % newton direction is always a direction of descent for inititesimal step
end
%% Compute lambda_min
xold_bar = abs(xold); xold_bar(xold_bar<1) = 1;
temp = abs(p)./xold_bar;
test = max(temp);
alamin = TOLX/test;
alam = 1.;                     % Always try Newton step first
%% Loop
while 1
    x = xold + alam*p;
    [f, fvec] = func(x,fvecfunc);    % undate fvec to output
    %% check return condition
    if alam < alamin                 % Return Condition 1: Convergence on delta x
        x = xold;
        check = true;
        return
    elseif f <= f_old + ALF*alam*slope % Return Condition 2: Sufficient function decrease
        return
    %% do something
    else                              % Subsequent backtracks
        if alam == 1                        % First backtrack
            tmplam = -slope/(2*(f - f_old - slope));
        else                                % Later backtrack
            rhs1 = f - f_old - alam*slope;
            rhs2 = f2 - f_old - alam2*slope;
            a = (rhs1/(alam*alam) - rhs2/(alam2*alam2))/(alam - alam2);
            b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2)/(alam - alam2));
            if a == 0
                tmplam = -slope/(2*b);
            else
                disc = b*b - 3*a*slope;
                if disc < 0
                    tmplam = 0.5*alam;
                elseif b <= 0
                        tmplam = (-b+sqrt(disc))/(3*a);
                else
                    tmplam = -slope/(b+sqrt(disc));
                end
            end
            if tmplam > 0.5*alam
                tmplam = 0.5*alam;            % lambda < 0.5*lambda1
            end
        end
    end
    %% update
    alam2 = alam;
    f2 = f;
    alam = max(tmplam, 0.1*alam);     % lambda > 0.1*lambda1                                               
end

end

