function [x, fvec, f, check] = newt(x,fvecfunc)
% Globally convergent Newton Raphson Method
% given starting point x and fvecfunc
% return new zeros of fvec at new x with high accuracy, meanwhile check if
% is in a valley
MAXITS = 200;
TOLF = 1E-9; TOLMIN = 1E-12; STPMX = 100;
TOLX = eps;
[f, fvec] = fmin(x,fvecfunc);
test = max(abs(fvec));
% disp(['iteration: 0'])
% disp(['error = ', num2str(test)] )
if test < 0.01*TOLF
    check = false;
    disp('convergent!')
    return
end
stpmax = STPMX*max(norm(x),length(x));
for its = 1:MAXITS
    fjac = NRfdjac(x,fvecfunc);
    g = fvec.'*fjac;          % g is a row vector  
    xold = x;
    f_old = f;
    p = -inv(fjac)*fvec;
    [x, f,fvec, check] = lnsrch(xold, f_old, g, p, stpmax, @fmin, fvecfunc);
    test = max(abs(fvec));
%     disp(['iteration:',num2str(its)])
%     disp(['error = ', num2str(test)] )
    if (test < TOLF)
        check = false;
        disp('convergent!')
        return
    end
    if check
        den = max(f, 0.5*length(x));
        x_bar = abs(x); x_bar(x_bar<1) = 1;
        temp = abs(g).*x_bar/den;
        test = max(temp);
        check = (test < TOLMIN);
        disp('convergent!')
        return
    end
    x_bar = abs(x); x_bar(x_bar<1) = 1;
    temp = abs(x - xold)./x_bar;
    test = max(temp);
    if test < TOLX
        return
    end   
end
error('MAXITS exceed in newt')
end

