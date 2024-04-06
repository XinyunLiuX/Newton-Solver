function [f, fvec] = fmin(x, func)
% f = 1/2fvec'*fvec
fvec = func(x);
f = sum(fvec.^2)/2;
end

