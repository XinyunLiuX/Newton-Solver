% <<<Xinyun Notes>>> 
% must choose odd number
N = 63;
x0 = zeros(N*2+1,1);
%x_plus  = (sqrt(1 - 2*m/g) + sqrt(1 + 2*m/g))/2;
%x_minus = (sqrt(1 - 2*m/g) - sqrt(1 + 2*m/g))/2;
x0(N) = 1/sqrt(2);
x0(N+1) = -1/sqrt(2);
x0(end) = m + g/2;
