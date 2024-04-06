function J = NRfdjac(x0,f)
% Newton Raphson finite difference Jacobian
% 《Numerical recipes: the art of scientific computing》chapter 9.7, p482
% Input: vector x0 (should be column vector)
% Output: Numerical approximation of Jacobian Matrix J of function f (results should be column vector)
% f: vector function defined elsewhere, could be easily modified
% N: Matrix size
% EPS: square root of machine precision; try default variable: eps
% h: step size for x0
% xh: a mall step from x0; xh = x0 + h*ej; Reset Caution!
% temp: storage for x0(j)
EPS = 1e-8;
n = length(x0);
J = zeros(n,n);
xh = x0;
for j = 1:n
    temp = x0(j);
    h = max(abs(temp)*EPS,EPS);
    xh(j) = temp + h;
    h = xh(j) - temp;           % Trick to reduce finite finite-precision error
    J(:,j) = (f(xh) - f(x0))/h; % Forward difference formular
    xh(j) = temp;               % Reset
end
end

