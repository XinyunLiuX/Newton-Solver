N = 63;
B = Amat(zeros(2*N+1,1));
B = -B(1:end-1,1:end-1);
[V, U] = eig(B);
plot(diag(U),'.')