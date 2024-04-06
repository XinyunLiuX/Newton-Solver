function A = Amat(x0)
    global m e g
    sideband = m*ones(length(x0)-1,1);
	sideband(2:2:end) = e;
	A = diag(x0(end) - g*x0.^2) + diag(sideband,-1) + diag(sideband,1);
	A(end,:) = x0.';
	A(:,end) = 0;
end

