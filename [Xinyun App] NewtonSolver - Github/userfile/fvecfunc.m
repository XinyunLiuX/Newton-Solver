function fvec = fvecfunc(x0)
	fvec = Amat(x0)*x0 - bvec(x0);
end
