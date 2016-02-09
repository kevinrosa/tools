function A = corc(xm,a,b,V,E)

% Calculate A matrix
% Call as 
%	A = corc(xm,a,b,V,E)
%
% where
%	xm = measurement locations
%	a, b = correlation function length scales
%	V,E = variances for phi

aa = a*a;
[n,m] = size(xm);
if n < m
	xa = xm'*ones(1,m);
else
   xa = xm*ones(1,n);
end

Ax = xa-xa';

A = cos(Ax/b).*exp(-(Ax.*Ax)/aa);
[n,m] = size(A);

A = V*A + E*diag(ones(n,1));
