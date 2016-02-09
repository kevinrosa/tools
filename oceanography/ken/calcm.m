function C = calcm(xm,x,a,b,V)

% Calculate C matrix
% Call as
%	C = calcm(xm,x,a,b,V)

aa = a*a;
[n,m] = size(xm);
[nn,mm] = size(x);

if n > m
	xm = xm';	
	[n,m] = size(xm);
end

if nn > mm
	x = x';
	[nn,mm] = size(x);
end

Ax = (xm'*ones(1,mm))';
Axx = x'*ones(1,m);
Ax  = Ax-Axx;

C = V*cos(Ax/b).*exp(-(Ax.*Ax)/aa);