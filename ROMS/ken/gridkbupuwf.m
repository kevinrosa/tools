


function gridkbupuwf(filename)

% construct grid for xshelf channel geometry, using the new netcdf
%   Arrangement with inflow and outflow ports    
%
%   Call as
%       gridkbupuwf(filename)
%   where
%       filename is the name of the grid file, e.g., 'gridkb01.nc'

%	originally by Jim Lerczak, 2005
%	Modified by K. Brink, 2/21/2006 6/29/2007 1/15/2010 1/2011
%          2/7/2014

% revised to modern netcdf

figure(1)



%3D
Lm = 800;    % x points
Mm = 300;    % y points

%iflat = 3  %   # of offfshore grid points with flattened bottom (iflat >= 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid parameters.

L = Lm+1;
M = Mm+1;
Lp= L +1;
Mp= M +1;

% grid sizes in km

dy_max = 0.15;
dy_min = 0.15;
dx_min = 0.15;
dx_max = 0.25;


%	You place the river location by chosing a y grid point

% opening on  east side 
ymlow = 900*1000;
%ymlow = 960*1000;  %
ymhigh = 950*1000;


ymlowlow = 250*1000;
ymhighlow = 300*1000;
ymlowlow = 350*1000;    %

jrm = 1;
ishore = 1;

%	Scales for tanh (grid #)
tanhxsp = 300;
tanhxsm = 300;
tanhysp = 30;
tanhysm = 30;

%	Depth information 

%hcoast = 1;		% Depth at coast (m)
%alph = 0.005;		% Bottom slope

%rwidth = 3000.0;		% River total width (m)
%rdepth = 10;	% River depth (m)

%jshore = 3;
%bdepth = 100;  % basin depth
%gwidth = 10000; % depth gradient width in m

%   background topo information:

%   2-slope base model

x1 = 140.*1000.;
x2 = 155.*1000.;

hz = 5.0;    % depth at coast



%hz = 11.5;
h1 = 90;   % depth at shelf break
h2 = 200;  % cdepth in deep basin;



al1 = (h1-hz)/x1;
al2 = (h2 - h1)/(x2-x1);
%al1 = 0.00035;

%al1 = 0.0014;
%al1 = 0.00071
%al2 = 0;

%h1 = hz + al1*x1;

%hhh = h1;

%	Coriolis parameter (1/sec)

fzero = 0.95e-04;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct grid, and calculate grid metrics.

% Construct dx_r and dy_r.  Integrate to get x_r and y_r.
i=[1:Lp];
j=[1:Mp];
mp=(Lp+1)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Mine

ddy = (dy_min - dy_max);
jjp = find(j > jrm);
jjm = find(j < jrm);

dy_r(jjp) = dy_min - ddy*tanh((jjp-jrm)/tanhysp);
dy_r(jrm) = dy_min;
dy_r(jjm) = dy_min + ddy*tanh((jjm-jrm)/tanhysm);
dy_r = dy_r*1000;


ddx = (dx_min - dx_max);
iip = find(i > ishore);
iim = find(i < ishore);

dx_r(iip) = dx_min - ddx*tanh((iip-ishore)/tanhxsp);
dx_r(ishore) = dx_min;
dx_r(iim) = dx_min + ddy*tanh((iim-ishore)/tanhxsm);
dx_r = dx_r*1000;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Jim's:

%dx_r = (dx_max-(dx_min/2).*tanh((i-(mp-25))/15)+...
%    (dx_min/2).*tanh((i-(mp+25))/15)).*1000;
%dx_r = 0*i + dx*1000 ;
%dy_r = 0.5.*(dy_max-dy_min.*tanh((j-(mp+3))/2)+...
%    dy_min.*tanh((j-(mp-3))/2)).*1000;
%dy_r = 0*j + dy*1000 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Minimum dx = ',num2str(min(dx_r))]);
disp(['Minimum dy = ',num2str(min(dy_r))]);

% Calculate pm and pn
[pm,pn]=meshgrid(1./dx_r,1./dy_r);

% Calculate x_r and y_r
dx=[dx_r(1)./2 0.5.*(dx_r(1:end-1)+dx_r(2:end))];
dy=[dy_r(1)./2 0.5.*(dy_r(1:end-1)+dy_r(2:end))];
[x_r,y_r]=meshgrid(cumsum(dx),cumsum(dy));

% Shift grid so x_u(:,1)=0 and y_v(1,:)=0.
y_r = y_r - y_r(1,1) - (y_r(2,1)-y_r(1,1))/2;
x_r = x_r - x_r(1,1) - (x_r(1,2)-x_r(1,1))/2;

% Calculate dmde and dndx.
dndx = zeros(Mp,Lp);
dmde = zeros(Mp,Lp);

dndx(2:M,2:L) = (1./pn(2:M,3:Lp) - 1./(pn(2:M,1:Lm)))/2;
dmde(2:M,2:L) = (1./pm(3:Mp,2:L) - 1./(pm(1:Mm,2:L)))/2;


% Calculate x_u, etc.
x_u = (x_r(:,1:L) + x_r(:,2:Lp))/2;
y_u = (y_r(:,1:L) + y_r(:,2:Lp))/2;

x_v = (x_r(1:M,:) + x_r(2:Mp,:))/2;
y_v = (y_r(1:M,:) + y_r(2:Mp,:))/2;

x_p = (x_r(1:M,1:L) + x_r(2:Mp,2:Lp))/2;
y_p = (y_r(1:M,1:L) + y_r(2:Mp,2:Lp))/2;

el = y_u(end,1);
xl = x_v(1,end);




plot(x_r/1000,y_r/1000,'+')
title('Unmasked grid point locations')
xlabel('x grid location (km)')
ylabel('y grid location (km)')
set(gca,'dataaspectratio',[1 1 1])
ymax = max(max(y_r));
xmax = max(max(x_r));
ylow = -dy_max;
yhigh = ymax/1000 + dy_max;
xlow = -dx_max;
xhigh = xmax/1000 + dx_max;
axis([xlow xhigh ylow yhigh])

pause



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct grid topography.


% Jim's:

ymp=mean(y_r);  ymp=ones(Mp,1)*ymp;
xmp=mean(x_r')';xmp=xmp*ones(1,Lp);

%h = 15.0 - 10*(y_r./ymp-1).^2 - 5.0.*exp(-(x_r-xmp).^2/2e8)...
%    -10.0.*exp(-(x_r-xl).^2/2e8);
%h = 15.0 - 10*(y_r./ymp-1).^2 - 2.0.*exp(-(x_r-xmp).^2/2e8);
%h = 15.0 - 10*(y_r./ymp-1).^2 ;
%h(1 ,:)=h(2,:);
%h(Mp,:)=h(M,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	My topographies

%h = 0*x_r + rdepth;
%h = 0*x_r + rdepth;
h = 0*x_r;


%h = h + h2;

%yrm = y_r(jrm,1)
%xshore = x_r(1,ishore)
%bwidth = max(x_r(1,:))
%xshore = bwidth;
%yshore = y_r(jshore,1)
%yflat = yshore + gwidth;
%channel_len = max(max(y_r))-yshore;


%Iocean = find( y_r >= yshore);

h = h + h2;
Ishelf = find(x_r < x1);
h(Ishelf) = hz + x_r(Ishelf)*al1;



Islope = find(x_r >= x1 & x_r <= x2);
h(Islope) = h1 + (x_r(Islope) - x1)*al2;

%h(Ishelf) = hhh;
%h(Islope) = h2;


%Iland = find(x_r <= (bwidth-rwidth) & y_r < yshore);
%h(Iocean) = bdepth;
%h(Iland) = 0.0;
%yrmin = yrm - rhwidth;
%yrmax = yrm + rhwidth;

%Jriv = find(y_r < yshore & x_r > (bwidth-rwidth));
%h(Jriv) = rdepth;

%jjj = find(y_r(:,1) >= yshore & y_r(:,1) <= yflat);
%h(jjj,:) = ((bdepth-rdepth)*(y_r(jjj,:) - yshore)/gwidth) + rdepth;

%Iriv = find((y_r >= yrmin) & ... 
%	(y_r <= yrmax) & ...
%	(x_r >= xshore));

%h(Iriv) = rdepth;

[nrr,mrr] = size(h);
jj = find(0*h(:,mrr) == 0);
nj = length(jj);

%for ji = 1:nj
%    iqq = find(h(jj(ji),:) < rdepth);
%    if isempty( iqq) ~= 1
%        h(jj(ji),iqq) = rdepth;
%    end
%end




% Flatten out the bottom at offshore points
%if iflat > 0.5
%    hhold = h(:,(iflat+1));
%    for jjj = 1:iflat
%        h(:,jjj) = hhold;
%    end
%end

cg = sqrt(9.8*h);
[nq,mq] = size(x_r);



cgy = y_r./cg;
dxx = x_r(:,2:mq)-x_r(:,1:(mq-1)); 
xx = dxx./(0.5*(cg(:,2:mq)+cg(:,1:(mq-1))));
dyy = y_r(2:nq,:)-y_r(1:(nq-1),:);
yy = dyy./(0.5*(cg(2:nq,:)+cg(1:(nq-1),:)));


dtx = min(min(xx));
dty = min(min(yy));


xmax = max(x_r(1,:));
eedge = xmax - 3000;
%eedge = xmax - 300;
%eedge = 75000;
%eedge = 69000;  % 6/3/2010
%rmask(Imask) = ones(size(h(Imask)));
rmask = ones(size(x_r));
Imasklow = find(y_r < ymlow & x_r > eedge);
rmask(Imasklow) = 0*ones(size(Imasklow));
Imaskhigh = find(y_r > ymhigh & x_r > eedge);
rmask(Imaskhigh) = 0*ones(size(Imaskhigh));

Imaskg = find(y_r > ymlowlow & y_r < ymhighlow & x_r > eedge);
rmask(Imaskg) = ones(size(Imaskg));

%
rmask = 0.0*rmask;  %remove mask
%
%rmask = ones(size(x_r));   %
[Iwet] = find(0*h == 0);

[C,hh] = contour(x_r(1,:)/1000,y_r(:,1)/1000,h);
clabel(C,hh)
hold on
contour(x_r(1,:)/1000, y_r(:,1)/1000,rmask,'k');
hold off
set(gca,'dataaspectratio',[1 1 1])


title('Contoured topography, with land mask')

pause



pcolor(x_r(1,:)/1000,y_r(:,1)/1000,h)
set(gca,'dataaspectratio',[1 1 1])
shading flat
colorbar
hold on
contour(x_r(1,:)/1000, y_r(:,1)/1000,rmask,'k');
hold off

pause


dx=ones(Mp,1)*dx_r;
dy=dy_r'*ones(1,Lp);
hbar=sum(sum(h(Iwet).*dx(Iwet).*dy(Iwet)))./ ...
	sum(sum(dx(Iwet).*dy(Iwet)));
disp(['Hbar = ',num2str(hbar)]);
%disp(['Along channel gravity wave time = ',...
%	num2str((channel_len/sqrt(9.8.*hbar))/86400,'%5.8f'), ...%
%	' days phase shift']);
%disp(['Western cross-sectional area = ',...
%	num2str(sum(h(2:end-1,1).*dy(2:end-1,1)))]);
%disp(['Eastern cross-sectional area = ',...
%	num2str(sum(h(2:end-1,end).*dy(2:end-1,end)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create masking grids, Coriolis grid, and angle.

% Masking at RHO-points.
%rmask=0*ones(size(x_r));

Imask = find(0*h == 0);




for i=2:Lp,
  for j=1:Mp,
    umask(j,i-1)=rmask(j,i)*rmask(j,i-1);
  end,
end,
for i=1:Lp,
  for j=2:Mp,
    vmask(j-1,i)=rmask(j,i)*rmask(j-1,i);
  end,
end,
for i=2:Lp,
  for j=2:Mp,
    pmask(j-1,i-1)=rmask(j,i)*rmask(j,i-1)*rmask(j-1,i)*rmask(j-1,i-1);
  end,
end,

% Coriolis parameter.

f= fzero*ones(size(x_r));

% Angle of the grid is zero.
angle=zeros(size(x_r));


plot(rmask.*x_r,rmask.*y_r,'+')


%
iim = find(0*h ~= 0);
%h(iim) = hcoast;
%h(iim) = rdepth;
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create NetCDF file.

% Open file

ncc = netcdf.create(filename,'SHARE');
%nc=netcdf(filename,'clobber');
%nc.Description = 'Testing shelf-estuary grid';
%nc.Author = 'Ken Brink';
%nc.Created = datestr(now);
%nc.type = 'ESTUARY GRD file';

% Dimensions
dimidxr = netcdf.defDim(ncc,'xi_rho',Lp);
dimidxu = netcdf.defDim(ncc,'xi_u',L);
dimidxv = netcdf.defDim(ncc,'xi_v',Lp);
dimidxps = netcdf.defDim(ncc,'xi_psi',L);

dimider = netcdf.defDim(ncc,'eta_rho',Mp);
dimideu = netcdf.defDim(ncc,'eta_u',Mp);
dimidev = netcdf.defDim(ncc,'eta_v',M);
dimideps = netcdf.defDim(ncc,'eta_psi',M);

dimido = netcdf.defDim(ncc,'one',1);

%nc('xi_rho')=Lp;
%nc('xi_u')  =L;
%nc('xi_v')  =Lp;
%nc('xi_psi')=L;

%nc('eta_rho')=Mp;
%nc('eta_u')  =Mp;
%nc('eta_v')  =M;
%nc('eta_psi')=M;

%nc('one')=1;

% Create variables

varidxr = netcdf.defVar(ncc,'x_rho','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varidxr,x_r')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'x_rho'}= ncdouble(dims);
%nc{'x_rho'}(:,:)=x_r;

varides = netcdf.defVar(ncc,'x_psi','double',[dimidxps,dimideps]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varides,x_p')
netcdf.reDef(ncc);
%dims = {'eta_psi'; 'xi_psi'};
%nc{'x_psi'}= ncdouble(dims);
%nc{'x_psi'}(:,:)=x_p;

varidxu = netcdf.defVar(ncc,'x_u','double',[dimidxu,dimideu]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varidxu,x_u')
netcdf.reDef(ncc);
%dims = {'eta_u'; 'xi_u'};
%nc{'x_u'}= ncdouble(dims);
%nc{'x_u'}(:,:)=x_u;

varid = netcdf.defVar(ncc,'x_v','double',[dimidxv,dimidev]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,x_v')
netcdf.reDef(ncc);
%dims = {'eta_v'; 'xi_v'};
%nc{'x_v'}= ncdouble(dims);
%nc{'x_v'}(:,:)=x_v;

varid = netcdf.defVar(ncc,'y_rho','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,y_r')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'y_rho'}= ncdouble(dims);
%nc{'y_rho'}(:,:)=y_r;


varid = netcdf.defVar(ncc,'y_psi','double',[dimidxps,dimideps]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,y_p')
netcdf.reDef(ncc);
%dims = {'eta_psi'; 'xi_psi'};
%nc{'y_psi'}= ncdouble(dims);
%nc{'y_psi'}(:,:)=y_p;

varid = netcdf.defVar(ncc,'y_u','double',[dimidxu,dimideu]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,y_u')
netcdf.reDef(ncc);
%dims = {'eta_u'; 'xi_u'};
%nc{'y_u'}= ncdouble(dims);
%nc{'y_u'}(:,:)=y_u;

varid = netcdf.defVar(ncc,'y_v','double',[dimidxv,dimidev]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,y_v')
netcdf.reDef(ncc);
%dims = {'eta_v'; 'xi_v'};
%nc{'y_v'}= ncdouble(dims);
%nc{'y_v'}(:,:)=y_v;

varid = netcdf.defVar(ncc,'pm','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,pm')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'pm'}= ncdouble(dims);
%nc{'pm'}(:,:)=pm;

varid = netcdf.defVar(ncc,'pn','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,pn')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'pn'}= ncdouble(dims);
%nc{'pn'}(:,:)=pn;

varid = netcdf.defVar(ncc,'dmde','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,dmde')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'dmde'}= ncdouble(dims);
%nc{'dmde'}(:,:)=dmde;

varid = netcdf.defVar(ncc,'dndx','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,dndx')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'dndx'}= ncdouble(dims);
%nc{'dndx'}(:,:)=dndx;

varid = netcdf.defVar(ncc,'angle','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,angle')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'angle'}= ncdouble(dims);
%nc{'angle'}(:,:)=angle;

varid = netcdf.defVar(ncc,'mask_rho','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,rmask')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'mask_rho'}= ncdouble(dims);
%nc{'mask_rho'}(:,:)=rmask;

varid = netcdf.defVar(ncc,'mask_psi','double',[dimidxps,dimideps]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,pmask')
netcdf.reDef(ncc);
%dims = {'eta_psi'; 'xi_psi'};
%nc{'mask_psi'}= ncdouble(dims);
%nc{'mask_psi'}(:,:)=pmask;

varid = netcdf.defVar(ncc,'mask_u','double',[dimidxu,dimideu]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,umask')
netcdf.reDef(ncc);
%dims = {'eta_u'; 'xi_u'};
%nc{'mask_u'}= ncdouble(dims);
%nc{'mask_u'}(:,:)=umask;

varid = netcdf.defVar(ncc,'mask_v','double',[dimidxv,dimidev]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,vmask')
netcdf.reDef(ncc);
%dims = {'eta_v'; 'xi_v'};
%nc{'mask_v'}= ncdouble(dims);
%nc{'mask_v'}(:,:)=vmask;

varid = netcdf.defVar(ncc,'h','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,h')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'h'}= ncdouble(dims);
%nc{'h'}(:,:)=h;

varid = netcdf.defVar(ncc,'f','double',[dimidxr,dimider]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,f')
netcdf.reDef(ncc);
%dims = {'eta_rho'; 'xi_rho'};
%nc{'f'}= ncdouble(dims);
%nc{'f'}(:,:)=f;

varid = netcdf.defVar(ncc,'el','double',[dimido]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,el)
netcdf.reDef(ncc);
%dims = {'one'};
%nc{'el'} = ncdouble(dims);
%nc{'el'}(:) = el;

varid = netcdf.defVar(ncc,'xl','double',[dimido]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,xl)
netcdf.reDef(ncc);
%dims = {'one'};
%nc{'xl'} = ncdouble(dims);
%nc{'xl'}(:) = xl;

varid = netcdf.defVar(ncc,'spherical','char',[dimido]);
netcdf.endDef(ncc)
netcdf.putVar(ncc,varid,'F')
netcdf.reDef(ncc);
%dims = {'one'};
%nc{'spherical'} = ncchar(dims);
%nc{'spherical'}(:) = 'F';


%   Set globals and close out

ti = 'Model grid topography, unsmoothed ' ;
au = 'K. Brink';
cr = datestr(now);
ty = 'Grid File';
netcdf.putAtt(ncc,netcdf.getConstant('NC_GLOBAL'),'Title',ti);
netcdf.putAtt(ncc,netcdf.getConstant('NC_GLOBAL'),'Name',filename);
netcdf.putAtt(ncc,netcdf.getConstant('NC_GLOBAL'),'Author',au);
netcdf.putAtt(ncc,netcdf.getConstant('NC_GLOBAL'),'Created',cr);
netcdf.putAtt(ncc,netcdf.getConstant('NC_GLOBAL'),'Type',ty);


netcdf.close(ncc)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting

figure(1);

clf
%rmask(find(rmask==0)) = nan;
pcolor (x_r(2:end,2:end)./1000,y_r(2:end,2:end)./1000,...
        h(2:end,2:end).*rmask(2:end,2:end))
hold on
contour(x_r(1,:)/1000, y_r(:,1)/1000,rmask,'k');
hold off

set(gca,'dataaspectratio',[1 1 1])
hc=colorbar;hcl=get(hc,'ylabel');
set(hcl,'string','Depth (m)');
xlabel ('xi distance (km)');
ylabel ('eta distance (km)');
title ('Topography');


%return

figure (3)
plot (y_r(2:end-1,:),h(2:end-1,:),'-o');

figure(2);clf
subplot(2,1,1);
plot (x_r./1000,dx_r./1000);hold on;
plot (x_r./1000,dx_r./1000,'o');
xlabel('Cross channel distance (km)');
ylabel('Cross channel grid spacing (km)');

subplot(2,1,2);
plot (y_r./1000,dy_r./1000);hold on;
plot (y_r./1000,dy_r./1000,'o');
xlabel('Along channel distance (km)');
ylabel('Along channel grid spacing (km)');
