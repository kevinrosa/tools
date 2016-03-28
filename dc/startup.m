fprintf('\n startup.m entered \n');

% update toolboxcache
rehash toolboxcache

% unicode support
feature('DefaultCharacterSet', 'UTF8');

warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
warning('off','MATLAB:lang:cannotClearExecutingFunction');

%% Make good figures

% good colormap
set(0,'DefaultFigureColormap',flipud(cbrewer('div', 'RdYlBu', 32)));
%set(0, 'DefaultFigureColormap', parula); close;
fontName = 'Helvetica';

% figure properties
%set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultFigureColor','w')
set(0,'DefaultFigureRenderer','opengl')
set(0,'DefaultFigurePaperPositionMode', 'auto');

set(0,'DefaultTextFontName', fontName);
set(0,'DefaultTextColor','k')
set(0,'DefaultTextFontSize',16);

% line style and color cycling
set(0,'DefaultAxesLineStyleOrder',{'-','--','-.'});
%set(0,'DefaultAxesColorOrder', cbrewer('qual', 'Dark2', 8));
set(0,'DefaultAxesColorOrder', linspecer(8, 'qualitative'));
%set(0,'DefaultAxesColorOrder', distinguishable_colors(15));

% axes
set(0,'DefaultAxesFontName',fontName)
set(0,'DefaultAxesFontWeight','normal')
set(0,'DefaultAxesTickLength'  , [.01 .01]);
set(0,'DefaultAxesFontSize',14)
set(0,'DefaultAxesBox','on')
set(0,'DefaultAxesTickDir','out')
%set(0,'DefaultAxesXMinorTick','on')
%set(0,'DefaultAxesYMinorTick','on')
%set(0,'DefaultAxesZMinorTick','on')
set(0,'DefaultAxesXColor',[.3 .3 .3])
set(0,'DefaultAxesYColor',[.3 .3 .3])
set(0,'DefaultAxesZColor',[.3 .3 .3])
set(0,'DefaultAxesLineWidth',1)
set(0,'DefaultLineLineWidth',2);

% multiple monitor support
%mm = get(0,'MonitorPositions');
%if size(mm,1) == 2,
%    disp('2 monitors detected.');
%    set(0,'DefaultFigurePosition',[0.725*(mm(1,3)+mm(2,3)) 0.5*mm(2,4) 560 420]);
%end


%% change to current working dir
if ~strcmpi(computer,'GLNXA64')
   cd('E:\Work\eddyshelf\');
else
   mach = evalc('system hostname');
   if strfind(mach,'poison')
     cd('/home/poison/deepak/ROMS/runs/eddyshelf/scripts/');
   else
     cd('/media/data/Work/eddyshelf/');
   end
end

