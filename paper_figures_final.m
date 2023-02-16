%% loading data
% Winnie's figures

% TO DO: 
% Fig 1 - single colorbar
% Fig 3 - change aspect ratio (wider)

%Adding paths + loading stuff for plotting/saving figures
addpath /home/mmazloff/ANALYSIS  %this has imagescnan and rdmds
addpath /home/mmazloff/ANALYSIS/export_fig/ %this has export fig

cd /home/winnie.chu/data
load pCO2update.mat pCO2res
load DCall_1deg_update.mat DCL

load /data/SO6/TPOSE/fwdtp6_bgc/grid_6_k51/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC

cd /data/SO6/TPOSE_diags/fwdtp6_bgc/diags 
ts=1460:1460:275940;
nt=length(ts);

%Loading pCO2 data
pCO2 = zeros(1128,336,nt,'single'); 
pCO2= rdmds(['diag_surf'],ts,'rec',2);
pCO2(pCO2==0)=nan;
pCO2 = pCO2(:,:,13:180); % cropping first year (but make sure the number of total months is divisible by 12)
nmonths = size(pCO2,3);
% should be saved under cd /home/winnie.chu/data
% called pCO2.mat
% load pCO2.mat should be good and contains pCO2 and pCO2res
% save pCO2update.mat pCO2 pCO2res nmonths

% atm to muatm
pCO2=pCO2*1e6;

cd /data/SO6/TPOSE/projects/pco2_correlations

% climatology
load /data/SO6/TPOSE/validation/gridded_products/landschutzer2018_tpose_grid6
load /data/SO6/TPOSE/tpose6/grid_6/grid XC YC
clim = pco2_interp(:,:,265:end); % 2004+

%% make figure 1

figure; set(gcf,'position',[50 50 1100 700],'paperpositionmode','auto','Color','w');
tiles=tiledlayout(2,1);

% model pCO2
nexttile
tmean_pCO2=nanmean(pCO2,3);
[c h]=contourf(x,y,tmean_pCO2',[250:650]); %colorbar;
t=title('a)');
t.FontSize = 18;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none')
set(gca,'color',[.8 .8 .8],'fontsize',20)
caxis([300 520])
%text(104,22.5,'a)','fontsize',16,'fontweight','bold')
ylim([-20 20]);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});
set(gca,'xtick',[])

% climatology pCO2
nexttile;
clim_mean=nanmean(clim,3);
[c h]=contourf(x,y,clim_mean',[250:650]); cb=colorbar;
t=title('b)');
t.FontSize = 18;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none')
set(gca,'color',[.8 .8 .8],'fontsize',20)
caxis([300 520])
%text(104,22.5,'b)','fontsize',16,'fontweight','bold')
ylim([-20 20]);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});

% shared colorbar
cb.Layout.Tile = 'east';
cb.Title.String='\muatm';
tiles.TileSpacing = 'tight';

% save figure
cd /home/winnie.chu/figures
eval(['export_fig  Fig1_v6.png']);

%% downloading figures to local computer (use gitbash)
scp winnie.chu@aha.ucsd.edu:/home/winnie.chu/figures/Fig1_v4.png C:/Users/wuc11/OneDrive/Desktop/surf

%% make figure 2

% plot pCO2 vs mooring data

cd /data/SO6/TPOSE_diags/fwdtp6_bgc/diags
pCO2fwd = rdmds('diag_surf',1460:1460:275940,'rec',2);
pCO2fwd = pCO2fwd*1e6;
pCO2fwd(pCO2fwd<=0)=NaN;

tmodelfwd=datenum(2004,1:size(pCO2fwd,3),15);
tclim = datenum(2004,1:168,15);

cd /data/SO6/TPOSE/constraints/profiles_bgc

list_vars = {'DIC','ALK','O2','NO3','PO4','PH','PCO'};
list_profs = {'GLO','SOC','UWO2'};

iprof=2;
ivar=7;

files_list = dir('SOC*');
profname = 'SOCATv6';

prof_obs=[]; prof_estim=[]; prof_weight=[]; prof_lon=[]; prof_lat=[]; prof_YYYYMMDD=[]; prof_date=[]; prof_HHMMSS=[]; prof_descr=[];
if length(files_list)>0 % ok?
    for m=1:length(files_list)
        ncid = netcdf.open(files_list(m).name,'nowrite');
        varfound=1;
        try
            ID = netcdf.inqVarID(ncid,['prof_' char(list_vars{ivar})]);
        catch exception
            if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                varfound=0;
            end
        end
        netcdf.close(ncid)
        if varfound==1
            prof_obs = [prof_obs ncread(files_list(m).name,['prof_' char(list_vars{ivar})])];
            prof_weight = [prof_weight ncread(files_list(m).name,['prof_' char(list_vars{ivar}) 'weight'])];
            prof_lon = [prof_lon; ncread(files_list(m).name,'prof_lon')];
            prof_lat = [prof_lat; ncread(files_list(m).name,'prof_lat')];
            prof_YYYYMMDD = [prof_YYYYMMDD; ncread(files_list(m).name,'prof_YYYYMMDD')];
            prof_HHMMSS = [prof_HHMMSS; ncread(files_list(m).name,'prof_HHMMSS')];
            prof_date = [prof_date; ncread(files_list(m).name,'prof_date')];
            if size(ncread(files_list(m).name,'prof_descr'),1)==length(ncread(files_list(m).name,'prof_lon'))
                prof_descr = [prof_descr ncread(files_list(m).name,'prof_descr')'];
            elseif size(ncread(files_list(m).name,'prof_descr'),2)==length(ncread(files_list(m).name,'prof_lon'))
                prof_descr = [prof_descr ncread(files_list(m).name,'prof_descr')];
            end
        end % if varfound==1
    end
prof_depth = ncread(files_list(m).name,'prof_depth');
prof_lon(prof_lon<0) = prof_lon(prof_lon<0)+360;
end

prof_obs(prof_weight<=0) = NaN;
prof_obs(prof_lon<148 & prof_lat<-8) = NaN;
prof_obs(prof_lon<131) = NaN;
prof_obs(abs(prof_lat)>17) = NaN;
wt=prof_weight;

load /data/averdy/datasets/SOCAT/v6/TAOexpocodes

for n=1:length(prof_obs)
    descr(n) = cellstr(char(prof_descr(:,n)'));
end
moorings = 0; 
for n=1:length(TAOexpocodes)
    moorings = moorings + strcmp(descr,TAOexpocodes{n});
    ind = find(strcmp(descr,TAOexpocodes{n}));
    if length(ind)>0
        for k=1:length(ind)
            descr{ind(k)}=TAOmooring{n};
        end
    end
end
descr=descr(logical(moorings));
prof_lon=prof_lon(logical(moorings));
prof_lat=prof_lat(logical(moorings));
prof_date=prof_date(logical(moorings));
prof_obs=prof_obs(logical(moorings));
[descrU,IA,IC]=unique(descr);
  
load ~/tpose/grid_6/grid XC YC hFacC

cd /data/SO6/TPOSE/constraints/profiles

list_vars = {'T'};
ivar=1;

files_list = dir('NODC*');
profname = 'MRB';

prof_obsT=[]; prof_estimT=[]; prof_weightT=[]; prof_lonT=[]; prof_latT=[]; prof_YYYYMMDDT=[]; prof_dateT=[]; prof_HHMMSST=[]; prof_descrT=[];
if length(files_list)>0 % ok?
    for m=1:length(files_list)
        ncid = netcdf.open(files_list(m).name,'nowrite');
        varfound=1;
        try
            ID = netcdf.inqVarID(ncid,['prof_' char(list_vars{ivar})]);
        catch exception
            if strcmp(exception.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
                varfound=0;
            end
        end
        netcdf.close(ncid)
        if varfound==1
            prof_obsT = [prof_obsT ncread(files_list(m).name,['prof_' char(list_vars{ivar})])];
            prof_weightT = [prof_weightT ncread(files_list(m).name,['prof_' char(list_vars{ivar}) 'weight'])];
            prof_lonT = [prof_lonT; ncread(files_list(m).name,'prof_lon')];
            prof_latT = [prof_latT; ncread(files_list(m).name,'prof_lat')];
            prof_YYYYMMDDT = [prof_YYYYMMDDT; ncread(files_list(m).name,'prof_YYYYMMDD')];
            prof_HHMMSST = [prof_HHMMSST; ncread(files_list(m).name,'prof_HHMMSS')];
            prof_dateT = [prof_dateT; ncread(files_list(m).name,'prof_date')];
            if size(ncread(files_list(m).name,'prof_descr'),1)==length(ncread(files_list(m).name,'prof_lon'))
                prof_descrT = [prof_descrT ncread(files_list(m).name,'prof_descr')'];
            elseif size(ncread(files_list(m).name,'prof_descr'),2)==length(ncread(files_list(m).name,'prof_lon'))
                prof_descrT = [prof_descrT ncread(files_list(m).name,'prof_descr')];
            end
        end % if varfound==1
    end
prof_depthT = ncread(files_list(m).name,'prof_depth');
prof_lonT(prof_lonT<0) = prof_lonT(prof_lonT<0)+360;
end

prof_obsT(prof_weightT<=0) = NaN;
prof_obsT(prof_lonT<148 & prof_latT<-8) = NaN;
prof_obsT(prof_lonT<131) = NaN;
prof_obsT(abs(prof_latT)>17) = NaN;
wt=prof_weightT;

load /data/averdy/datasets/SOCAT/v6/TAOexpocodes

for n=1:length(prof_obsT)
    descrT(n) = cellstr(char(prof_descrT(13:17,n)'));
end
[descrUT,IA,IC]=unique(descrT);
 
for n=[14 15 17 21 53 54 55]
    ind=strcmp(descr,descrUT{n});
    lonmean(n)=mean(prof_lonT(ind)); 
    min(prof_lonT(ind));
    max(prof_lonT(ind));
    latmean(n)=mean(prof_latT(ind));
    min(prof_latT(ind));
    max(prof_latT(ind));
end

aa=find(abs(latmean)<0.5)
bb=find(abs(lonmean-236)<0.1)
bb=find(abs(lonmean-190)<0.1)


cd /data/SO6/TPOSE_diags/fwdtp6_bgc/diags
ts=1460:1460:249840;
for t=1:length(ts)
    tmp = rdmds('diag_state',ts(t),'rec',1);
    Tfwd5(:,:,t) = tmp(:,:,1);
    Tfwd5_20m(:,:,t) = mean(tmp(:,:,4:5),3);
end
Tfwd5(Tfwd5<=0)=NaN;
tmodelfwd5 = datenum(2004,1:171,15);

figure;
set(gcf,'position',[441   210   823   565],'paperpositionmode','auto','color','w');
clf

axes;
set(gca,'position',[.1 .55 .4 .4]);
n=7;
ind=strcmp(descr,descrU{n});

plot(prof_date(ind),prof_obs(ind)*1e6,'.','color',[.5 .5 .5]); 

lonm=mean(prof_lon(ind)); 
min(prof_lon(ind))
max(prof_lon(ind))

latm=mean(prof_lat(ind));
min(prof_lat(ind))
max(prof_lat(ind))

[~,ii]=min(abs(lonm-XC(:,1)));
[~,jj]=min(abs(latm-YC(1,:)));

hold on;plot(tclim,squeeze(clim(ii,jj,:)),'b','linewidth',2);
hold on;plot(tmodelfwd,squeeze(pCO2fwd(ii,jj,:)),'linewidth',2,'color',[.9 .5 .2]);

set(gca,'xtick',[],'fontsize',12,'xlim',[datenum(2004,1,1) datenum(2019,1,1)]);
ylabel('pCO_2 (\muatm)')

legend('TAO','Landschutzer','model')

ylim([300 600])

text(datenum(2004,1,1),620,['a) ' num2str(round(360-lonm)) '^\circW'],'fontsize',14,'fontweight','bold');

axes;
set(gca,'position',[.55 .55 .4 .4]);
n=2;
ind=strcmp(descr,descrU{n});

plot(prof_date(ind),prof_obs(ind)*1e6,'.','color',[.5 .5 .5]); 

lonm=mean(prof_lon(ind)); 
min(prof_lon(ind))
max(prof_lon(ind))

latm=mean(prof_lat(ind));
min(prof_lat(ind))
max(prof_lat(ind))

[~,ii]=min(abs(lonm-XC(:,1)));
[~,jj]=min(abs(latm-YC(1,:)));

hold on;plot(tclim,squeeze(clim(ii,jj,:)),'b','linewidth',2);
hold on;plot(tmodelfwd,squeeze(pCO2fwd(ii,jj,:)),'linewidth',2,'color',[.9 .5 .2]);

set(gca,'xtick',[],'ytick',[],'fontsize',12,'xlim',[datenum(2004,1,1) datenum(2019,1,1)]);

ylim([300 600])

text(datenum(2004,1,1),620,['b) ' num2str(round(360-lonm)) '^\circW'],'fontsize',14,'fontweight','bold');

n=14;
ind=strcmp(descrT,descrUT{n});

axes;
set(gca,'position',[.1 .08 .4 .4]);

plot(prof_dateT(ind),prof_obsT(2,ind),'.','color',[.5 .5 .5]); 
if n==15
    plot(prof_dateT(ind),prof_obsT(6,ind),'.','color',[.5 .5 .5]); 
end

lonm=mean(prof_lonT(ind)); 
min(prof_lonT(ind))
max(prof_lonT(ind))

latm=mean(prof_latT(ind));
min(prof_latT(ind))
max(prof_latT(ind))

[~,ii]=min(abs(lonm-XC(:,1)));
[~,jj]=min(abs(latm-YC(1,:)));

if n==15
    hold on;plot(tmodelfwd5,squeeze(Tfwd5_20m(ii,jj,:)),'linewidth',2,'color',[.9 .5 .2]);
else
    hold on;plot(tmodelfwd5,squeeze(Tfwd5(ii,jj,:)),'linewidth',2,'color',[.9 .5 .2]);
end

ylabel('T (^\circC)')
ylim([18 32])

set(gca,'fontsize',12,'xlim',[datenum(2004,1,1) datenum(2019,1,1)],'xtick',datenum(2006:4:2019,1,1),'xticklabel',2006:4:2020);

text(datenum(2004,1,1),33,['c) '],'fontsize',14,'fontweight','bold');

axes;
set(gca,'position',[.55 .08 .4 .4]);

n=15;
ind=strcmp(descrT,descrUT{n});

plot(prof_dateT(ind),prof_obsT(2,ind),'.','color',[.5 .5 .5]); 
if n==15
    plot(prof_dateT(ind),prof_obsT(6,ind),'.','color',[.5 .5 .5]); 
end
% hold on; plot(prof_estim(ind));

lonm=mean(prof_lonT(ind)); 
min(prof_lonT(ind))
max(prof_lonT(ind))

latm=mean(prof_latT(ind));
min(prof_latT(ind))
max(prof_latT(ind))

[~,ii]=min(abs(lonm-XC(:,1)));
[~,jj]=min(abs(latm-YC(1,:)));

if n==15
    hold on;plot(tmodelfwd5,squeeze(Tfwd5_20m(ii,jj,:)),'linewidth',2,'color',[.9 .5 .2]);
else
    hold on;plot(tmodelfwd5,squeeze(Tfwd5(ii,jj,:)),'linewidth',2,'color',[.9 .5 .2]);
end

datetick('x')
ylim([18 32])

set(gca,'ytick',[],'fontsize',12,'xlim',[datenum(2004,1,1) datenum(2019,1,1)],'xtick',datenum(2006:4:2019,1,1),'xticklabel',2006:4:2020);

text(datenum(2004,1,1),33,['d) '],'fontsize',14,'fontweight','bold');

cd /data/SO6/TPOSE/projects/pco2_correlations
eval(['export_fig  Fig2.png']); 

%%
figure; set(gcf,'position',[50 50 880 1200],'paperpositionmode','auto','Color','w');
tiledlayout(3,1);
nexttile;
imagesc(1,1,1); cb1=colorbar;
t=title('a)');
t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
nexttile;
imagesc(2,2,2); cb2=colorbar;
t=title('b)');
t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
nexttile;
imagesc(3,3,3); cb3=colorbar;
t=title('c)');
t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';

%% make figure 3

% load data
cd /data/SO6/TPOSE/projects/pco2_correlations
load data_fig3

figure; set(gcf,'position',[50 50 880 1200],'paperpositionmode','auto','Color','w');
tiledlayout(3,1);

% Plotting total variance in pCO2 (pre-processing)
nexttile
pCO2var = var(pCO2,0,3);
[c h]=contourf(x,y,pCO2var',[0:5:1200]); cb3a=colorbar; 
t=title('a)');
t.FontSize = 18;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8],'fontsize',18)
caxis([100 1000]);
ylim([-20 20]);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});
set(gca,'xtick',[])
cb3a.Title.String = '\muatm^2';

% Plotting percent variance explained by annual cycle
nexttile
[c h]=contourf(x,y,pCO2VE(:,:,1)',[0:100]);
t=title('b)');
t.FontSize = 18;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8],'fontsize',18)
ylim([-20 20]);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$E','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});
set(gca,'xtick',[])
cb3b=colorbar;
caxis([0 100]);
cb3b.Title.String = '%';

% Plotting percent residual variance 
nexttile
[c h]=contourf(x,y,pCO2VE(:,:,5)',[0:100]); %colorbar; 
t=title('c)');
t.FontSize = 18;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8],'fontsize',18)
ylim([-20 20]);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});
cb3c=colorbar;
caxis([0 100]);
cb3c.Title.String = '%';

colormap(nexttile(1),'default');
col=[[256,256,256]'/256,[75,84,202]'/256]'; % two points - first is white, second is blue but can be ANY two colors
col=interp1([0:1],col,[0:1/999:1]); % interps on to 100 point grid…or 1000 or whatever you want, 100 is usually more than enough! 
colormap(nexttile(2),col);
colormap(nexttile(3),col);

cd /home/winnie.chu/figures
eval(['export_fig Fig3_v3.png']);

%Finding and plotting total variance in pCO2res (fig. 3d which was not included)
pCO2resvar = var(pCO2res,0,3);
clf
imagescnan(x,y(34:277),pCO2resvar(:,34:277)','NanColor',[.8 .8 .8]);axis xy; a=colorbar; 
caxis([100 1000])
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$N','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
text(291,22.5,'\muatm^2','fontsize',14);
set(gca,'fontsize',16)
set(gcf,'Position',[50 50 1100 500],'Color', 'w');
cd /data/SO6/TPOSE/projects/pco2_correlations
%eval(['export_fig Fig3d.png']);


%Plotting residual variance (not included)
clf
[c h]=contourf(x,y,100-pCO2VE(:,:,5)',[0:100]); colorbar; 
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8])
ylim([-20 20]);
colormap 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
col=[[256,256,256]'/256,[75,84,202]'/256]'; % two points - first is white, second is blue but can be ANY two colors
col=interp1([0:1],col,[0:1/999:1]); % interps on to 100 point grid…or 1000 or whatever you want, 100 is usually more than enough! 
colormap(col);
a=colorbar;
%a.Label.String = 'Percent residual variance [%]';
%a.FontSize=16;
text(295,22.5,'%','fontsize',14);
text(104,22.5,'c)','fontsize',16,'fontweight','bold')
set(gca,'fontsize',16)

%% make figure 4

figure; set(gcf,'position',[50 50 900 350],'paperpositionmode','auto','Color','w');
cd /data/SO6/TPOSE/projects/pco2_correlations
load data_fig6

clf;
[c h]=contourf(X2(37:401),Y2(14:149),DC2(37:401,14:149,:)',-.1:.01:1); cb4=colorbar;
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8]);
caxis([0 1])
colormap('parula')
[testX1,testX2] = meshgrid(X2, Y2);
hold on; 
mjr=DCL(337,157,1); % for 0 N, 160 E
mnr=DCL(337,157,2); % for 0 N, 160 E
thta=DCL(337,157,3); % for 0 N, 160 E
F2=exp(-((testX1.*cos(thta)+testX2.*sin(thta)).^2/mjr.^2+(testX1.*sin(thta)+testX2.*cos(thta)).^2./mnr^2));
v=[1/exp(1),1/exp(1)];
contour(X2,Y2,F2,v,'black','LineWidth',1); %show the contour line w/ s-major length, s-minor length
xlabel('zonal distance from center (km)');
ylabel('meridional distance from center (km)');
set(gca,'fontsize',16)
xlim([-3000 3000])
ylim([-1250 1250])
cb4.Title.String = 'r';


eval(['export_fig Fig4_v3.png']);

%% make figure 5

load /home/winnie.chu/data/DCall_1deg

figure; set(gcf,'position',[50 50 880 1200],'paperpositionmode','auto','Color','w');
tiledlayout(3,1);

%Plotting semi-major lengths (replace 1 index with 2 or 3 or 4 for mnr or thta or r^2)
nexttile
[row tmp]=find(~isnan(DCL(:,:,1))); % row is row, col is linear index
row=unique(row);
[col d]=ind2sub([336,1],tmp); % turn linear index into column and dimension were d is just a bunch of 1s 
col=unique(col);
QQ2=DCL(row,col,1); % just to make sure we don’t lose the data
QQ2(QQ2<10)=nan; % bc there’s a bunch of random min that are 1.681 which
%imagescnan(x(row),y(col),real(QQ2')); axis xy; a=colorbar;
[c h]=contourf(x(row),y(col),real(QQ2'),1000:10:6000); cb5a=colorbar;
t=title('a) semi-major');
t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8],'fontsize',18)
ylim([-20 20]);
caxis([1.0 4.5]*1e3) % find based on extrema
%text(291,22.5,'km','fontsize',14);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});
set(gca,'xtick',[]);
cb5a.Title.String = 'km';
%text(104,22.5,'a) semi-major','fontsize',16,'fontweight','bold')
%eval(['export_fig  Fig7a.png']); 

nexttile
%Plotting semi-minor lengths (replace 1 index with 2 or 3 or 4 for mnr or thta or r^2)
[row tmp]=find(~isnan(DCL(:,:,2))); % row is row, col is linear index
row=unique(row);
[col d]=ind2sub([336,1],tmp); % turn linear index into column and dimension were d is just a bunch of 1s 
col=unique(col);
QQ2=DCL(row,col,2); % just to make sure we don’t lose the data
QQ2(QQ2<10)=nan; % bc there’s a bunch of random min that are 1.681 which
%imagescnan(x(row),y(col),real(QQ2')); axis xy; a=colorbar;
[c h]=contourf(x(row),y(col),real(QQ2'),0:10:2000); cb5b=colorbar;
t=title('b) semi-minor');
t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8],'fontsize',18)
ylim([-20 20]);
caxis([100 1000]) % find based on extrema
%text(291,22.5,'km','fontsize',14);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});
set(gca,'xtick',[]);
cb5b.Title.String = 'km';
%text(104,22.5,'b) semi-minor','fontsize',16,'fontweight','bold')

%eval(['export_fig  Fig7b.png']); 

%Plotting r^2 (replace 1 index with 2 or 3 or 4 for mnr or thta or r^2)
nexttile
[row tmp]=find(~isnan(DCL(:,:,4))); % row is row, col is linear index
row=unique(row);
[col d]=ind2sub([336,1],tmp); % turn linear index into column and dimension were d is just a bunch of 1s 
col=unique(col);
QQ2=DCL(row,col,4); % just to make sure we don’t lose the data
%imagescnan(x(row),y(col),real(QQ2')); axis xy; a=colorbar;
[c h]=contourf(x(row),y(col),(QQ2'),0:.001:1); cb5c=colorbar;
t=title('c) r^2');
t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
set(gca,'color',[.8 .8 .8],'fontsize',18)
ylim([-20 20]);
caxis([0 1]) % find based on extrema
%text(291,22.5,'km','fontsize',14);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$N'});
%text(104,22.4,'c) r^2','fontsize',16,'fontweight','bold')
cb5c.Title.String = 'r^2';

col=[[256,256,256]'/256,[75,84,202]'/256]'; % two points - first is white, second is blue but can be ANY two colors
col=interp1([0:1],col,[0:1/999:1]); % interps on to 100 point grid…or 1000 or whatever you want, 100 is usually more than enough! 
colormap(nexttile(3),col);

eval(['export_fig  Fig5.png']);

%% make figure 6

figure; set(gcf,'position',[50 50 900 350],'paperpositionmode','auto','Color','w');
cd /data/SO6/TPOSE/projects/pco2_correlations
load data_fig6b

clf;
[c h]=contourf(X2(37:401),Y2(14:149),DC2(37:401,14:149,:)',-.1:.01:1); cb6=colorbar;
set(h,'linestyle','none');
caxis([0 1])
set(gca,'color',[.8 .8 .8]);
colormap('parula')
[testX1,testX2] = meshgrid(X2, Y2);
hold on; 
mjr=DCL(907,265,1); % for 18 N, 110 W
mnr=DCL(907,265,2); % for 18 N, 110 W
thta=DCL(907,265,3); % for 18 N, 110 W
F2=exp(-((testX1.*cos(thta)+testX2.*sin(thta)).^2/mjr.^2+(testX1.*sin(thta)+testX2.*cos(thta)).^2./mnr^2));
v=[1/exp(1),1/exp(1)];
contour(X2,Y2,F2,v,'black','LineWidth',1); %show the contour line w/ s-major length, s-minor length
xlabel('zonal distance from center (km)');
ylabel('meridional distance from center (km)');
set(gca,'fontsize',16)
xlim([-3000 3000])
ylim([-1250 1250])
cb6.Title.String = 'r';

cd /home/winnie.chu/figures
eval(['export_fig Fig6_v3.png']);

%%

figure; set(gcf,'position',[50 50 1200 500],'paperpositionmode','auto','Color','w');
t=tiledlayout(2,2);
nexttile
imagesc(1,1,1)
nexttile
imagesc(2,2,2)
nexttile
imagesc(3,3,3)
nexttile
imagesc(4,4,4); cbtest = colorbar;
cbtest.Layout.Tile = 'east';
%t.TileSpacing = 'tight';

%%
major=DCL(:,:,1);
minor=DCL(:,:,2);
r2=DCL(:,:,4);

[row, col] = find(imag(major)~=0); % finding the weird points
%length(row)/(size(major,1)*size(major,2)) % proportion of points that are weird


%% make figure 7

load /home/winnie.chu/data/pCO2_100floats_all_v3

figure; set(gcf,'position',[50 50 1200 500],'paperpositionmode','auto','Color','w');
tiles=tiledlayout(2,2);

nexttile
n=25;
%axes;
%set(gca,'position',[.1 .55 .4 .36]);
[c h]=contourf(x(Irange),y(Jrange),PlotMapErr(:,:,n)',0:.001:1);
t=title('\bf a) n=25:\rm 38.32% coverage');
ax = gca;
ax.TitleHorizontalAlignment = 'left';
%t.FontWeight='bold';
set(h,'linestyle','none');
caxis([0 1])
ylim([-20 20])
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wo','MarkerSize',12,'LineWidth',2.5)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wx','MarkerSize',12,'LineWidth',2.5)
set(gca,'xtick',[],'ytick',-20:10:20,'Color', [.8 .8 .8],'fontsize', 18);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'10$^\circ$S','0$^\circ$','10$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$S'});
%set(gca, 'Color', [.8 .8 .8],'fontsize', 16)
%text(130,22.4,'a) n=25','fontsize',14,'fontweight','bold')
t.FontSize = 18;

nexttile
n=50;
%axes;
%set(gca,'position',[.55 .55 .4 .36]);
[c h]=contourf(x(Irange),y(Jrange),PlotMapErr(:,:,n)',0:.001:1);
t=title('\bf b) n=50:\rm 53.83% coverage');
%t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
%t.FontWeight='bold';
set(h,'linestyle','none');
caxis([0 1])
ylim([-20 20])
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wo','MarkerSize',12,'LineWidth',2.5)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wx','MarkerSize',12,'LineWidth',2.5)
set(gca,'xtick',[],'ytick',[],'Color', [.8 .8 .8],'fontsize', 16);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%text(130,22.4,'b) n=50','fontsize',14,'fontweight','bold')
t.FontSize = 18;

nexttile
n=75;
%axes;
%set(gca,'position',[.1 .1 .4 .36]);
[c h]=contourf(x(Irange),y(Jrange),PlotMapErr(:,:,n)',0:.001:1);
t=title('\bf c) n=75:\rm 61.35% coverage');
%t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
caxis([0 1])
ylim([-20 20])
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wo','MarkerSize',12,'LineWidth',2.5)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wx','MarkerSize',12,'LineWidth',2.5)
set(gca,'xtick',[150:30:290],'ytick',-20:10:20,'Color', [.8 .8 .8],'fontsize', 16);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'150$^\circ$E','180$^\circ$','150$^\circ$W','120$^\circ$W','90$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'10$^\circ$S','0$^\circ$','10$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$S'});
%text(130,22.4,'c) n=75','fontsize',14,'fontweight','bold')
t.FontSize = 18;

nexttile
n=100;
%axes;
%set(gca,'position',[.55 .1 .4 .36]);
[c h]=contourf(x(Irange),y(Jrange),PlotMapErr(:,:,n)',0:.001:1); cb7=colorbar;
t=title('\bf d) n=100:\rm 66.60% coverage');
%t.FontSize = 16;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
t.FontWeight='bold';
set(h,'linestyle','none');
caxis([0 1])
ylim([-20 20])
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wo','MarkerSize',12,'LineWidth',2.5)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wx','MarkerSize',12,'LineWidth',2.5)
set(gca,'ytick',[],'xtick',[150:30:290],'Color', [.8 .8 .8],'fontsize', 16);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'150$^\circ$E','180$^\circ$','150$^\circ$W','120$^\circ$W','90$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%text(130,22.4,'d) n=100','fontsize',14,'fontweight','bold')
t.FontSize = 18;

% shared colorbar
cb7.Layout.Tile = 'east';
cb7.Label.String = 'Normalized mapping error';
cb7.FontSize=16;
tiles.TileSpacing = 'tight';

eval(['export_fig Fig7_v4.png']); 

%% make figure 8

load /home/winnie.chu/data/pCO2_100floats_all_v3 % if not already loaded

PlotCoverage0 = [0 PlotCoverage];

figure; set(gcf,'position',[50 50 1200 400],'paperpositionmode','auto','Color','w');

scatter(0:100,PlotCoverage0,200,[75/256 84/256 202/256],'.','LineWidth',3);
set(gca,'fontsize', 18);
xlabel('Number of floats deployed')
ylabel('Percent coverage')
cd /home/winnie.chu/figures
eval(['export_fig Fig8.png']); 
%% make figure 9 

figure; set(gcf,'position',[50 50 1100 700],'paperpositionmode','auto','Color','w');
tiles=tiledlayout(2,1);

% track P4
load /home/winnie.chu/data/pCO2_20float10_130_10_290.mat

nexttile;
n=10;
%axes;
%set(gca,'position',[.1 .55 .4 .36]);
[c h]=contourf(x(Irange),y(Jrange),PlotMapErr(:,:,n)',0:.001:1);
t=title('\bf a) n=10:\rm 13.99% coverage');
%t.FontSize = 18;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(h,'linestyle','none');
set(gca, 'Color', [.8 .8 .8],'FontSize', 22,'xtick',[],'ytick',-20:10:20)
caxis([0 1])
ylim([-20 20])
hold on; plot([130 290],[y(float_locations(1,2)) y(float_locations(1,2))],'k--','LineWidth',2)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wo','MarkerSize',15,'LineWidth',2.5)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wx','MarkerSize',15,'LineWidth',2.5)
hold off;
%set(gca,'xtick',[],'ytick',-10:10:10);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$S'});
%text(130,22.4,'a) ','fontsize',14,'fontweight','bold')
%text(216,22.4,'n=20: 15.3% coverage','fontsize',14,'fontweight','normal')
t.FontSize = 18;

% track P16
load /home/winnie.chu/data/pCO2_20float-20_207_20_207.mat

nexttile;
n=10;
%axes;
%set(gca,'position',[.1 .1 .4 .36]);
[c h]=contourf(x(Irange),y(Jrange),PlotMapErr(:,:,n)',0:.001:1); cb=colorbar;
t=title('\bf b) n=10:\rm 16.71% coverage');
%t.FontSize = 18;
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(h,'linestyle','none');
caxis([0 1])
ylim([-20 20])
hold on; plot([x(float_locations(1)) x(float_locations(1))],[-20 20],'k--','LineWidth',2)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wo','MarkerSize',15,'LineWidth',2.5)
hold on; plot(x(float_locations(1:n,1)),y(float_locations(1:n,2)),'wx','MarkerSize',15,'LineWidth',2.5)
hold off;
set(gca, 'Color', [.8 .8 .8],'FontSize', 22,'xtick',150:30:290,'ytick',-20:10:20);
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'150$^\circ$E','180$^\circ$','150$^\circ$W','120$^\circ$W','90$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
%yticklabels({'10$^\circ$S','0$^\circ$','10$^\circ$N'});
yticklabels({'20$^\circ$S','10$^\circ$S','0$^\circ$','10$^\circ$N','20$^\circ$S'});
%text(130,22.4,'b) ','fontsize',14,'fontweight','bold')
%text(216,22.4,'n=20: 17.9% coverage','fontsize',14,'fontweight','normal')
t.FontSize = 18;

% shared colorbar
cb.Layout.Tile = 'east';
cb.Label.String='Normalized mapping error';
tiles.TileSpacing = 'tight';

eval(['export_fig Fig9_v5.png']); 





