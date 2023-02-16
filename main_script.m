%When you first enter Putty
screen %to save your work
matlab -softwareopengl -nodesktop

%Adding paths + loading stuff for plotting/saving figures
addpath /home/mmazloff/ANALYSIS  %this has imagescnan and rdmds
addpath /home/mmazloff/ANALYSIS/export_fig/ %this has export fig
cd /data/SO6/TPOSE/bgc_tpose6_k51/2004_fwd_run/diags 
%cd /data/SO6/TPOSE/bgc_tpose6/2004_fwd_run_SAVE/diags/ %new model doesn’t work
load  /data/SO6/TPOSE/bgc_tpose6_k51/grid/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC
cd /home/winnie.chu/data
%load velocityfields.mat U V Usurf Vsurf
load pCO2update.mat pCO2res
load DCall_1deg_update.mat DCL

%Loading velocity field data
U = rdmds('diag_state',1460,'rec',3); % east-west (zonal)?
V = rdmds('diag_state',1460,'rec',4); % north-south (meridional)?
W = rdmds('diag_state',1460,'rec',5); % vertical velocity? 

%Loading pCO2 data
pCO2 = zeros(1128,336,189,'single'); % 189 will change depending on the runs done
pCO2= rdmds(['diag_surf'],nan,'rec',2);
pCO2(pCO2==0)=nan;
pCO2 = pCO2(:,:,13:180); % cropping first year (but make sure the number of total months is divisible by 12)
nmonths = size(pCO2,3);
% should be saved under cd /home/winnie.chu/data
% called pCO2.mat
% load pCO2.mat should be good and contains pCO2 and pCO2res
% save pCO2update.mat pCO2 pCO2res nmonths

%Loading surface O2 data (only works on old model, NEED UPDATE)
O2 = zeros(1128,336,51,'single'); % third dimension is depth, with 51 profiles
for i = 1:7
O2(:,:,:,(i-1)*12+1:i*12) = rdmds(['run_' num2str(i) '/diag_bgc'],nan,'rec',3);
end
O2(O2==0)=nan;

%Making movie of 7 yr pCO2 data
cd /home/winnie.chu/figures % goes to the right directory to save file
v = VideoWriter('pCO2_NEW.avi','Uncompressed AVI');
v.FrameRate = 3; %set to 3 frames per second
open(v)
lower=min(pCO2,[],'all');
upper=max(pCO2,[],'all');
for t=1:nmonths
imagescnan(x, y, pCO2(:,:,t)'); axis xy; a=colorbar;
set(gcf, 'Color', 'w')
caxis([lower,upper]) %so that everything shows up in the same scale
a.Label.String = 'pressure (atm)';
title('pCO2 over 7 years');
xlabel('longitude');
ylabel('latitude');
drawnow
frame= getframe(gcf);
writeVideo(v,frame);
close;
end
close(v);


%Making movie of 7 yr surface O2 data (NEEDS UPDATE)
for t=1:84
imagescnan(x, y, O2(:,:,1,t)');axis xy; a=colorbar;
set(gcf, 'Color', 'w')
caxis([0.15 0.32]) %so that everything shows up in the same scale
a.Label.String = 'concentration (mol/m^3)';
title('surface O2 over 7 years');
xlabel('longitude');
ylabel('latitude');
drawnow
M(t)= getframe(gcf);
close;
end
movie(M) % to play the movie
cd /home/winnie.chu/figures % goes to the right directory to save file
movie2avi(M,'surface_O2_7yrs.avi', 'fps', 2); % to save the movie

%Saving movie/figures to local computer (do in command window/gitbash not in Putty)
scp winnie.chu@aha.ucsd.edu:/home/winnie.chu/data/DCall_1deg_update.mat /c/Users/16262/desktop/SURF

% Pushing files to the remote server
scp /c/Users/16262/desktop/SURF/plot_float_coverage.m winnie.chu@aha.ucsd.edu:/home/winnie.chu/functions



%Convert the avi file to gif using this website!
https://www.onlineconverter.com/avi-to-gif 

%Plotting time mean of pCO2 data
tmean_pCO2=nanmean(pCO2,3);
imagescnan(x,y,tmean_pCO2'); axis xy; a=colorbar; %automatically in 2d
% sometimes does weird things, would run this in two sections, figure then labels
a.Label.String = 'pressure (atm)';
title(sprintf('pCO2 time mean over %d months', nmonths))
xlabel('longitude');
ylabel('latitude');
set(gcf, 'Color', 'w')
eval(['export_fig  pCO2_time_mean.png']); %export as any file type, make sure correct directory 

%Plotting time STD of pCO2 data 
time_std_pCO2=nanstd(pCO2,0,3);
% size(time_std_pCO2) %checks that mean is 2d array
imagescnan(x, y, time_std_pCO2');axis xy; a=colorbar; % creates 2d colormap 
a.Label.String = 'pressure (atm)';
title(sprintf('pCO2 time standard deviation over %d months', nmonths))
xlabel('longitude');
ylabel('latitude');
set(gcf, 'Color', 'w')
eval(['export_fig  pCO2_time_STD.png']);

%Plotting spatial mean of 7 yr pCO2 data
for t=1:nmonths
tmp = pCO2(:,:,t);
pCO2mean(t) = nanmean(tmp(:));
end
plot(1:nmonths, pCO2mean(1:nmonths))
title(sprintf('pCO2 spatial mean over %d months', nmonths))
xlabel('months')
ylabel('pressure (atm)')
set(gcf, 'Color', 'w')
eval(['export_fig  pCO2_spatial_mean.png']);

%Plotting spatial STD of 7 yr pCO2 data
for t=1:nmonths
tmp = pCO2(:,:,t);
pCO2std(t) = nanstd(tmp(:));
end
plot(1:nmonths, pCO2std(1:nmonths))
title(sprintf('pCO2 spatial STD over %d months', nmonths))
xlabel('months')
ylabel('pressure (atm)')
set(gcf, 'Color', 'w')
eval(['export_fig  pCO2_spatial_STD.png']);

%Isolating signal without annual cycles
% initializing paths and arrays
addpath /home/mmazloff/ANALYSIS/
addpath /home/mmazloff/ANALYSIS/OSSE
addpath /home/mmazloff/ANALYSIS/SOSE
addpath /home/mmazloff/ANALYSIS/export_fig
addpath /data/mmazloff/ANALYSIS/m_map
load /home/mmazloff/ANALYSIS/cmap
fpath =  '/data/SO6/TPOSE/bgc_tpose6_k51/2004_fwd_run/diags/';
cd /data/SO6/TPOSE/bgc_tpose6_k51/2004_fwd_run/diags 
load  /data/SO6/TPOSE/bgc_tpose6_k51/grid/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC
NT = nmonths; 
NX = 1128;
NY = 336;
Time = [1:365.25/12:365.25*NT/12];Time = Time' - mean(Time); 
pCO2res = zeros(1128,336,nmonths,'single');
pCO2VE = zeros(1128,336,5,'single');

% simple pCO2 calculations + wavenumber stuff
pCO2_bar = mean(pCO2,3);
pCO2_std = std(pCO2,[],3);
pCO2cos = zeros(NX,NY,4,'single'); pCO2sin = pCO2cos;
pvar = 1;dvar= .01;ifwd = 1;
kx1 = (2*pi)*[1]/[365.25];%wavenumbers to fit to: 1 yr fit to pCO2
kx2 = (2*pi)*[2]/[365.25];%wavenumbers to fit to: 6 mnt fit to pCO2
kx3 = (2*pi)*[3]/[365.25];%wavenumbers to fit to: 4 mnt fit to pCO2
kx4 = (2*pi)*[4]/[365.25];%wavenumbers to fit to: 3 mnt fit to pCO2
TimeT = Time;
clear Time %to be safe

% loop for variance and residual signal
for j = 1:NY
    NY-j
    for i = 1:NX
      tmp2 = squeeze(pCO2(i,j,:));tmp2a = detrend(tmp2 - mean(tmp2)); % remember to crop first year
      [p2,y2] = gft(TimeT,tmp2a,[kx1 kx2 kx3 kx4]',pvar,dvar,ifwd); % gft in home/ANALYSIS-- least squares fit of 4 harmonics to time series
      pCO2cos(i,j,:) = p2(1:4);%  cosin amp of kx1
      pCO2sin(i,j,:) = p2(5:8);%    sin amp of kx1
      y11 = p2(1)*cos(TimeT*kx1)+p2(5)*sin(TimeT*kx1); %annual cycle
      y12 = p2(2)*cos(TimeT*kx2)+p2(6)*sin(TimeT*kx2); %semi annual
      y13 = p2(3)*cos(TimeT*kx3)+p2(7)*sin(TimeT*kx3); % 4 month
      y14 = p2(4)*cos(TimeT*kx4)+p2(8)*sin(TimeT*kx4); % 3 month
      pCO2VE(i,j,1) = 100 - 100*var(tmp2-y11)/var(tmp2); % variance explained by annual cycle
      pCO2VE(i,j,2) = 100 - 100*var(tmp2-y12)/var(tmp2); % variance explained by semi annual
      pCO2VE(i,j,3) = 100 - 100*var(tmp2-y13)/var(tmp2); % variance explained by 4 month
      pCO2VE(i,j,4) = 100 - 100*var(tmp2-y14)/var(tmp2); % variance explained by 3 month
      pCO2VE(i,j,5) = 100 - 100*var(y2)/var(tmp2); %res var, how much variance is leftover from original signal
      pCO2res(i,j,:) = tmp2 - y2; % time series residual, outside of trends/cycles
      pCO2H_std(i,j) = std(y2);
      pCO2R_std(i,j) = std(tmp2 - y2);
   end
end
clear tmp*;

%Finding and plotting total variance in pCO2 (pre-processing)
pCO2var = var(pCO2,0,3);
imagescnan(x,y(34:277),pCO2var(:,34:277)','NanColor',[.8 .8 .8]);axis xy; a=colorbar; 
caxis([1e-10 10e-10])
a.Label.String = 'Total variance [atm^2]';
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$E','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$N','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
a.FontSize=14;
xlabel('Longitude');
ylabel('Latitude');
set(gca,'fontsize',16)
set(gcf,'Position',[50 50 1100 500],'Color', 'w');
cd /home/winnie.chu/figures
eval(['export_fig  pCO2totalvar_v4.png']);


%Plotting variance explained by annual cycle
imagescnan(x,y(34:277),pCO2VE(:,34:277,1)','NanColor',[.8 .8 .8]);axis xy;
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$E','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
yticklabels({'20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$N','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N'});
xlabel('Longitude');
ylabel('Latitude');
set(gca,'fontsize',16)
set(gcf,'Position',[50 50 1100 500],'Color', 'w');
col=[[256,256,256]'/256,[75,84,202]'/256]'; % two points - first is white, second is blue but can be ANY two colors
col=interp1([0:1],col,[0:1/999:1]); % interps on to 100 point grid…or 1000 or whatever you want, 100 is usually more than enough! 
colormap(col)
a=colorbar;
a.Label.String = 'Percent variance explained [%]';
a.FontSize=16;
cd /home/winnie.chu/figures
eval(['export_fig  pCO2VE_annual_v4.png']);

%Plotting residual variance 
imagescnan(x,y,pCO2VE(:,:,5)','NanColor',[.8 .8 .8]);axis xy; % creates 2d colormap 
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; xticklabels({'120$^\circ$E','140$^\circ$E','160$^\circ$E','180$^\circ$E','160$^\circ$W','140$^\circ$W','120$^\circ$W','100$^\circ$W','80$^\circ$W'});
yaxisproperties= get(gca, 'YAxis');
yaxisproperties.TickLabelInterpreter = 'latex';  
yticklabels({'25$^\circ$S','20$^\circ$S','15$^\circ$S','10$^\circ$S','5$^\circ$S','0$^\circ$N','5$^\circ$N','10$^\circ$N','15$^\circ$N','20$^\circ$N','25$^\circ$N'});
col=[[256,256,256]'/256,[75,84,202]'/256]'; % two points - first is white, second is blue but can be ANY two colors
col=interp1([0:1],col,[0:1/999:1]); % interps on to 100 point grid…or 1000 or whatever you want, 100 is usually more than enough! 
colormap(col);
a=colorbar;
a.Label.String = 'Percent residual variance [%]';
a.FontSize=16;
xlabel('Longitude');
ylabel('Latitude');
set(gca,'fontsize',16)
set(gcf,'Position',[50 50 1100 500],'Color', 'w');
cd /home/winnie.chu/figures
eval(['export_fig   pCO2_residualvariance_v5.png']);

% Getting phase of cycles
cycle_phases_y11 = zeros(1128,336,'single'); 
cycle_phases_y12=cycle_phases_y11; cycle_phases_y13=cycle_phases_y11; cycle_phases_y14=cycle_phases_y11;

for j = 1:NY
    NY-j
    for i = 1:NX
      tmp2 = squeeze(pCO2(i,j,:));tmp2a = detrend(tmp2 - mean(tmp2)); % remember to crop first year
      [p2,y2] = gft(TimeT,tmp2a,[kx1 kx2 kx3 kx4]',pvar,dvar,ifwd); % gft in home/ANALYSIS-- least squares fit of 4 harmonics to time series
      pCO2cos(i,j,:) = p2(1:4);%  cosin amp of kx1
      pCO2sin(i,j,:) = p2(5:8);%    sin amp of kx1
      y11 = p2(1)*cos(TimeT*kx1)+p2(5)*sin(TimeT*kx1); %annual cycle
      y12 = p2(2)*cos(TimeT*kx2)+p2(6)*sin(TimeT*kx2); %semi annual
      y13 = p2(3)*cos(TimeT*kx3)+p2(7)*sin(TimeT*kx3); % 4 month
      y14 = p2(4)*cos(TimeT*kx4)+p2(8)*sin(TimeT*kx4); % 3 month
	cycle_phases_y11(i,j)= atan2(p2(5),p2(1));
	cycle_phases_y12(i,j)= atan2(p2(6),p2(2));
	cycle_phases_y13(i,j) = atan2(p2(7),p2(3));
cycle_phases_y14(i,j) = atan2(p2(8),p2(4));
   end
end
clear tmp*;

%Saving phase plot for annual cycle
%test_cycle_phases_y11= cycle_phases_y11*(12/(2*pi))+6; % minimum pCO2
test_cycle_phases_y11= cycle_phases_y11*(12/(2*pi)); % maximum pCO2
test_cycle_phases_y11(test_cycle_phases_y11<=0)=test_cycle_phases_y11(test_cycle_phases_y11<=0)+12; % maximum pCO2
imagescnan(x,y,test_cycle_phases_y11');axis xy; a=colorbar('Ticks',[0.01, 1,2,3,4,5,6,7,8,9,10,11],...
         'TickLabels', {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'});
a.Label.String = 'month of peak pCO2';
%title('Month of maximum pCO2 values')
set(gcf, 'Color', 'w')
xlabel('longitude');
ylabel('latitude');
colormap(hsv); %gets a wrapped colorbar
cd /home/winnie.chu/figures
eval(['export_fig  pCO2phase_annual_v5.png']);

%Visualizing the correlation map
cd /home/winnie.chu/data
xk = x*deg2km(1);yk = y*deg2km(1); 
NX = 1128; NY = 336; NT = 180;
amp=1;
dstx = 240; %max grid points to look in x
dsty = 72; %max grid points to look in y matching aspect ratio of grid 1128x336
Q=pCO2res; %NOTE USING PCO2 RESIDUAL SIGNAL
QQ = Q(:,:,1:6)*nan;
NoiseFac = 10;
SNR = 1/NoiseFac^2;
NoiseFac1=10;NoiseFac2 = 1;
NSEC = 18; %ONLY GO OUT TO A THRESHOLD. DO IN 360/NSEC deg SECTORS
thrshcntr = .8;
tic 
i=906; %336 for 160 E
j=264; %156 for 0 N
loop=0;
mjr=0; %throwaway
mnr=0; %throwaway
errx=550;
erry=550;
        while errx>0.01 && erry>0.01 && loop<6
            prev_mjr=mjr;
            prev_mnr=mnr;
            if loop==0
                    LXU=200; LYU=200;
            else
                    LXU=20; LYU=20;
            end
            A_p = [0.5*(LXU^-2) 0 1.5*(LYU^-2)]';
            P_p = [A_p(1)^2 A_p(1)*A_p(3) A_p(3)^2]';
            Pinv = diag(P_p.^-1); % how much we trust the prior, weird bc this means the larger the uncertainty the more we trust
            I = i-dstx:i+dstx; I(I>NX)=[]; I(I<1)=[]; 
            J = j-dsty:j+dsty; J(J>NY)=[]; J(J<1)=[];
            snx = length(I);
            sny = length(J);
            DC = nan(snx*sny,1,'single'); X = DC; Y = DC; % decorrelation field initialization
        %MATT NOTES: WANT TO ALSO HAVE ON 2D grid. Make 2D gid
            %X2 = linspace(-4200,4200,snx/1.1); %used to be 1.5 instead of 1.1 for both x2 and y2
            %Y2 = linspace(-3500,3500,sny/1.1);
            X2 = linspace(-4200,4200,snx/1.1); %1.1
            Y2 = linspace(-1500,1500,sny/0.89); %0.89 works!!!
            DC2 = nan(length(X2),length(Y2),'single');
            tmp1 = squeeze(Q(i,j,:));
            count = 0;
            if isfinite(tmp1(1))
                for ii = I
                    for jj = J
                        tmp2 = squeeze(Q(ii,jj,:));
                        if isfinite(tmp2(1))
                            a = (xk(ii) - xk(i));
                            dx = a*cosd((y(jj) + y(j))*.5);% scaled distance in X so 1 degree is equal everywhere despite latitude
                            dy =  yk(jj) - yk(j);
                            if (dx.^2 + dy.^2) < 16000000 %max distance
                                count = count + 1;                
                                X(count) = dx; 
                                DC(count) = corr(tmp1,tmp2); 
                                Y(count) = dy;
        %MATT NOTES 2: PUT ON 2D GRID
             iii = min(find(dx<X2));
             jjj = min(find(dy<Y2));
             DC2(iii,jjj) = DC(count);
                            end
                        end
                    end
                end
                NP = count;
                count = 0;K=0;
                [p,r] = cart2pol(X(1:NP),Y(1:NP)); p = p + pi; %MAP TO POLAR
                p = ceil(p*NSEC/2/pi); %SPLIT INTO DISCRETE SECTORS, ceil means round up
                [trash,I] = sort(r); %SORT SO LOOK OUT IN RADIUS, looks for the radius’ indices
                maxr = ones(1,NSEC)*999999;%INIT OK RADIUS 

                for ti = I'
                    if r(ti) < maxr(p(ti))
                        if DC(ti)>thrshcntr %THIS IS RADIUS WILL USE FOR PRIOR
                            count = count + 1; 
                            K(count) = ti;
                        else
                            maxr(p(ti)) = r(ti);%WE ARE BAD. THIS IS NEW MAX RADIUS FOR SECTOR
                        end
                    end 
                end
        %  CHOOSE PRIOR BASED ON maxr
                %if 1 %CHANGE PRIOR BASED ON maxr
                if loop==0
                    LXP = 1/sqrt(-2*log(thrshcntr))*mean(maxr([1 2 8:11 17:18])); %zonal
                    LYP = 1/sqrt(-2*log(thrshcntr))*mean(maxr([3:7  12:16])); %meridional
                else 
                    LXP = prev_mjr; %zonal
                    LYP = prev_mnr; %meridional (how come we always set y as minor?)
                end
                A_p = [0.5*(LXP^-2) 0 0.5*(LYP^-2)]'; %reassigns prior to the same as before
                %end
                %NOW FIT TO GAUSSIAN
                K = find((DC>0)==1);%THIS IS OK AS LONG AS NoiseFac2 > 0, points where the correlation is greater than 0
                if length(K)>24
                    ym = -log(DC(K));                   
                    H = [X(K).^2, X(K).*Y(K), Y(K).^2];
                    Rinv  = (NoiseFac1*( 1.01-exp(NoiseFac2*(-A_p(1)*X(K).^2 - A_p(2)*X(K).*Y(K) - A_p(3)*Y(K).^2) ) )).^-1 ; % seems like original code has noisefac2 only affecting the LXP guess, not theta or LYP, and 1.01=1+SNR
                %Rinv=10 implementing the no prior guess
                    R2inv = repmat(Rinv.^2,1,3); % how much we trust data
                  %least squares fit  
                    Hty = (H.*R2inv)'*(ym - H*A_p);
                    A  = A_p + ((H.*R2inv)'*H + Pinv)\Hty; % A_p + R2inv/(R2inv+Pinv)*(ym-A_P) in scalar form 
                    R2  = 1 - var(Rinv.*(exp(-ym) - exp(-H*A)))/var(Rinv.*exp(-ym));
                    [V,L] = eig([A(1), A(2)/2; A(2)/2, A(3)]);
                    mjr = (2*L(1,1))^-.5;
                    mnr = (2*L(2,2))^-.5;
                    thta = atan2(V(1,1),V(2,1)) - pi/2;
                    if thta>pi/2 
                        thta = thta-pi;%rotate 180 so goes from -pi/2 to pi/2
                    elseif thta<-pi/2 
                        thta = thta+pi;
                    end %thta*360/2/pi
                    QQ(i,j,:) = [mjr mnr thta R2 LXP LYP];
                end
            end
            errx=abs(mjr-prev_mjr)/mjr;
            erry=abs(mnr-prev_mnr)/mnr;
		loop = loop + 1
end



imagescnan(X2(37:401),Y2(14:149),DC2(37:401,14:149,:)','NanColor',[.8 .8 .8]);axis xy; a=colorbar;
a.Label.String = 'Correlation (r)';
a.FontSize=16;
[testX1,testX2] = meshgrid(X2, Y2);
hold on; 
F2=exp(-((testX1.*cos(thta)+testX2.*sin(thta)).^2/mjr.^2+(testX1.*sin(thta)+testX2.*cos(thta)).^2./mnr^2));
v=[1/exp(1),1/exp(1)];
contour(X2,Y2,F2,v,'black','LineWidth',1); %show the contour line w/ s-major length, s-minor length
xlabel('zonal distance from center (km)');
ylabel('meridional distance from center (km)');
set(gca,'fontsize',16)
set(gcf,'Position',[50 50 1000 400],'Color', 'w');
cd /home/winnie.chu/figures
eval(['export_fig   gfit_160E0N_v2.png']);


%Calculating convergent gaussian fit (still being modified, last working version)
cd /home/winnie.chu/data
%fileID = fopen('QQ_array','a+');

% from the other code Annual CyclePaper_GausFit_R1_V3
xk = x*deg2km(1);yk = y*deg2km(1); 
NX = 1128; NY = 336; NT = 72;
amp=1;
dstx = 240; %max grid points to look in x
dsty = 72; %max grid points to look in y matching aspect ratio of grid 1128x336
Q=pCO2res; %NOTE USING PCO2 RESIDUAL SIGNAL
QQ = Q(:,:,1:6)*nan;
NoiseFac = 10;
SNR = 1/NoiseFac^2;
NoiseFac1=10;NoiseFac2 = 1;
NSEC = 18; %ONLY GO OUT TO A THRESHOLD. DO IN 360/NSEC deg SECTORS
thrshcntr = .8;
mjr=0; %throwaway
mnr=0; %throwaway
tic 
for i=276:12:996 % x, where spacing matches aspect ratio of x to y grid
    % between 150 to 270 degrees E
    i
for j=90:12:210 % y
    % between -10 to 10 degrees N
    j
	loop=0;
    mjr=0; %throwaway
    mnr=0; %throwaway
    errx=550;
    erry=550;
while errx>0.01 && erry>0.01 && loop<6
    prev_mjr=mjr;
	prev_mnr=mnr;
    if loop==0
			LXU=200; LYU=200;
    else
			LXU=20; LYU=20;
    end
    A_p = [0.5*(LXU^-2) 0 1.5*(LYU^-2)]';
    P_p = [A_p(1)^2 A_p(1)*A_p(3) A_p(3)^2]';
    Pinv = diag(P_p.^-1); % how much we trust the prior, weird bc this means the larger the uncertainty the more we trust
      I = i-dstx:i+dstx; I(I>NX)=[]; I(I<1)=[]; 
      J = j-dsty:j+dsty; J(J>NY)=[]; J(J<1)=[];
      snx = length(I); sny = length(J);
      DC = nan(snx*sny,1,'single'); X = DC; Y = DC; % decorrelation field initialization
%MATT NOTES: WANT TO ALSO HAVE ON 2D grid. Make 2D gid
      %X2 = linspace(-4200,4200,snx/1.5);
      %Y2 = linspace(-3500,3500,sny/1.5);
	  %DC2 = nan(length(X2),length(Y2),'single');
      tmp1 = squeeze(Q(i,j,:));
      count = 0;
      if isfinite(tmp1(1))
        for ii = I
          for jj = J
            tmp2 = squeeze(Q(ii,jj,:));
            if isfinite(tmp2(1))
              a = (xk(ii) - xk(i));
              dx = a*cosd((y(jj) + y(j))*.5);% scaled distance in X so 1 degree is equal everywhere despite latitude
              dy =  yk(jj) - yk(j);
              if (dx.^2 + dy.^2) < 16000000 %max distance
                count = count + 1;                
                X(count) = dx; 
                DC(count) = corr(tmp1,tmp2); 
                Y(count) = dy;
%MATT NOTES 2: PUT ON 2D GRID
           %iii = min(find(dx<X2));
     %jjj = min(find(dy<Y2));
     %DC2(iii,jjj) = DC(count);
              %Y(count) = dy;
              end
            end
          end
        end
        NP = count;
        count = 0;K=0;
        [p,r] = cart2pol(X(1:NP),Y(1:NP)); p = p + pi; %MAP TO POLAR
        p = ceil(p*NSEC/2/pi); %SPLIT INTO DISCRETE SECTORS, ceil means round up
        [trash,I] = sort(r); %SORT SO LOOK OUT IN RADIUS, looks for the radius’ indices
        maxr = ones(1,NSEC)*999999;%INIT OK RADIUS 

        for ti = I'
          if r(ti) < maxr(p(ti))
            if DC(ti)>thrshcntr %THIS IS RADIUS WILL USE FOR PRIOR
              count = count + 1; 
              K(count) = ti;
            else
              maxr(p(ti)) = r(ti);%WE ARE BAD. THIS IS NEW MAX RADIUS FOR SECTOR
            end
          end 
        end
%  CHOOSE PRIOR BASED ON maxr
        %if 1 %CHANGE PRIOR BASED ON maxr
            if loop==0
                LXP = 1/sqrt(-2*log(thrshcntr))*mean(maxr([1 2 8:11 17:18])); %zonal
                LYP = 1/sqrt(-2*log(thrshcntr))*mean(maxr([3:7  12:16])); %meridional
            else 
                LXP = prev_mjr; %zonal
                LYP = prev_mnr; %meridional (how come we always set y as minor?)
            end
           A_p = [0.5*(LXP^-2) 0 0.5*(LYP^-2)]'; %reassigns prior to the same as before
        %end
        %NOW FIT TO GAUSSIAN
        K = find((DC>0)==1);%THIS IS OK AS LONG AS NoiseFac2 > 0, points where the correlation is greater than 0
        if length(K)>24
          ym = -log(DC(K));                   
          H = [X(K).^2, X(K).*Y(K), Y(K).^2];
          Rinv  = (NoiseFac1*( 1.01-exp(NoiseFac2*(-A_p(1)*X(K).^2 - A_p(2)*X(K).*Y(K) - A_p(3)*Y(K).^2) ) )).^-1 ; % seems like original code has noisefac2 only affecting the LXP guess, not theta or LYP, and 1.01=1+SNR
	    %Rinv=10 implementing the no prior guess
          R2inv = repmat(Rinv.^2,1,3); % how much we trust data
          %least squares fit  
          Hty = (H.*R2inv)'*(ym - H*A_p);
          A  = A_p + ((H.*R2inv)'*H + Pinv)\Hty; % A_p + R2inv/(R2inv+Pinv)*(ym-A_P) in scalar form 
          R2  = 1 - var(Rinv.*(exp(-ym) - exp(-H*A)))/var(Rinv.*exp(-ym));
          [V,L] = eig([A(1), A(2)/2; A(2)/2, A(3)]);
          mjr = (2*L(1,1))^-.5;
          mnr = (2*L(2,2))^-.5;
          thta = atan2(V(1,1),V(2,1)) - pi/2;
          if thta>pi/2 
             thta = thta-pi;%rotate 180 so goes from -pi/2 to pi/2
          elseif thta<-pi/2 
             thta = thta+pi;
          end %thta*360/2/pi
          QQ(i,j,:) = [mjr mnr thta R2 LXP LYP];
        end
      end
      errx=abs(mjr-prev_mjr)/mjr;
      erry=abs(mnr-prev_mnr)/mnr;
      loop=loop+1
    end
    %fprintf(fileID,'%g %g %g %g %g\n',i,j,mjr,mnr,thta);
  end
end

%fclose('all') % or fclose(fileID) or fclose(‘all’) if that doesnt work
t=toc;
save qq_fifth.mat QQ t

%Plotting semi-major lengths (replace 1 index with 2 or 3 or 4 for mnr or thta or r^2)
[row tmp]=find(~isnan(DCL(:,:,1))); % row is row, col is linear index
[col d]=ind2sub([336,1],tmp); % turn linear index into column and dimension were d is just a bunch of 1s 
QQ2=DCL; % just to make sure we don’t lose the data
QQ2(QQ2(:,:,1)<10)=nan; % bc there’s a bunch of random min that are 1.681 which
imagescnan(x(row),y(col),real(QQ2(row,col,1)')); axis xy; a=colorbar;
caxis([1.3 5.2]*1e3) % find based on extrema
a.Label.String = 'decorrelation length (km)';
title('semimajor length scales')
xlabel('longitude');
ylabel('latitude');
set(gcf, 'Color', 'w')
set(gca,'fontsize', 16);
cd /home/winnie.chu/figures
eval(['export_fig  mjr_2x2deg.png']); %export as any file type

%Covariance plotting
%% testing the code
% a simple model of the covariance code
prompt = 'Input number of floats deployed: ';
num_iter = input(prompt);
% need to run the above portion first and enter input before subsequent part

% Initialization
QQ=pCO2res;
%QQ(isnan(QQ))=0;
%CC=zeros(size(QQ,1),size(QQ,2));
%float_QQ=zeros(size(QQ,1),size(QQ,2),size(QQ,3),num_iter);
%count=0;
max_cov=1e-9;
plot_var=zeros(size(QQ,1),size(QQ,2)); %initialization
oneQQ=ones(size(QQ,1),size(QQ,2)); % ones matrix the same size as QQ
float_locations=zeros(num_iter,2);
%oldmax_ij=[0 0];
xk = x*deg2km(1);yk = y*deg2km(1); 
indices = find(isnan(pCO2(:,:,1)));
plot_var(indices)=nan; %remove points on land
imin=276;
imax=996; %range of i points
igap=6;
jmin=96;
jmax=216; %range of j points
jgap=6;
CC=zeros(size(QQ,1),size(QQ,2));
float_I=[];
float_J=[];
dcx=3500; % just general estimate of decorrelation lengths in x direction
dcy=1000; % estimate of decorrelation lengths in y direction
% replace this later once you get all decorrelation lengths
% make sure that when you do correlation analysis it matches this one so
% that we can have a decorrelation length for each of the points  
PlotMapErr=zeros(length(imin:igap:imax),length(jmin:jgap:jmax),num_iter);
PlotCoverage=zeros(num_iter);

% Finding model-model, model-float, and float-float covariance
%MATT’S ATTEMPT AT Cmm
%SUBSAMPLE DATA
Qss = QQ(imin:igap:imax,jmin:jgap:jmax,:);
[NX,NY,NT] = size(Qss); % getting the dimensions of the compressed matrix
Qssv = reshape(Qss,NX*NY,NT); % turning the xy into one dimension vs. time
[IsOcean,junk] = find(isnan(Qssv(:,1))==0); % finding all the indices that are not land
Qsssv = Qssv(IsOcean,:); % compressing further to only ocean points
NM = length(IsOcean); % new length of xy that is on ocean
Cmm = Qsssv*(Qsssv'); % Cmm is the local variances size should be NMxNM

%TO VISUALIZE:
Qvar = zeros(NX*NY,1); Qvar(IsOcean) = diag(Cmm); Qvar = reshape(Qvar,NX,NY);
% Qvar is the local variances of each of the ocean points in the model
%NEXT IS TO FIND Cmd which is covariance at data points
%Qsssv is signal at every model point
%NOW HAVE DECIDED TO DO THIS IN DATA SPACE, SO YOU CAN DELETE CMM IF YOU WANT

%HERE GOES WINNIE’S LOOP TO FIND THE DATA LOCATIONS I AND J
% starting the actual iteration
tic
for iter=1:num_iter
    max_ij=[0,0]; %initializing location of first float
    count=0; %initializing count
    for i=imin:igap:imax
        for j=jmin:jgap:jmax
            tmp1=squeeze(QQ(i,j,:)); %getting the time series of first point
            if isnan(tmp1(1))==1
                %continue
                %count=0; 
            elseif count > 0
            %if count > 0
                % if not the first run, then find the total covariance and
                % compare to previous covariances
                %cov_ij=sum(abs(CC),'all','omitnan');
                cov_ij=sum(CC,'all','omitnan');
                %test_ij=isequal([i j], oldmax_ij); %is this an existing float location
                if cov_ij > max_cov %&& test_ij==0
                    % if not an existing float location and bigger
                    % covariance
                    max_cov=cov_ij; %store the highest covariances
                    all_covs=CC; %store covariances over entire area
                    max_ij=[i,j] %store float location
                    %CC=zeros(size(QQ,1),size(QQ,2));
                end
            end
            count=count+1;
            for I=imin:igap:imax
                for J=jmin:jgap:jmax
                    tmp2=squeeze(QQ(I,J,:)); %get the time series for 2nd point
                    if isnan(tmp2(1))==1
                        continue
                    end
                    %covmatrix=cov(tmp1,tmp2);
                    a = (xk(I) - xk(i));
                    dx = a*cosd((y(J) + y(j))*.5);% scaled distance in X 
                    dy =  yk(J) - yk(j);
                    % (admittedly large) radius of influence
                    %if (dx.^2 + dy.^2) < 1000000 %max distance
                    %if dx <= 1500 && dy <= 350 && (dx.^2 + dy.^2) < 1000000
                        % remember to half dx and dy if you're using
                        % correlation lengths
                        covmatrix=cov(tmp1,tmp2);
                        %CC(I,J)=round((1/(dist+1)).*covmatrix(2,1),3,'significant');
                        CC(I,J)=round(covmatrix(2,1)*exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))),3,'significant'); %covariance of point ij with point IJ
                    %else
                        %CC(I,J)=0;
                    %end
                end
            end
        end
    end
    float_locations(iter,:)=max_ij;
    float_I=[float_I round((max_ij(1)-imin)/igap)]; %getting the float location indices in compressed QQ
    float_J=[float_J round((max_ij(2)-jmin)/jgap)];
    ND=length(float_I);
    %NOW THIS GOES AFTER YOUR ITERATIVE LOOP. RUN EACH TIME TO GET “FORMAL MAP ERROR”
    %Need signal just at every data point
    MASK = zeros(NX,NY);
    DIST = zeros(NX,NY,ND);
    for num = 1:ND %NDATA
        MASK(float_I(num),float_J(num)) = 1;
    end
    for num = 1:ND
        I=round(float_I(num)*igap+imin);
        J=round(float_J(num)*jgap+jmin);
        for nx=1:NX
            for ny=1:NY
                i=round(nx*igap+imin);
                j=round(ny*jgap+jmin);
                a = (xk(I) - xk(i));
                dx = a*cosd((y(J) + y(j))*.5);% scaled distance in X
                dy = yk(J) - yk(j);
                DIST(nx,ny,num)=round(exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))),3,'significant'); %gaussian distances
            end
        end
    end
    %TO TEST
    %I = round(([390 492 948 324 846] - imin)/6);
    %J = round(([180 174 168 162 174] - jmin)/6);
    %ND = 5;
    %NOW PLOT VarUnexp
%     figure()
%     title(sprintf('percent coverage of %d floats is %.2f%%',iter,coverage_percent))
%     imagescnan(x(imin:igap:imax),y(jmin:jgap:jmax),MapErr'); axis xy; a=colorbar;
%     caxis([0,1])
%     hold on
%     text(x(float_locations(:,1)),y(float_locations(:,2)),'\otimes','Color','white','FontSize',16,'HorizontalAlignment','center','FontWeight','bold')
%     hold off
    if max_ij(1)==0
        % no max_ij means we found nothing to report
        continue
    else 
        DIST=reshape(DIST,NX*NY,ND); % turn distances into a linear vector for each float
        %DISTsss=DIST(IsOcean,:); %only needs floats in ocean should be lenght of IsOcean x ND
        MASK = MASK(:); % turned into a linear vector rather than matrix
        %MASKsss = MASK(IsOcean); %only need the points in ocean
        %[Is] = find(MASKsss==1); %where the floats are
        %MASKsss = MASK(IsOcean); %only need the points in ocean
        [Is] = find(MASK==1); %where the floats are
        D = Qssv(Is,:); % getting just the float points
        float_ind = sub2ind(size(Qss(:,:,1)),float_I,float_J);
        DDist =DIST(float_ind,:); %use float_ind instead of IS to make symmetric dist
        %Cmd = Qsssv*(D'); % finding the covariance between all those float points and the model
        Cmd = Qssv*(D').*DIST; % number of NXNY points that are ocean x number of floats
        %NOW GET DATA DATA COVARIANCE
        Cmd = Cmd(IsOcean,:);
        Cdd = D*(D').*DDist; % get covariance between float locations
        %ADD NOISE TO DATA BECAUSE FLOATS AREN’T PERFECT
        NoiseValue = 1e-5; %Noise variance in obs keep in mind order of magnitude of pCO2 values is around 10^-4
        %NoiseValue = 0;
        Cdd = Cdd + NoiseValue*eye(ND); % adding noise to the variance of float locations
        %FORMAL MAPPING ERROR
        %FME = diag(Cmm) - diag(Cmd*(Cdd)^-1*Cmd');
        FME = diag(Cmm) - diag(Cmd*(Cdd)^-1*(Cmd'));
        %NORMALIZED VERSION
        NME = FME./diag(Cmm);
        %%OLD   %VarExpVss = (Cmd/Cdd)*ones(ND,1);
        MapErr = zeros(NX*NY,1); % keeping it in the linear form
        MapErr(MapErr==0)=nan;
        MapErr(IsOcean) = NME; MapErr = reshape(MapErr,NX,NY);
        coverage_percent=100-100*sum(NME,'all')/length(IsOcean);
        % old stuff
        gamma=all_covs./all_covs(max_ij(1),max_ij(2)); %normalize covariances by variance of float location
        QQtmp=round(QQ-gamma.*QQ(max_ij(1),max_ij(2),:),3,'significant'); %remove influence of float
        %float_QQ(:,:,:,iter)=QQtmp;
        QQ=QQtmp; %use the new pCO2 without the influence of float
        count=0; %reinitialize count for next run
        max_cov=1e-9; %reinitialize maximum covariance for next run
        %plot_var=plot_var+gamma.*oneQQ; %add onto plot from last iteration
        %float_locations(iter,:)=max_ij;
        %figure()
        %imagescnan(x(imin:igap:imax),y(jmin:jgap:jmax),plot_var(imin:igap:imax,jmin:jgap:jmax)');axis xy; a=colorbar; %caxis([-1,20]); % creates 2d colormap 
        %pause;
        %oldmax_ij=max_ij; %for making sure future float locations are different
        PlotMapErr(:,:,iter)=MapErr;
        PlotCoverage(iter)=coverage_percent;
    end
end
toc

% Plotting the uncertainty 
cd /home/winnie.chu/figures % goes to the right directory to save file
v = VideoWriter('pCO2_20floats.avi','Uncompressed AVI');
v.FrameRate = 1; %set to 3 frames per second
open(v)

for iter=1:num_iter
    figure()
    title(sprintf('percent coverage of %d floats is %.2f%%',iter,PlotCoverage(iter)))
    tmpplot=PlotMapErr(:,:,iter);
    imagescnan(x(imin:igap:imax),y(jmin:jgap:jmax),tmpplot'); axis xy; a=colorbar;
    caxis([0,1])
    xlabel('longitude');
    ylabel('latitude');
    set(gcf, 'Color', 'w')
    hold on
    text(x(float_locations(1:iter,1)),y(float_locations(1:iter,2)),'\otimes','Color','white','FontSize',16,'HorizontalAlignment','center','FontWeight','bold')
    hold off
    drawnow
    frame= getframe(gcf);
    writeVideo(v,frame);
    close;
end

close(v);

% Running the floats code
% in order to get the mapping data for the floats
cd /home/winnie.chu/functions
tic
[PlotMapErr,PlotCoverage,float_locations,Irange,Jrange]=find_floats(100);
toc


%Covariance plotting (LATEST WORKING VERSION)
function [PlotMapErr,PlotCoverage,float_locations,Irange,Jrange]=find_floats(num_iter)
%% for absolute best float locations
% LAST WORKING VERSION USE THIS ONE 

addpath /home/mmazloff/ANALYSIS  %this has imagescnan and rdmds
addpath /home/mmazloff/ANALYSIS/export_fig/ %this has export fig
cd /data/SO6/TPOSE/bgc_tpose6/2004_fwd_run/diags %new model
load  /data/SO6/TPOSE/bgc_tpose6/grid/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC
cd /home/winnie.chu/data
load pCO2update.mat pCO2res
load DCall_1deg_update.mat DCL

%% initialization of variables
%num_iter=15; % how many floats to find
xk = x*deg2km(1);yk = y*deg2km(1); % km version of grid
% range in longitude-direction
% imin=276;
% imax=996; 
% igap=6;
imin=1;
imax=1128;
igap=6;
Irange=imin:igap:imax;
% range in latitude-direction
jmin=1;
jmax=336; 
jgap=6;
Jrange=jmin:jgap:jmax;
PlotMapErr=zeros(length(Irange),length(Jrange),num_iter); % for plotting error map
PlotCoverage=zeros(1,num_iter); % for recording float coverage of field
float_ind=zeros(1,num_iter); % linear indices of float locations
float_locations=zeros(num_iter,2); % index in 2d of float locations
best_gain=0; % for recursion to compare information gained by extra float
NoiseValue = 1e-10; % constant representing noise of float data

%% retrieving DCL (decorrelation lengths matrix)
% making modifications so that the areas without a convergent gaussian fit
% are omitted and the average is taken from the surroundings
% semimajor lengths
% major=real(DCL(:,:,1));
% major(isnan(major)==1)=0;
% major(major==0)=mean(major(major>0),'all');
% % semiminor lengths
% minor=real(DCL(:,:,2));
% minor(isnan(minor)==1)=0;
% minor(minor==0)=mean(minor(minor>0),'all');
% % angle of rotation
% theta=real(DCL(:,:,3));
% theta(isnan(theta)==1)=0;

major=DCL(:,:,1);
minor=DCL(:,:,2);
theta=DCL(:,:,3);

[row, col] = find(imag(major)~=0); % finding the weird points
%length(row)/(size(major,1)*size(major,2)) % proportion of points that are weird

% first make all the complex areas nan
for err_ind = 1:length(row)
    major(row(err_ind),col(err_ind))=nan;
    minor(row(err_ind),col(err_ind))=nan;
    theta(row(err_ind),col(err_ind))=nan;
end

% then replace the complex areas with the mean of all the real lengths
% surrounding it 
for err_ind = 1:length(row)
    if row(err_ind) <= 60
        if col(err_ind) <= 60
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind):row(err_ind)+60,col(err_ind):col(err_ind)+60),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind):row(err_ind)+60,col(err_ind):col(err_ind)+60),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind):row(err_ind)+60,col(err_ind):col(err_ind)+60),'all','omitnan');
        elseif col(err_ind) > 60 && col(err_ind) <= 276
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind):row(err_ind)+60,col(err_ind)-60:col(err_ind)+60),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind):row(err_ind)+60,col(err_ind)-60:col(err_ind)+60),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind):row(err_ind)+60,col(err_ind)-60:col(err_ind)+60),'all','omitnan');
        elseif col(err_ind) > 276
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind):row(err_ind)+60,col(err_ind)-60:col(err_ind)),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind):row(err_ind)+60,col(err_ind)-60:col(err_ind)),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind):row(err_ind)+60,col(err_ind)-60:col(err_ind)),'all','omitnan');
        end
    elseif row(err_ind) > 60 && row(err_ind) <= 1068
        if col(err_ind) <= 60
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind)-60:row(err_ind)+60,col(err_ind):col(err_ind)+60),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind)-60:row(err_ind)+60,col(err_ind):col(err_ind)+60),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind)-60:row(err_ind)+60,col(err_ind):col(err_ind)+60),'all','omitnan');
        elseif col(err_ind) > 60 && col(err_ind) <= 276
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind)-60:row(err_ind)+60,col(err_ind)-60:col(err_ind)+60),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind)-60:row(err_ind)+60,col(err_ind)-60:col(err_ind)+60),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind)-60:row(err_ind)+60,col(err_ind)-60:col(err_ind)+60),'all','omitnan');
        elseif col(err_ind) > 276
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind)-60:row(err_ind)+60,col(err_ind)-60:col(err_ind)),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind)-60:row(err_ind)+60,col(err_ind)-60:col(err_ind)),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind)-60:row(err_ind)+60,col(err_ind)-60:col(err_ind)),'all','omitnan');
        end
    elseif row(err_ind) > 1068
        if col(err_ind) <= 60
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind)-60:row(err_ind),col(err_ind):col(err_ind)+60),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind)-60:row(err_ind),col(err_ind):col(err_ind)+60),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind)-60:row(err_ind),col(err_ind):col(err_ind)+60),'all','omitnan');
        elseif col(err_ind) > 60 && col(err_ind) <= 276
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind)-60:row(err_ind),col(err_ind)-60:col(err_ind)+60),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind)-60:row(err_ind),col(err_ind)-60:col(err_ind)+60),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind)-60:row(err_ind),col(err_ind)-60:col(err_ind)+60),'all','omitnan');
        elseif col(err_ind) > 276
            major(row(err_ind),col(err_ind)) = mean(major(row(err_ind)-60:row(err_ind),col(err_ind)-60:col(err_ind)),'all','omitnan');
            minor(row(err_ind),col(err_ind)) = mean(minor(row(err_ind)-60:row(err_ind),col(err_ind)-60:col(err_ind)),'all','omitnan');
            theta(row(err_ind),col(err_ind)) = mean(theta(row(err_ind)-60:row(err_ind),col(err_ind)-60:col(err_ind)),'all','omitnan');
        end
    end
end


%% retriving pCO2 data
QQ=pCO2res; % 1128 x 336 x 168 matrix of x, y and time data
Qss = QQ(Irange,Jrange,:); % subsampling of pCO2 field to save space
[NX,NY,NT] = size(Qss); % getting the dimensions of the compressed matrix
Qssv = reshape(Qss,NX*NY,NT); % turning the xy into one dimension vs. time
[IsOcean,~] = find(isnan(Qssv(:,1))==0); % finding all the indices that are not land
Qsssv = Qssv(IsOcean,:); % compressing further to only ocean points NM x NT
NM = length(IsOcean); % new length of xy that is on ocean
normQv = Qsssv-mean(Qsssv,2); % subtract the average pCO2 value from time series at each point
% newQv = normQv;
newQv = Qsssv; % newQv needed for iterations

%% building Gaussian distances using decorrelation lengths and angles
DIST=zeros(NM,NM); % gaussian distance matrix for all ocean points in Irange Jrange
for ind1=1:NM % in paper ind 1 = location K 
    [I,J]=ind2sub([NX NY],IsOcean(ind1));
    I=round((I-1)*igap+imin);
    J=round((J-1)*jgap+jmin);
    dca=major(I,J); % semimajor length of Gaussian at point IJ
    dcb=minor(I,J); % semiminor length of Gaussian at point IJ
    dctheta=theta(I,J); % angle of rotation of Gaussian
    for ind2=1:NM % in  paper ind 2 = location k 
        [i,j]=ind2sub([NX NY],IsOcean(ind2));
        i=round((i-1)*igap+imin);
        j=round((j-1)*jgap+jmin);
        a = (xk(I) - xk(i));
        dx = a*cosd((y(J) + y(j))*.5); % scaled distance in X
        dy = yk(J) - yk(j);
        DIST(ind1,ind2)=exp(-((dx.*cos(dctheta)+dy.*sin(dctheta))^2/dca.^2+(dx.*sin(dctheta)+dy.*cos(dctheta))^2/dcb.^2)); %gaussian distances
    end
end

%% recursive method to find best locations

% PlotMapErr=zeros(length(Irange),length(Jrange),num_iter); % for plotting error map
% PlotCoverage=zeros(1,num_iter); % for recording float coverage of field
% float_ind=zeros(1,num_iter); % linear indices of float locations
% float_locations=zeros(num_iter,2); % index in 2d of float locations
% best_gain=0; % for recursion to compare information gained by extra float
% %newQv = normQv; % newQv needed for iterations
% newQv = normQv;

for iter=1:num_iter
    best_gain=0;
    iter
    for test_ind=1:NM
        Cmnew=newQv*newQv'; % covariance of model to model points
        Dtmp = newQv(test_ind,:); % time series of test location for float
        Cddtmp = Dtmp*Dtmp'+NoiseValue; % variance of time series of test location with added noise
        Cmdtmp = newQv*Dtmp'.*DIST(:,test_ind); % covariance of model and test location
        gaintmp = diag(Cmdtmp*(Cddtmp)^-1*(Cmdtmp')); %gain matrix for each model point
        gain=sum(gaintmp,'all'); % integrated gain total
        if gain > best_gain
            best_float_ind=test_ind;
            best_gain=gain;
        end
    end
    float_ind(iter)=best_float_ind; % record linear index of best float location
    [max_i,max_j]=ind2sub([NX NY], IsOcean(best_float_ind)); % convert linear index to 2d coordinates
    max_i=round((max_i-1)*igap+imin); % project 2d coordinates to Irange 
    max_j=round((max_j-1)*jgap+jmin); % project 2d coordinates to Jrange 
    float_locations(iter,:)=[max_i,max_j]; % best float location in 2d
    gamma=Cmnew(:,best_float_ind)/Cmnew(best_float_ind,best_float_ind); %regression coefficient of covariance 
    newQv=newQv-gamma.*newQv(best_float_ind,:); % remove influence of best float from all model time series
    % iterate through again! 
end


%% creating matrices to plot the locations and coverage of floats

Cmm = normQv*normQv'.*DIST; % model to model covariance localized by gaussian distance
for iter=1:num_iter
    floats=float_ind(1:iter); % simulate iterative deployment of floats
    ND=length(floats); % number of floats currently deployed 
    D = normQv(floats,:); % getting just the float points
    Cmd = normQv*(D').*DIST(:,floats); % model to float covariance localized by gaussian distance
    Cdd = D*(D').*DIST(floats,floats); % float to float covariance localized by gaussian distance
    Cdd = Cdd + NoiseValue*eye(ND); % adding noise to the variance of float locations
    FME = diag(Cmm) - diag(Cmd*(Cdd)^-1*(Cmd')); % formal mapping error
    NME = FME./diag(Cmm); % normalized mapping error 0 to 1
    MapErr = zeros(NX*NY,1); % in linear form
    MapErr(MapErr==0)=nan; % start with all NaN for areas that are land
    MapErr(IsOcean) = NME; % populate ocean points with normalized mapping error
    MapErr = reshape(MapErr,NX,NY); % turn into 2d map 
    coverage_percent=100-100*sum(NME,'all')/NM; % total amount of certainty from float observations
    PlotMapErr(:,:,iter)=MapErr; % saving mapping error in matrix to be used in plot_float_coverage function
    PlotCoverage(iter)=coverage_percent; 
end

cd /home/winnie.chu/functions
end
%% test against random selection 

%float_ind=randi(NM,1,num_iter);

% to find the index corresponding to latitude/longitude
min(find(abs(x-250)==min(abs(x-250)))); % 250 is long (E) we are trying to find
min(find(abs(y-18)==min(abs(y-18)))); % 18 is lat (N) we are trying to find



