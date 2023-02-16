function [PlotMapErr,PlotCoverage,float_locations,waypoints,Irange,Jrange]=float_boat(num_iter, wp1, wp2)
% LAST WORKING VERSION

% need to provide 2 wavepoints for a predetermined cruise path
% wavepoints are in [lat long] format
% latitude/longtiude of wavepoints
%wp1_lat=-5; wp2_lat=10; wp1_long=280; wp2_long=140;
%y=-20:2:20;
%x=0:10:300;
addpath /home/mmazloff/ANALYSIS  %this has imagescnan and rdmds
addpath /home/mmazloff/ANALYSIS/export_fig/ %this has export fig
cd /data/SO6/TPOSE/bgc_tpose6/2004_fwd_run/diags %new model
load  /data/SO6/TPOSE/bgc_tpose6/grid/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC
cd /home/winnie.chu/data
load pCO2update.mat pCO2res
load DCall_1deg_update.mat DCL

%% initialization of variables
xk = x*deg2km(1);yk = y*deg2km(1); % km version of grid
% range in longitude-direction
imin=157;
imax=1128;
igap=6;
Irange=imin:igap:imax;
% range in latitude-direction
jmin=37;
jmax=277; 
jgap=6;
Jrange=jmin:jgap:jmax;
PlotMapErr=zeros(length(Irange),length(Jrange),num_iter); % for plotting error map
PlotCoverage=zeros(1,num_iter); % for recording float coverage of field
float_ind=zeros(1,num_iter); % linear indices of float locations
float_locations=zeros(num_iter,2); % index in 2d of float locations
best_gain=0; % for recursion to compare information gained by extra float
NoiseValue = 1e-10; % constant representing noise of float data

%% retrieving DCL (decorrelation lengths matrix)
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
newQv = normQv; % newQv needed for iterations

%% find waypoints from lat/long start and endpoints input
wp1_lat=wp1(1); wp1_long=wp1(2); wp2_lat=wp2(1); wp2_long=wp2(2);
wp1_y = Jrange(min(find(abs(y(Jrange)-wp1_lat)==min(abs(y(Jrange)-wp1_lat)))));
wp2_y = Jrange(min(find(abs(y(Jrange)-wp2_lat)==min(abs(y(Jrange)-wp2_lat)))));
wp1_x = Irange(min(find(abs(x(Irange)-wp1_long)==min(abs(x(Irange)-wp1_long)))));
wp2_x = Irange(min(find(abs(x(Irange)-wp2_long)==min(abs(x(Irange)-wp2_long)))));
waypoints = [wp1_x, wp1_y; wp2_x, wp2_y];

% finding the equation of line connecting two wavepoints
% coefficients = polyfit([wp1_x, wp2_x], [wp1_y, wp2_y], 1);
% a = coefficients (1);
% b = coefficients (2);
% where y=a*x+b
% plot([x(wp1_x), x(wp2_x)], [y(wp1_y), y(wp2_y)])
% assumes straight line between waypoints rather than curvature of earth

% getting two arrays describing x and y along line
if wp1_x<wp2_x
    line_x=wp1_x:wp2_x;
    if wp1_y==wp2_y
        line_y=wp1_y*ones(1,length(line_x));
    else
        coefficients = polyfit([wp1_x, wp2_x], [wp1_y, wp2_y], 1);
        a = coefficients (1);
        b = coefficients (2);
        line_y=round(a*line_x+b);
    end
elseif wp1_x>wp2_x
    line_x=wp2_x:wp1_x;
    if wp1_y==wp2_y
        line_y=wp1_y*ones(1,length(line_x));
    else
        coefficients = polyfit([wp1_x, wp2_x], [wp1_y, wp2_y], 1);
        a = coefficients (1);
        b = coefficients (2);
        line_y=round(a*line_x+b);
    end
elseif wp1_x==wp2_x
    if wp1_y<wp2_y
        line_y=wp1_y:wp2_y;
        line_x=wp1_x*ones(1,length(line_y));
    elseif wp2_y<wp1_y
        line_y=wp2_y:wp1_y;
        line_x=wp1_x*ones(1,length(line_y));
    end
end
%line_y=round(a*line_x+b);
line_xy=cat(1,line_x,line_y); % is a 2xlength(line_x)vector where row1=x,row2=y
% need decorrelation matrix for all points within the line

% find all points in line within Irange and Jrange
line_xyss=line_xy(:,find(ismember(line_x,Irange)==1)); % contains x points from line
line_xyss=line_xyss(:,find(ismember(line_xyss(2,:),Jrange)==1)); % contains y points from line too

% convert to subscripts using Irange Jrange
line_xyss(1,:)=(line_xyss(1,:)-imin)/igap + 1;
line_xyss(2,:)=(line_xyss(2,:)-jmin)/jgap + 1;

% convert to linear indices
line_xyss=sub2ind([NX,NY],line_xyss(1,:),line_xyss(2,:));

% convert to linear indices from IsOcean
line_xyss=find(ismember(IsOcean,line_xyss)==1);

% total 
NL=length(line_xyss);

%% building Gaussian distances using decorrelation lengths and angles
DIST=zeros(NM,NM); % gaussian distance matrix for all ocean points in Irange Jrange
for ind1=1:NM
    [I,J]=ind2sub([NX NY],IsOcean(ind1));
    I=round((I-1)*igap+imin);
    J=round((J-1)*jgap+jmin);
    dca=major(I,J); % semimajor length of Gaussian at point IJ
    dcb=minor(I,J); % semiminor length of Gaussian at point IJ
    dctheta=theta(I,J); % angle of rotation of Gaussian
    for ind2=1:NM
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
% newQv = normQv; % newQv needed for iterations

for iter=1:num_iter
    best_gain=0;
    iter
    for line_ind=1:NL
        test_ind=line_xyss(line_ind);
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

%% OLD STUFF
% CC=zeros(size(QQ,1),size(QQ,2));
% float_I=[];
% float_J=[];
% 
% %MATT’S ATTEMPT AT Cmm
% %SUBSAMPLE DATA
% Qss = QQ(imin:igap:imax,jmin:jgap:jmax,:);
% [NX,NY,NT] = size(Qss); % getting the dimensions of the compressed matrix
% Qssv = reshape(Qss,NX*NY,NT); % turning the xy into one dimension vs. time
% [IsOcean,~] = find(isnan(Qssv(:,1))==0); % finding all the indices that are not land
% Qsssv = Qssv(IsOcean,:); % compressing further to only ocean points
% %NM = length(IsOcean); % new length of xy that is on ocean
% Cmm = Qsssv*(Qsssv'); % Cmm is the local variances size should be NMxNM
% 
% %TO VISUALIZE:
% % Qvar = zeros(NX*NY,1); Qvar(IsOcean) = diag(Cmm); Qvar = reshape(Qvar,NX,NY);
% % Qvar is the local variances of each of the ocean points in the model
% %NEXT IS TO FIND Cmd which is covariance at data points
% %Qsssv is signal at every model point
% 
% %HERE GOES WINNIE’S LOOP TO FIND THE DATA LOCATIONS I AND J
% % starting the actual iteration
% tic
% for iter=1:num_iter
%     max_ij=[0,0]; %initializing location of first float
%     count=0; %initializing count
%     for line=1:length(line_xyss)
%         i=line_xyss(1,line);
%         j=line_xyss(2,line);
%         dcx=DC(i,j,1);
%         dcy=DC(i,j,2);
%         tmp1=squeeze(QQ(i,j,:)); %getting the time series of first point
%         if isnan(tmp1(1))==1
%             %count=0; 
%             continue
%         elseif count > 0
%             %cov_ij=sum(abs(CC),'all','omitnan');
%             %if line==1
%                 %cov_ij=sum(CC.*CC,'all','omitnan')/var(QQ(i,j,:));
%                 %test_repeat=ismember([i,j],float_locations,'rows'); % check if float location is already found
%             %else
%                 prev_i=line_xyss(1,line-1);
%                 prev_j=line_xyss(2,line-1);
%                 cov_ij=sum(CC.*CC,'all','omitnan')/var(QQ(prev_i,prev_j,:));
%                 test_repeat=ismember([prev_i,prev_j],float_locations,'rows'); % check if float location is already found
%             if cov_ij > max_cov && test_repeat==0
%             %if cov_ij > max_cov
%                 % if not an existing float location and bigger
%                 % covariance
%                 max_cov=cov_ij; %store the highest covariances
%                 all_covs=CC; %store covariances over entire area
%                 max_ij=[prev_i,prev_j] %store float location
%                 %CC=zeros(size(QQ,1),size(QQ,2));
%             end
%         end
%         count=count+1;
%         for I=imin:igap:imax
%             for J=jmin:jgap:jmax
%                 tmp2=squeeze(QQ(I,J,:)); %get the time series for 2nd point
%                 if isnan(tmp2(1))==1
%                     continue
%                 end
%                 a = (xk(I) - xk(i));
%                 dx = a*cosd((y(J) + y(j))*.5);% scaled distance in X 
%                 dy =  yk(J) - yk(j);
%                 % (admittedly large) radius of influence
%                 %if (dx.^2 + dy.^2) < 1000000 %max distance
%                 %if dx <= 1500 && dy <= 350 && (dx.^2 + dy.^2) < 1000000
%                 % remember to half dx and dy if you're using
%                 % correlation lengths
%                 covmatrix=cov(tmp1,tmp2);
%                 CC(I,J)=covmatrix(2,1)*exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))); %covariance of point ij with point IJ
%                 % ASK IF THIS SHOULD DCX SHOULD BE DOUBLED FOR THE OVERALL
%                 % WIDTH OF GAUSSIAN OR IF SEMIMAJOR/SEMIMINOR LENGTH IS
%                 % GOOD
%             end
%         end
%     end
%     %TO TEST
%     %I = round(([390 492 948 324 846] - imin)/6);
%     %J = round(([180 174 168 162 174] - jmin)/6);
%     %ND = 5;
%     %NOW PLOT VarUnexp
% %     figure()
% %     title(sprintf('percent coverage of %d floats is %.2f%%',iter,coverage_percent))
% %     imagescnan(x(imin:igap:imax),y(jmin:jgap:jmax),MapErr'); axis xy; a=colorbar;
% %     caxis([0,1])
% %     hold on
% %     text(x(float_locations(:,1)),y(float_locations(:,2)),'\otimes','Color','white','FontSize',16,'HorizontalAlignment','center','FontWeight','bold')
% %     hold off
%     if max_ij(1)==0
%         % no max_ij means we found nothing to report
%         continue
%     else 
%         float_locations(iter,:)=max_ij;
%         % oof this won't work when max_ij isn't within imin:igap:imax
%         float_I=[float_I round((max_ij(1)-imin)/igap+1)]; %getting the float location indices in compressed QQ
%         float_J=[float_J round((max_ij(2)-jmin)/jgap+1)];
%         ND=length(float_I);
%         %NOW THIS GOES AFTER YOUR ITERATIVE LOOP. RUN EACH TIME TO GET “FORMAL MAP ERROR”
%         %Need signal just at every data point
%         MASK = zeros(NX,NY);
%         DIST = zeros(NX,NY,ND);
%         for num = 1:ND %NDATA
%             MASK(float_I(num),float_J(num)) = 1;
%         end
%         for num = 1:ND
%             I=round((float_I(num)-1)*igap+imin);
%             J=round((float_J(num)-1)*jgap+jmin);
%             dcx=DC(I,J,1);
%             dcy=DC(I,J,2);
%             for nx=1:NX
%                 for ny=1:NY
%                     i=round(nx*igap+imin);
%                     j=round(ny*jgap+jmin);
%                     a = (xk(I) - xk(i));
%                     dx = a*cosd((y(J) + y(j))*.5);% scaled distance in X
%                     dy = yk(J) - yk(j);
%                     %DIST(nx,ny,num)=round(exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))),3,'significant'); %gaussian distances
%                     DIST(nx,ny,num)=exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))); %gaussian distances
%                 end
%             end
%         end
%         DIST=reshape(DIST,NX*NY,ND); % turn distances into a linear vector for each float
%         %DISTsss=DIST(IsOcean,:); %only needs floats in ocean should be lenght of IsOcean x ND
%         MASK = MASK(:); % turned into a linear vector rather than matrix
%         %MASKsss = MASK(IsOcean); %only need the points in ocean
%         %[Is] = find(MASKsss==1); %where the floats are
%         %MASKsss = MASK(IsOcean); %only need the points in ocean
%         [Is] = find(MASK==1); %where the floats are
%         D = Qssv(Is,:); % getting just the float points
%         float_ind = sub2ind(size(Qss(:,:,1)),float_I,float_J);
%         DDist =DIST(float_ind,:); %use float_ind instead of IS to make symmetric dist
%         %Cmd = Qsssv*(D'); % finding the covariance between all those float points and the model
%         Cmd = Qssv*(D').*DIST; % number of NXNY points that are ocean x number of floats
%         %NOW GET DATA DATA COVARIANCE
%         Cmd = Cmd(IsOcean,:);
%         Cdd = D*(D').*DDist; % get covariance between float locations
%         %ADD NOISE TO DATA BECAUSE FLOATS AREN’T PERFECT
%         NoiseValue = 1e-5; %Noise variance in obs keep in mind order of magnitude of pCO2 values is around 10^-4
%         %NoiseValue = 0;
%         Cdd = Cdd + NoiseValue*eye(ND); % adding noise to the variance of float locations
%         %FORMAL MAPPING ERROR
%         %FME = diag(Cmm) - diag(Cmd*(Cdd)^-1*Cmd');
%         FME = diag(Cmm) - diag(Cmd*(Cdd)^-1*(Cmd'));
%         %NORMALIZED VERSION
%         NME = FME./diag(Cmm);
%         %%OLD   %VarExpVss = (Cmd/Cdd)*ones(ND,1);
%         MapErr = zeros(NX*NY,1); % keeping it in the linear form
%         MapErr(MapErr==0)=nan;
%         MapErr(IsOcean) = NME; MapErr = reshape(MapErr,NX,NY);
%         coverage_percent=100-100*sum(NME,'all')/length(IsOcean);
%         % old stuff
%         gamma=all_covs./var(QQ(max_ij(1),max_ij(2),:)); %normalize covariances by variance of float location
%         %QQtmp=round(QQ-gamma.*QQ(max_ij(1),max_ij(2),:),3,'significant'); %remove influence of float
%         QQ=QQ-gamma.*QQ(max_ij(1),max_ij(2),:); %remove influence of float
%         %float_QQ(:,:,:,iter)=QQtmp;
%         %QQ=QQtmp; %use the new pCO2 without the influence of float
%         %count=0; %reinitialize count for next run
%         %max_cov=1e-9; %reinitialize maximum covariance for next run
%         max_cov=0;
%         CC=zeros(size(QQ,1),size(QQ,2));
%         %plot_var=plot_var+gamma.*oneQQ; %add onto plot from last iteration
%         %float_locations(iter,:)=max_ij;
%         %figure()
%         %imagescnan(x(imin:igap:imax),y(jmin:jgap:jmax),plot_var(imin:igap:imax,jmin:jgap:jmax)');axis xy; a=colorbar; %caxis([-1,20]); % creates 2d colormap 
%         %pause;
%         %oldmax_ij=max_ij; %for making sure future float locations are different
%         PlotMapErr(:,:,iter)=MapErr;
%         PlotCoverage(iter)=coverage_percent;
%     end
% end
% toc
% cd /home/winnie.chu/functions
% end
% 
