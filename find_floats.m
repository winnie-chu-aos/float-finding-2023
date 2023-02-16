function [PlotMapErr,PlotCoverage,float_locations,Irange,Jrange]=find_floats(num_iter)
%% for absolute best float locations

addpath /home/mmazloff/ANALYSIS  %this has imagescnan and rdmds
addpath /home/mmazloff/ANALYSIS/export_fig/ %this has export fig
cd /data/SO6/TPOSE/bgc_tpose6/2004_fwd_run/diags/ %old model
load  /data/SO6/TPOSE/bgc_tpose6/grid/grid XC YC Depth; 
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
