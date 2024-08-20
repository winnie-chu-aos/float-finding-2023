
%% Loading pCO2 data
pCO2 = zeros(1128,336,189,'single');           % 1128 longitude pts, 336 latitude pts, 189 time steps (may change depending on iteration of TPOSE model)
pCO2= rdmds(['diag_surf'],nan,'rec',2);        % read files with rdmds fxn
pCO2(pCO2==0)=nan;      
pCO2 = pCO2(:,:,13:180);                       % crop first year (make sure the number of total months is divisible by 12)
nmonths = size(pCO2,3);

%% Data preprocessing
% initializing paths and arrays
load  /data/SO6/TPOSE/bgc_tpose6_k51/grid/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC                           % get coordinates in x and y
NT = nmonths; 
NX = 1128;
NY = 336;
Time = [1:365.25/12:365.25*NT/12];Time = Time' - mean(Time); 
pCO2res = zeros(1128,336,nmonths,'single');                     % initialize residual pCO2 matrix
pCO2VE = zeros(1128,336,5,'single');                            % initialize variance explained pCO2 matrix

% simple pCO2 calculations + wavenumber stuff
pCO2_bar = mean(pCO2,3);                                      
pCO2_std = std(pCO2,[],3);                                      
pCO2cos = zeros(NX,NY,4,'single'); pCO2sin = pCO2cos;
pvar = 1;dvar= .01;ifwd = 1;
kx1 = (2*pi)*[1]/[365.25];                                      % wavenumbers to fit to: 1 yr fit to pCO2
kx2 = (2*pi)*[2]/[365.25];                                      % wavenumbers to fit to: 6 mnt fit to pCO2
kx3 = (2*pi)*[3]/[365.25];                                      % wavenumbers to fit to: 4 mnt fit to pCO2
kx4 = (2*pi)*[4]/[365.25];                                      % wavenumbers to fit to: 3 mnt fit to pCO2
TimeT = Time;
clear Time 

% loop for variance and residual signal
for j = 1:NY
    NY-j
    for i = 1:NX
      tmp2 = squeeze(pCO2(i,j,:));tmp2a = detrend(tmp2 - mean(tmp2)); % get rid of trend in pCO2
      [p2,y2] = gft(TimeT,tmp2a,[kx1 kx2 kx3 kx4]',pvar,dvar,ifwd);   % fourier transform fxn gft - least squares fit of 4 harmonics to time series
      pCO2cos(i,j,:) = p2(1:4);                                       % cosin amp of kx1
      pCO2sin(i,j,:) = p2(5:8);                                       % sin amp of kx1
      y11 = p2(1)*cos(TimeT*kx1)+p2(5)*sin(TimeT*kx1);                % annual cycle
      y12 = p2(2)*cos(TimeT*kx2)+p2(6)*sin(TimeT*kx2);                % semi annual
      y13 = p2(3)*cos(TimeT*kx3)+p2(7)*sin(TimeT*kx3);                % 4 month
      y14 = p2(4)*cos(TimeT*kx4)+p2(8)*sin(TimeT*kx4);                % 3 month
      pCO2VE(i,j,1) = 100 - 100*var(tmp2-y11)/var(tmp2);              % variance explained by annual cycle
      pCO2VE(i,j,2) = 100 - 100*var(tmp2-y12)/var(tmp2);              % variance explained by semi annual
      pCO2VE(i,j,3) = 100 - 100*var(tmp2-y13)/var(tmp2);              % variance explained by 4 month
      pCO2VE(i,j,4) = 100 - 100*var(tmp2-y14)/var(tmp2);              % variance explained by 3 month
      pCO2VE(i,j,5) = 100 - 100*var(y2)/var(tmp2);                    % res var, how much variance is leftover from original signal
      pCO2res(i,j,:) = tmp2 - y2;                                     % time series residual, outside of trends/cycles
      pCO2H_std(i,j) = std(y2);
      pCO2R_std(i,j) = std(tmp2 - y2);
   end
end
clear tmp*;

% getting phase of cycles
cycle_phases_y11 = zeros(1128,336,'single'); 
cycle_phases_y12=cycle_phases_y11; cycle_phases_y13=cycle_phases_y11; cycle_phases_y14=cycle_phases_y11;

for j = 1:NY
    NY-j
    for i = 1:NX
      tmp2 = squeeze(pCO2(i,j,:));tmp2a = detrend(tmp2 - mean(tmp2)); % remember to crop first year
      [p2,y2] = gft(TimeT,tmp2a,[kx1 kx2 kx3 kx4]',pvar,dvar,ifwd);   % gft in home/ANALYSIS-- least squares fit of 4 harmonics to time series
      pCO2cos(i,j,:) = p2(1:4);                                       % cosin amp of kx1
      pCO2sin(i,j,:) = p2(5:8);                                       % sin amp of kx1
      y11 = p2(1)*cos(TimeT*kx1)+p2(5)*sin(TimeT*kx1);                % annual cycle
      y12 = p2(2)*cos(TimeT*kx2)+p2(6)*sin(TimeT*kx2);                % semi annual
      y13 = p2(3)*cos(TimeT*kx3)+p2(7)*sin(TimeT*kx3);                % 4 month
      y14 = p2(4)*cos(TimeT*kx4)+p2(8)*sin(TimeT*kx4);                % 3 month
	cycle_phases_y11(i,j)= atan2(p2(5),p2(1));
	cycle_phases_y12(i,j)= atan2(p2(6),p2(2));
	cycle_phases_y13(i,j) = atan2(p2(7),p2(3));
    cycle_phases_y14(i,j) = atan2(p2(8),p2(4));
   end
end
clear tmp*;

% saving phase plot for annual cycle
test_cycle_phases_y11= cycle_phases_y11*(12/(2*pi)); % maximum pCO2
test_cycle_phases_y11(test_cycle_phases_y11<=0)=test_cycle_phases_y11(test_cycle_phases_y11<=0)+12; % maximum pCO2

%% Correlation length estimation
% obtaining correlation map for a single point

xk = x*deg2km(1);yk = y*deg2km(1); 
NX = 1128; NY = 336; NT = 180;                                                    % number of position and time points
amp=1;
dstx = 240;                                                                       % max grid points to look in x
dsty = 72;                                                                        % max grid points to look in y matching aspect ratio of grid 1128x336
Q=pCO2res;                                                                        % use pCO2 residual signal
QQ = Q(:,:,1:6)*nan;                                                              % initialize matrix storing correlation lengths
NoiseFac = 10;
SNR = 1/NoiseFac^2;                                                               % signal to noise ratio
NoiseFac1=1;NoiseFac2 = 1;
NSEC = 18;                                                                        % only go out to a threshold, do in 360/NSEC deg sectors
thrshcntr = .8;
tic 

for i=276:12:996                                                                  % between 150 to 270 degrees E
    i
    for j=90:12:210                                                               % between -10 to 10 degrees N
        j
	    loop=0;
        mjr=0;                                                                    % throwaway
        mnr=0;                                                                    % throwaway
        errx=550;
        erry=550;
        while errx>0.01 && erry>0.01 && loop<6                                    % error tolerance is 0.01 and convergence is 6 loops or less
            prev_mjr=mjr;
            prev_mnr=mnr;
            if loop==0
                    LXU=200; LYU=200;                           
            else
                    LXU=20; LYU=20;
            end
            A_p = [0.5*(LXU^-2) 0 1.5*(LYU^-2)]';                                 % a priori guess
            P_p = [A_p(1)^2 A_p(1)*A_p(3) A_p(3)^2]';
            Pinv = diag(P_p.^-1);                                                 % how much we trust the prior, weird bc this means the larger the uncertainty the more we trust
            I = i-dstx:i+dstx; I(I>NX)=[]; I(I<1)=[]; 
            J = j-dsty:j+dsty; J(J>NY)=[]; J(J<1)=[];
            snx = length(I);
            sny = length(J);
            DC = nan(snx*sny,1,'single'); X = DC; Y = DC;                         % decorrelation field initialization
            X2 = linspace(-4200,4200,snx/1.1); 
            Y2 = linspace(-1500,1500,sny/0.89);
            DC2 = nan(length(X2),length(Y2),'single');
            tmp1 = squeeze(Q(i,j,:));
            count = 0;
            if isfinite(tmp1(1))
                for ii = I
                    for jj = J
                        tmp2 = squeeze(Q(ii,jj,:));
                        if isfinite(tmp2(1))
                            a = (xk(ii) - xk(i));
                            dx = a*cosd((y(jj) + y(j))*.5);                       % scaled distance in X so 1 degree is equal everywhere despite latitude
                            dy =  yk(jj) - yk(j);
                            if (dx.^2 + dy.^2) < 16000000                         % max distance
                                count = count + 1;                
                                X(count) = dx; 
                                DC(count) = corr(tmp1,tmp2); 
                                Y(count) = dy;
             iii = min(find(dx<X2));                                              % put on 2d grid to visualize
             jjj = min(find(dy<Y2));                                              % put on 2d grid to visualize
             DC2(iii,jjj) = DC(count);
                            end
                        end
                    end
                end
                NP = count;
                count = 0;K=0;
                [p,r] = cart2pol(X(1:NP),Y(1:NP)); p = p + pi;                    % map to polar
                p = ceil(p*NSEC/2/pi);                                            % split into discrete sectors, ceil means round up
                [trash,I] = sort(r);                                              % sort so look out in radius, looks for the radius indices
                maxr = ones(1,NSEC)*999999;                                       % init ok radius
        
                for ti = I'
                    if r(ti) < maxr(p(ti))
                        if DC(ti)>thrshcntr                                       % radius used for prior
                            count = count + 1; 
                            K(count) = ti;
                        else
                            maxr(p(ti)) = r(ti);                                  % we are bad, this is new max radius for sector
                        end
                    end 
                end
        %  choose prior based on maxr
                if loop==0                                                        % first loop more uncertainty
                    LXP = 1/sqrt(-2*log(thrshcntr))*mean(maxr([1 2 8:11 17:18])); % zonal
                    LYP = 1/sqrt(-2*log(thrshcntr))*mean(maxr([3:7  12:16]));     % meridional
                else                                                              % guess correlation lenths of previous iteration
                    LXP = prev_mjr;                                               % zonal
                    LYP = prev_mnr;                                               % meridional
                end
                A_p = [0.5*(LXP^-2) 0 0.5*(LYP^-2)]';                             % reassigns prior to the same as before
        
                %NOW FIT TO GAUSSIAN
                K = find((DC>0)==1);                                              % ok as long as NoiseFac2 > 0, points where the correlation is greater than 0
                if length(K)>24
                    ym = -log(DC(K));                   
                    H = [X(K).^2, X(K).*Y(K), Y(K).^2];
                    Rinv  = (NoiseFac1*( 1.01-exp(NoiseFac2*(-A_p(1)*X(K).^2 - A_p(2)*X(K).*Y(K) - A_p(3)*Y(K).^2) ) )).^-1 ;
                    R2inv = repmat(Rinv.^2,1,3);                                  % how much we trust data
                    %least squares fit  
                    Hty = (H.*R2inv)'*(ym - H*A_p);
                    A  = A_p + ((H.*R2inv)'*H + Pinv)\Hty;                        % A_p + R2inv/(R2inv+Pinv)*(ym-A_P) in scalar form 
                    R2  = 1 - var(Rinv.*(exp(-ym) - exp(-H*A)))/var(Rinv.*exp(-ym));
                    [V,L] = eig([A(1), A(2)/2; A(2)/2, A(3)]);
                    mjr = (2*L(1,1))^-.5;
                    mnr = (2*L(2,2))^-.5;
                    thta = atan2(V(1,1),V(2,1)) - pi/2;
                    if thta>pi/2 
                        thta = thta-pi;                                           % rotate 180 so goes from -pi/2 to pi/2
                    elseif thta<-pi/2 
                        thta = thta+pi;
                    end                                                         
                    QQ(i,j,:) = [mjr mnr thta R2 LXP LYP];                        % semimajor lengths (mjr), semiminor lengths (mnr), angle of rotation (thta), r^2 (R2)
                end
            end
            errx=abs(mjr-prev_mjr)/mjr;                                           % calculate error as a percentage of the previous iteration
            erry=abs(mnr-prev_mnr)/mnr;
        loop = loop + 1
        end
    end 
end

t=toc;
%% getting variables ready for float finding algorithm

% initialization of variables
xk = x*deg2km(1);yk = y*deg2km(1);                                      % km version of grid
imin=1;                                                                 % range in longitude-direction
imax=1128;
igap=6;
Irange=imin:igap:imax;
jmin=1;                                                                 % range in latitude-direction
jmax=336; 
jgap=6;
Jrange=jmin:jgap:jmax;
PlotMapErr=zeros(length(Irange),length(Jrange),num_iter);               % for plotting error map
PlotCoverage=zeros(1,num_iter);                                         % for recording float coverage of field
float_ind=zeros(1,num_iter);                                            % linear indices of float locations
float_locations=zeros(num_iter,2);                                      % index in 2d of float locations
best_gain=0;                                                            % for recursion to compare information gained by extra float
NoiseValue = 1e-10;                                                     % constant representing noise of float data

%% retrieving DCL (correlation lengths matrix)
major=DCL(:,:,1);
minor=DCL(:,:,2);
theta=DCL(:,:,3);

[row, col] = find(imag(major)~=0);                                      % finding the weird points
% length(row)/(size(major,1)*size(major,2)) % proportion of points that are weird

% first make all the complex areas nan
for err_ind = 1:length(row)
    major(row(err_ind),col(err_ind))=nan;
    minor(row(err_ind),col(err_ind))=nan;
    theta(row(err_ind),col(err_ind))=nan;
end

% then replace the complex areas with the mean of all the real lengths surrounding it 
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
QQ=pCO2res;                                                             % 1128 x 336 x 168 matrix of x, y and time data
Qss = QQ(Irange,Jrange,:);                                              % subsampling of pCO2 field to save space
[NX,NY,NT] = size(Qss);                                                 % getting the dimensions of the compressed matrix
Qssv = reshape(Qss,NX*NY,NT);                                           % turning the xy into one dimension vs. time
[IsOcean,~] = find(isnan(Qssv(:,1))==0);                                % finding all the indices that are not land
Qsssv = Qssv(IsOcean,:);                                                % compressing further to only ocean points NM x NT
NM = length(IsOcean);                                                   % new length of xy that is on ocean
normQv = Qsssv-mean(Qsssv,2);                                           % subtract the average pCO2 value from time series at each point
newQv = Qsssv;                                                          % newQv needed for iterations

%% building Gaussian distances using decorrelation lengths and angles
DIST=zeros(NM,NM);                                                      % gaussian distance matrix for all ocean points in Irange Jrange
for ind1=1:NM                                                           % in paper ind 1 = location K 
    [I,J]=ind2sub([NX NY],IsOcean(ind1));   
    I=round((I-1)*igap+imin);
    J=round((J-1)*jgap+jmin);
    dca=major(I,J);                                                     % semimajor length of Gaussian at point IJ
    dcb=minor(I,J);                                                     % semiminor length of Gaussian at point IJ
    dctheta=theta(I,J);                                                 % angle of rotation of Gaussian
    for ind2=1:NM                                                       % in  paper ind 2 = location k 
        [i,j]=ind2sub([NX NY],IsOcean(ind2));
        i=round((i-1)*igap+imin);
        j=round((j-1)*jgap+jmin);
        a = (xk(I) - xk(i));
        dx = a*cosd((y(J) + y(j))*.5);                                  % scaled distance in X
        dy = yk(J) - yk(j);
        DIST(ind1,ind2)=exp(-((dx.*cos(dctheta)+dy.*sin(dctheta))^2/dca.^2+(dx.*sin(dctheta)+dy.*cos(dctheta))^2/dcb.^2)); % gaussian distances based on correlation lengths
    end
end

%% recursive method to find best locations

for iter=1:num_iter
    best_gain=0;
    iter
    for test_ind=1:NM
        Cmnew=newQv*newQv';                                             % covariance of model to model points
        Dtmp = newQv(test_ind,:);                                       % time series of test location for float
        Cddtmp = Dtmp*Dtmp'+NoiseValue;                                 % variance of time series of test location with added noise
        Cmdtmp = newQv*Dtmp'.*DIST(:,test_ind);                         % covariance of model and test location
        gaintmp = diag(Cmdtmp*(Cddtmp)^-1*(Cmdtmp'));                   % gain matrix for each model point
        gain=sum(gaintmp,'all');                                        % integrated gain total
        if gain > best_gain
            best_float_ind=test_ind;
            best_gain=gain;
        end
    end
    float_ind(iter)=best_float_ind;                                     % record linear index of best float location
    [max_i,max_j]=ind2sub([NX NY], IsOcean(best_float_ind));            % convert linear index to 2d coordinates
    max_i=round((max_i-1)*igap+imin);                                   % project 2d coordinates to Irange 
    max_j=round((max_j-1)*jgap+jmin);                                   % project 2d coordinates to Jrange 
    float_locations(iter,:)=[max_i,max_j];                              % best float location in 2d
    gamma=Cmnew(:,best_float_ind)/Cmnew(best_float_ind,best_float_ind); %regression coefficient of covariance 
    newQv=newQv-gamma.*newQv(best_float_ind,:);                         % remove influence of best float from all model time series
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

%% Float finding algorithm
% a simple model of the covariance code
prompt = 'Input number of floats deployed: ';
num_iter = input(prompt);
% need to run the above portion first and enter input before subsequent part

% Initialization
QQ=pCO2res;
max_cov=1e-9;
plot_var=zeros(size(QQ,1),size(QQ,2));                                  % initialization
oneQQ=ones(size(QQ,1),size(QQ,2));                                      % ones matrix the same size as QQ
float_locations=zeros(num_iter,2);
xk = x*deg2km(1);yk = y*deg2km(1); 
indices = find(isnan(pCO2(:,:,1)));
plot_var(indices)=nan;                                                  % remove points on land
imin=276;
imax=996;                                                               % range of i points
igap=6;
jmin=96;
jmax=216;                                                               % range of j points
jgap=6;
CC=zeros(size(QQ,1),size(QQ,2));
float_I=[];
float_J=[];
dcx=3500;                                                               % just general estimate of decorrelation lengths in x direction
dcy=1000;                                                               % estimate of decorrelation lengths in y direction
PlotMapErr=zeros(length(imin:igap:imax),length(jmin:jgap:jmax),num_iter);
PlotCoverage=zeros(num_iter);

% Finding model-model, model-float, and float-float covariance
% subsample data
Qss = QQ(imin:igap:imax,jmin:jgap:jmax,:);
[NX,NY,NT] = size(Qss);                                                 % getting the dimensions of the compressed matrix
Qssv = reshape(Qss,NX*NY,NT);                                           % turning the xy into one dimension vs. time
[IsOcean,junk] = find(isnan(Qssv(:,1))==0);                             % finding all the indices that are not land
Qsssv = Qssv(IsOcean,:);                                                % compressing further to only ocean points
NM = length(IsOcean);                                                   % new length of xy that is on ocean
Cmm = Qsssv*(Qsssv');                                                   % Cmm is the local variances size should be NMxNM
    
% To visualize
Qvar = zeros(NX*NY,1); Qvar(IsOcean) = diag(Cmm); Qvar = reshape(Qvar,NX,NY);
% Qvar is the local variances of each of the ocean points in the model
%NEXT IS TO FIND Cmd which is covariance at data points
%Qsssv is signal at every model point

% starting the actual iteration
tic
for iter=1:num_iter
    max_ij=[0,0];                                                       % initializing location of first float
    count=0;                                                            % initializing count
    for i=imin:igap:imax
        for j=jmin:jgap:jmax
            tmp1=squeeze(QQ(i,j,:));                                    % getting the time series of first point
            if isnan(tmp1(1))==1
            elseif count > 0
                cov_ij=sum(CC,'all','omitnan');
                if cov_ij > max_cov 
                    % if not an existing float location and bigger covariance
                    max_cov=cov_ij;                                     % store the highest covariances
                    all_covs=CC;                                        % store covariances over entire area
                    max_ij=[i,j]                                        % store float location
                end
            end
            count=count+1;
            for I=imin:igap:imax
                for J=jmin:jgap:jmax
                    tmp2=squeeze(QQ(I,J,:));                            % get the time series for 2nd point
                    if isnan(tmp2(1))==1
                        continue
                    end
                    a = (xk(I) - xk(i));
                    dx = a*cosd((y(J) + y(j))*.5);                      % scaled distance in X 
                    dy =  yk(J) - yk(j);
                    covmatrix=cov(tmp1,tmp2);
                    CC(I,J)=round(covmatrix(2,1)*exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))),3,'significant'); %covariance of point ij with point IJ
                end
            end
        end
    end
    float_locations(iter,:)=max_ij;
    float_I=[float_I round((max_ij(1)-imin)/igap)];                     % getting the float location indices in compressed QQ
    float_J=[float_J round((max_ij(2)-jmin)/jgap)];
    ND=length(float_I);
    % run each time for formal map error
    % need signal just at every data point
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
                dx = a*cosd((y(J) + y(j))*.5);                          % scaled distance in X
                dy = yk(J) - yk(j);
                DIST(nx,ny,num)=round(exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))),3,'significant'); %gaussian distances
            end
        end
    end

    if max_ij(1)==0
        % no max_ij means we found nothing to report
        continue
    else 
        DIST=reshape(DIST,NX*NY,ND);                                    % turn distances into a linear vector for each float
        MASK = MASK(:);                                                 % turned into a linear vector rather than matrix
        [Is] = find(MASK==1);                                           % where the floats are
        D = Qssv(Is,:);                                                 % getting just the float points
        float_ind = sub2ind(size(Qss(:,:,1)),float_I,float_J);
        DDist =DIST(float_ind,:);                                       % use float_ind instead of IS to make symmetric dist
        Cmd = Qssv*(D').*DIST;                                          % number of NXNY points that are ocean x number of floats
        %NOW GET DATA DATA COVARIANCE
        Cmd = Cmd(IsOcean,:);
        Cdd = D*(D').*DDist;                                            % get covariance between float locations
        %ADD NOISE TO DATA BECAUSE FLOATS AREN’T PERFECT
        NoiseValue = 1e-5;                                              % noise variance in obs keep in mind order of magnitude of pCO2 values is around 10^-4
        %NoiseValue = 0;
        Cdd = Cdd + NoiseValue*eye(ND);                                 % adding noise to the variance of float locations
        %FORMAL MAPPING ERROR
        FME = diag(Cmm) - diag(Cmd*(Cdd)^-1*(Cmd'));
        %NORMALIZED VERSION
        NME = FME./diag(Cmm);
        MapErr = zeros(NX*NY,1); % keeping it in the linear form
        MapErr(MapErr==0)=nan;
        MapErr(IsOcean) = NME; MapErr = reshape(MapErr,NX,NY);
        coverage_percent=100-100*sum(NME,'all')/length(IsOcean);
        % old stuff
        gamma=all_covs./all_covs(max_ij(1),max_ij(2)); 			  %normalize covariances by variance of float location
        QQtmp=round(QQ-gamma.*QQ(max_ij(1),max_ij(2),:),3,'significant'); %remove influence of float
        %float_QQ(:,:,:,iter)=QQtmp;
        QQ=QQtmp; %use the new pCO2 without the influence of float
        count=0; %reinitialize count for next run
        max_cov=1e-9; %reinitialize maximum covariance for next run
        PlotMapErr(:,:,iter)=MapErr;
        PlotCoverage(iter)=coverage_percent;
    end
end
toc

% Running the floats code
% in order to get the mapping data for the floats
cd /home/winnie.chu/functions
tic
[PlotMapErr,PlotCoverage,float_locations,Irange,Jrange]=find_floats(100);
toc

%Covariance plotting
function [PlotMapErr,PlotCoverage,float_locations,Irange,Jrange]=find_floats(num_iter)
