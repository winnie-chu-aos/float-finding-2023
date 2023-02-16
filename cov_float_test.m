function [PlotMapErr,PlotCoverage,float_locations,Irange,Jrange]=cov_float_test(num_iter)

% testing the code
%a simple model of the covariance code
% prompt = 'Input number of floats deployed: ';
% num_iter = input(prompt);

addpath /home/mmazloff/ANALYSIS  %this has imagescnan and rdmds
addpath /home/mmazloff/ANALYSIS/export_fig/ %this has export fig
cd /data/SO6/TPOSE/bgc_tpose6/2004_fwd_run/diags %new model
load  /data/SO6/TPOSE/bgc_tpose6/grid/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC
cd /home/winnie.chu/data
load pCO2.mat pCO2res
%load DCall_1deg.mat DCL

% initialization
QQ=pCO2res;
max_cov=0;
%plot_var=zeros(size(QQ,1),size(QQ,2)); %initialization
%oneQQ=ones(size(QQ,1),size(QQ,2)); % ones matrix the same size as QQ
float_locations=zeros(num_iter,2);
%xk = x*deg2km(1);yk = y*deg2km(1); 
imin=276;
imax=996; %range of i points
igap=6;
jmin=96;
jmax=216; %range of j points
jgap=6;
Irange=imin:igap:imax;
Jrange=jmin:jgap:jmax;
CC=zeros(size(QQ,1),size(QQ,2));
float_I=[];
float_J=[];
%dcx=3500; % just general estimate of decorrelation lengths in x direction
%dcy=1000; % estimate of decorrelation lengths in y direction
% replace this later once you get all decorrelation lengths
% make sure that when you do correlation analysis it matches this one so
% that we can have a decorrelation length for each of the points  
PlotMapErr=zeros(length(imin:igap:imax),length(jmin:jgap:jmax),num_iter);
PlotCoverage=zeros(num_iter);

%MATT’S ATTEMPT AT Cmm
%SUBSAMPLE DATA
Qss = QQ(imin:igap:imax,jmin:jgap:jmax,:);
[NX,NY,NT] = size(Qss); % getting the dimensions of the compressed matrix
Qssv = reshape(Qss,NX*NY,NT); % turning the xy into one dimension vs. time
[IsOcean,~] = find(isnan(Qssv(:,1))==0); % finding all the indices that are not land
Qsssv = Qssv(IsOcean,:); % compressing further to only ocean points
%NM = length(IsOcean); % new length of xy that is on ocean
Cmm = Qsssv*(Qsssv'); % Cmm is the local variances size should be NMxNM

%TO VISUALIZE:
%Qvar = zeros(NX*NY,1); Qvar(IsOcean) = diag(Cmm); Qvar = reshape(Qvar,NX,NY);
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
            %dcx=DC(i,j,1);
            %dcy=DC(i,j,2);
            tmp1=squeeze(QQ(i,j,:)); %getting the time series of first point
            if isnan(tmp1(1))==1
                %count=0; 
                continue
            elseif count > 0
            %if count > 0
                % if not the first run, then find the total covariance and
                % compare to previous covariances get total magnitude of
                % covariance
                %cov_ij=sum(abs(CC),'all','omitnan');
                % CC is the matrix of the previous point
                if j>jmin
                    cov_ij=sum(CC.*CC,'all','omitnan')/var(QQ(i,j-jgap,:));
                else
                    cov_ij=sum(CC.*CC,'all','omitnan')/var(QQ(i-igap,jmax,:));
                end
                %test_ij=isequal([i j], oldmax_ij); %is this an existing float location
                %test_repeat=ismember([i,j],float_locations,'rows'); % check if float location is already found
                if cov_ij > max_cov %&& test_repeat==0
                    % if not an existing float location and bigger
                    % covariance
                    max_cov=cov_ij; %store the highest covariances
                    all_covs=CC; %store covariances over entire area
                    if j>jmin
                        max_ij=[i,j-jgap] %store float location
                    else
                        max_ij=[i-igap,jmax]
                    end
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
                    %a = (xk(I) - xk(i));
                    %dx = a*cosd((y(J) + y(j))*.5);% scaled distance in X 
                    %dy =  yk(J) - yk(j);
                    % (admittedly large) radius of influence
                    %if (dx.^2 + dy.^2) < 1000000 %max distance
                    %if dx <= 1500 && dy <= 350 && (dx.^2 + dy.^2) < 1000000
                        % remember to half dx and dy if you're using
                        % correlation lengths
                        covmatrix=cov(tmp1,tmp2);
                        %CC(I,J)=round((1/(dist+1)).*covmatrix(2,1),3,'significant');
                        CC(I,J)=covmatrix(2,1);%*exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))); %covariance of point ij with point IJ
                    %else
                        %CC(I,J)=0;
                    %end
                end
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
        float_locations(iter,:)=max_ij;
        float_I=[float_I round((max_ij(1)-imin)/igap)+1]; %getting the float location indices in compressed QQ
        float_J=[float_J round((max_ij(2)-jmin)/jgap)+1];
        ND=length(float_I);
        %NOW THIS GOES AFTER YOUR ITERATIVE LOOP. RUN EACH TIME TO GET “FORMAL MAP ERROR”
        %Need signal just at every data point
        MASK = zeros(NX,NY);
        DIST = zeros(NX,NY,ND);
        for num = 1:ND %NDATA
            MASK(float_I(num),float_J(num)) = 1;
        end
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
%                     %DIST(nx,ny,num)=exp(-(dx.^2/(2.*dcx.^2)+dy.^2/(2.*dcy.^2))); %gaussian distances
% 
%                 end
%             end
%             DIST=ones(NX,NY,ND);
%         end
        DIST=ones(NX,NY,ND);
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
        NoiseValue = 1; %Noise variance in obs keep in mind order of magnitude of pCO2 values is around 10^-4
        %NoiseValue = 0;
        Cdd = Cdd + rand*NoiseValue*eye(ND); % adding noise to the variance of float locations
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
        %gamma=all_covs./all_covs(max_ij(1),max_ij(2)); %normalize covariances by variance of float location
        gamma=all_covs./var(QQ(max_ij(1),max_ij(2),:));
        %QQtmp=round(QQ-gamma.*QQ(max_ij(1),max_ij(2),:),3,'significant'); %remove influence of float
        QQ=QQ-gamma.*QQ(max_ij(1),max_ij(2),:); %remove influence of float
        %float_QQ(:,:,:,iter)=QQtmp;
        %QQ=QQtmp; %use the new pCO2 without the influence of float
        %count=0; %reinitialize count for next run
        %max_cov=1e-9; %reinitialize maximum covariance for next run
        max_cov=0;
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
cd /home/winnie.chu/functions
end
