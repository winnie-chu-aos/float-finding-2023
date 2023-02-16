% This function produces figures showing the uncertainty in pCO2, float locations in white, and cruise path if necessary
% Format is plot_float_coverage(PlotMapErr,PlotCoverage,float_locations,waypoints,Irange,Jrange,savemovie,boatpath)
% PlotMapErr, PlotCoverage, float_locations, waypoints, Irange, and Jrange should all come out of the float_boat function 
% savemovie option is a string, 'N' for no save or 'NameOfMovie' for save to an avi file in /home/winnie.chu/figures 
% boatpath option is a string, 'Y' for displaying cruise path or 'N' for no path

% LAST WORKING VERSION USE THIS TO PLOT THE OUTPUT OF FIND FLOATS
function plot_float_coverage(PlotMapErr,PlotCoverage,float_locations,waypoints,Irange,Jrange,savemovie,boatpath)

addpath /home/mmazloff/ANALYSIS  %this has imagescnan and rdmds
addpath /home/mmazloff/ANALYSIS/export_fig/ %this has export fig
cd /data/SO6/TPOSE/bgc_tpose6_k51/2004_fwd_run/diags  %old model
load  /data/SO6/TPOSE/bgc_tpose6_k51/grid/grid XC YC; 
x = XC(:,1); y = YC(1,:); clear XC YC

if savemovie=='N'
    movieon=0;
else
    cd /home/winnie.chu/figures % goes to the right directory to save file
    v = VideoWriter([savemovie '.avi'],'Uncompressed AVI');
    v.FrameRate = 1; %set to 3 frames per second
    open(v)
    movieon=1;
end

for iter=1:length(PlotCoverage)
    figure()
    title(sprintf('percent coverage of %d floats is %.2f%%',iter,PlotCoverage(iter)))
    tmpplot=PlotMapErr(:,:,iter);
    imagescnan(x(Irange),y(Jrange),tmpplot'); axis xy; a=colorbar;
    caxis([0,1])
    a.Label.String = 'normalized mapping error';
    xlabel('longitude');
    ylabel('latitude');
    set(gcf, 'Color', 'w')
    hold on
    text(x(float_locations(1:iter,1)),y(float_locations(1:iter,2)),'\otimes','Color','white','FontSize',16,'HorizontalAlignment','center','FontWeight','bold')
    if boatpath=='Y'
        plot([x(waypoints(1,1)), x(waypoints(2,1))], [y(waypoints(1,2)), y(waypoints(2,2))], '--ok', 'LineWidth', 1, 'MarkerFaceColor', 'k')
    end
    hold off
    if movieon==1
        drawnow
        frame= getframe(gcf);
        writeVideo(v,frame);
        close;
    end
end

if movieon==1
    close(v);
end
cd /home/winnie.chu/functions
end