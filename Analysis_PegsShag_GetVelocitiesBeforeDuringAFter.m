%%%Started 02/06/2018 by Perrin
%%%grab all poppy seed 2in peg trials and calculate variables.

%  load('F:\Dropbox (GaTech)\Research\scattering\Data\WKSP_shagPegs_cut128_091218.mat')

% load('F:\Dropbox (GaTech)\Research\scattering\Data\WKSP_shagPegs_no128_angle15Kappamax15_233Runs.mat')
addpath 'F:\Dropbox\Research\Chionactis_tracking\functions'

% numruns = length(xtail);
numruns = length(allX);

cutoff = 5;    %%changes how much is cut off of ends of run when finding mean
%%larger number = more cut off
Fs = 200;
dt = 1/Fs;

count = 1;
numsplinepts = 200;
% pxperm = 1364;
pxperm = 845; %%measured 121817 in imagej from
% {external 1}\scattering\Shag\RAWDATA\prelimData\071116\cal\cal_bed_071116.avi
lengths = zeros(132,1);
lengths(120) = 36.4;
lengths(121) = 36.5;
lengths(122) = 39.2;
lengths(123) = 40.1;
lengths(124) = 38.1;
lengths(125) = 36.6;
lengths(128) = 38.0;
lengths(129) = 39.3;
lengths(130) = 37.1;
lengths(131) = 33.1;
lengths(132) = 38.9;
% M = struct('snakenum',nan(1,numruns),'pegxing',nan(1,numruns),'vcom',nan(numsplinepts,numruns),'T',nan(numsplinepts,numruns),'time',{{}},'vseg',{{}},'lambdas',{{}},...
%     'lambda',{{}},'KappaMax',{{}},'dist',{{}},'kpld',{{}},'ksi',{{}},'thetamax',{{}},'vcomTimeResolved',{{}});
Vbefore = nan(1,numruns);Vduring = nan(1,numruns);Vafter = nan(1,numruns);
VSEGbefore = nan(1,numruns);VSEGduring = nan(1,numruns);VSEGafter = nan(1,numruns);
Frequency = nan(1,numruns);TimeDuring = nan(1,numruns);
sn(sn==128) = [];
close all
for wz=1:numruns
% for wz  =45
    clear tempx tempy xpos ypos Curvature x y curvemaxstuff pointd;
    fprintf(1,'Run Number %i \n',wz);

%         ntimes = length(xSort);

% xss = xtail{wz}';yss = ztail{wz}';
xss = allX{wz};
 yss = allY{wz};
 if isempty(xss)==0 %|| ii==48 %bad tracking at beginning of 48
[numpts,ntimes] = size(xss);
        if xss(end,round(ntimes/2))>xss(1,round(ntimes/2))
            x = flipud(xss)./pxperm;
            y = flipud(yss)./pxperm;
        else
            x = xss./pxperm;
            y = yss./pxperm;
        end
        
        
%         anginc=5;
%         [Curvature] = spatialCurvature(x,y,anginc);
%         Curvature(1:anginc,:) = [];
% %         Curvature = Curvature(anginc+1:end,:);
%         cmax = max(abs(Curvature));
%         
        [TemporCurvature] = temporalCurvature(x,y,2);
% %         TemporCurvature = TemporCurvature(:,isfinite(TemporCurvature(1,:)));
%         seginc = 3;
%         theta = tangentangle3(x,y,seginc,1);   %%%last 1 means this will subtract the mean theta to correct for snake moving on an angle
        [frequency,vcom] = FrequencyVelocity(x,y,TemporCurvature,Fs,cutoff);
%         vseg = sqrt(diff(x,[],2).^2 + diff(y,[],2).^2).*Fs;
%         PeriodinFrames = round(Fs/mean(frequency,'omitnan'));
% %         [lambdas, lambda,kappamaxfrompeaks,kplds] = KappaLambdas(x,y,Curvature);
%         distance = DistancePerUndulation(x,y,PeriodinFrames,TemporCurvature);
        
%         pegxing = find(x(1,:) < (mean(pegXY(:,1),'omitnan')-mean(pegR))./pxperm,1,'last');
        
%         if strcmp(fnamealanlysis(1),'1')==0
%             snnum = str2double(fnamealanlysis(1:2));
%             snnum = 100+snnum;
%         else
%             snnum = str2double(fnamealanlysis(1:3));
%         end
        velSEG =  mean(sqrt(diff(x,[],2).^2+diff(y,[],2).^2)./dt,'omitnan');
        vel = sqrt(diff(mean(x,'omitnan')).^2+diff(mean(y,'omitnan')).^2)./dt;
%         pegxing = find(x(1,:)>0,1,'first');
pegxing = find(x(end,:)>0,1,'first');
        pegxinglast = find(x(end,:)>0,1,'first');
        if pegxing < 50 || (length(vel)-pegxinglast)<50
            fprintf(1,'Not enough data run %i',wz);
        else
        vbefore = vel(1:pegxing);
        vduring = vel(pegxing:pegxinglast);
        vafter = vel(pegxinglast:end);
       vSEGbefore = velSEG(1:pegxing);
        vSEGduring = velSEG(pegxing:pegxinglast);
        vSEGafter = velSEG(pegxinglast:end);
%         subplot(4,4,[1,2,5,6]);
%         hold off;
%         plot(cell2mat(xSort')./pxperm,cell2mat(ySort')./pxperm,'ok','MarkerFaceColor','k','MarkerEdgeColor','none');hold on;
%         plot(x,y,'Color',[0.6 0 0.6 0.2],'LineWidth',2);
%         viscircles(pegXY./pxperm,ones(1,size(pegXY,1)).*0.0064,'EdgeColor','k');
%         plot(x(:,pegxing),y(:,pegxing),'Color','c');
%         axis equal tight
%         subplot(4,4,[3,4,7,8]);
%         pcolor(Curvature);shading flat;colorbar;colormap(colors);
%         subplot(4,4,[9,10]);
%         t = -pegxing:length(kplds)-pegxing-1;
%         plot(t,kplds,'o','MarkerFaceColor',[0.2 0.7 0.7]);title([num2str(mean(kplds,'omitnan')),'\pm',num2str(std(kplds,'omitnan'))]);
%         subplot(4,4,[13,14]);
%         plot(t,max(abs(theta)),'o','MarkerFaceColor',[1 0.5 1]);
%          subplot(4,4,[11,12]);
%         plot(abs(vcom),'LineWidth',3);title([num2str(mean(vcom,'omitnan')),'\pm',num2str(std(vcom,'omitnan'))]);
%         subplot(4,4,[15,16]);
%         t = -pegxing:length(vseg)-pegxing-1;
%         errorbar(t,mean(vseg,'omitnan'),std(vseg,'omitnan'),'o');title([num2str(mean(mean(vseg,'omitnan'))),'\pm',num2str(std(mean(vseg,'omitnan')))]);
%         drawnow;

%         M.snakenum(count) = sn(wz);    %%Snake number
%         M.pegxing(count) = pegxing;
%         M.vcom(:,count) = vcom;  %%%average com speed
%         M.vcomTimeResolved{count} = vel;
%         M.T(:,count) = 1./frequency;   %%%%Temporal period
%         M.time{count} = t;
%         M.vseg{count} = vseg;
%         M.lambdas{count} = {lambdas};  %%%% found average arclength, exlude values outside one std
%         M.lambda{count} = {lambda};          %%%% found average wavelength, exlude values outside one std
%         M.KappaMax{count} = {kappamaxfrompeaks};
%         M.dist{count} = {distance};   %%%distance travelled per undulation--measured as euclidean distance between same point going through zero curvature
%         M.kpld{count} = {kplds};
%         M.ksi{count} ={ lengths(snnum)./lambdas};  %%%ksi
%         M.thetamax{count} = {max(rad2deg(theta))};
        Vbefore(wz) = mean(vbefore,'omitnan');
        Vduring(wz) = mean(vduring,'omitnan');
        Vafter(wz) = mean(vafter,'omitnan');
        VSEGbefore(wz) = mean(vSEGbefore,'omitnan');
        VSEGduring(wz) = mean(vSEGduring,'omitnan');
        VSEGafter(wz) = mean(vSEGafter,'omitnan');
        TimeDuring(wz) = VSEGduring(wz)./lengths(sn(wz));
        Frequency(wz) = mean(frequency,'omitnan');
        count = count+1;
%         subplot(3,1,1);
%         plot(wz,Vbefore(wz),'o','MarkerFaceColor',[0.8 0 0.8],'MarkerSize',10);hold on;
%         plot(wz,Vduring(wz),'o','MarkerFaceColor',[0 0.8 0.8],'MarkerSize',10);plot(wz,Vafter(wz),'o','MarkerFaceColor',[0.8 0.8 0],'MarkerSize',10);
%        plot(wz,VSEGbefore(wz),'*','Color',[0.8 0 0.8],'MarkerSize',10);hold on;
%         plot(wz,VSEGduring(wz),'*','Color',[0 0.8 0.8],'MarkerSize',10);plot(wz,VSEGafter(wz),'*','Color',[0.8 0.8 0],'MarkerSize',10);
%         subplot(2,1,2);
%         plot(wz,Frequency(wz),'o','MarkerFaceColor',[0.1 0.1 0.1]);hold on;
%         drawnow;
        end
    end
end
% p = WelchsStudents_ttest(mean(Vbefore),std(Vbefore),mean(Vafter),std(Vafter),numruns,numruns);
% fprintf(1,'Test vcom %f \n',p);
% p = WelchsStudents_ttest(mean(VSEGbefore),std(VSEGbefore),mean(VSEGafter),std(VSEGafter),numruns,numruns);
% fprintf(1,'Test vseg %f \n',p);
%%
figure
[vals,bins] = histcounts(Vbefore.*100,'Normalization','pdf');
plot(bins(1:end-1)-(bins(2)-bins(1))/2,vals,'LineWidth',3,'Color','k');hold on;
vals = histcounts(Vduring.*100,bins,'Normalization','pdf');
plot(bins(1:end-1)-(bins(2)-bins(1))/2,vals,'LineWidth',3,'Color','r');hold on;
vals = histcounts(Vafter.*100,bins,'Normalization','pdf');
plot(bins(1:end-1)-(bins(2)-bins(1))/2,vals,'LineWidth',3,'Color','b');hold on;


set(gca,'FontSize',24,'LineWidth',3);
set(gca,'XTick',0:10:100,'YTick',0:0.02:0.04);
xlabel('vcom (cm/s)');
ylabel('Probability Density');
legend('Vbefore','Vduring','Vafter');
set(gcf,'Position',[10 10 1500 800])
plotname = 'F:\Dropbox (GaTech)\Research\scattering\plots\Figures_supporting\VbeforeVduringVafter_pdfs_233runs_vcom';
saveas(gcf,[plotname,'.fig']);
saveas(gcf,[plotname,'.emf']);
saveas(gcf,[plotname,'.jpg']);
%%
load('F:\Dropbox (GaTech)\Research\scattering\snake_codes_SHAG\WKSP_110118_ForceAnalysis_angle15Kappa15');
Fs = 200;
figure
tforce = forcetouch./Fs;
tpass  = totaltouch./Fs;
% ind = find(tforce<10);
% tforce = tforce(ind);tpass = tpass(ind);

[vals,bins] = histcounts(tpass,'Normalization','pdf');
plot([bins(1),bins+(bins(2)-bins(1))/2],[0,vals,0],'LineWidth',3,'Color','k');hold on;

[vals,bins] = histcounts(tforce,'Normalization','pdf');
plot([bins(1),bins+(bins(2)-bins(1))/2],[0,vals,0],'LineWidth',3,'Color','r');hold on;
% [vals,bins] = histcounts(1./Frequency,'Normalization','pdf');
% L = 0.37;
% [vals,bins] = histcounts(VSEGduring./L,'Normalization','pdf');

 [vals,bins] = histcounts(0.7.*tpass,'Normalization','pdf');
plot([bins(1),bins+(bins(2)-bins(1))/2],[0,vals,0],'LineWidth',3,'Color','b');hold on;
legend('time snake is passing the array','time forces are above threshhold','transit time (vseg/L)');
set(gca,'FontSize',24,'LineWidth',3);
xlabel('time (s)');ylabel('Probability Density');
set(gca,'XTick',0:0.25:3,'XTickLabel',{'0','','','','1','','','','2','','','','3'},'YTick',0:0.3:1.2);
set(gcf,'Position',[10 10 1500 800])

annotation('textbox',[0.5 0.5, 0.5 0.5],'String','at 70% of the time for the tracked points to pass the array compared to time forces are above thresh 0.005 p=0.26 paired t-test','LineStyle','none','FontSize',18)

% plotname = 'F:\Dropbox (GaTech)\Research\scattering\plots\Figures_supporting\time forces measured and time transit_pdfs_233runs_vcom';
% saveas(gcf,[plotname,'.fig']);
% saveas(gcf,[plotname,'.emf']);
% saveas(gcf,[plotname,'.jpg']);
%  p = WelchsStudents_ttest(mean(totaltouch./Fs),std(totaltouch./Fs),mean(forcetouch./Fs),std(forcetouch./Fs),numruns,numruns);

[h,p,ci,stats] = ttest(0.7*tpass,tforce);
 fprintf(1,'t-test of force times %f \n',p)