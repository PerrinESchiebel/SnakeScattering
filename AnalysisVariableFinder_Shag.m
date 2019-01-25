%%%Started 02/06/2018 by Perrin
%%%grab all poppy seed 2in peg trials and calculate variables.
%
addpath 'F:\Dropbox\Research\Chionactis_tracking\functions'
directoryresults2 = 'F:\Dropbox (GaTech)\Research\scattering\Data\snake_mats_for_combinsResults\';
listresults2 = dir([directoryresults2,'*.mat']);
%
numruns = size(listresults2,1);
% M = nan(numruns,15);
cutoff = 5;    %%changes how much is cut off of ends of run when finding mean
%%larger number = more cut off
Fs = 200;
anglemax = 15;   %%%%%%%MAXIMUM ANGLE OF TRAJECTORY BEFORE PEGS
kappamin = 30;
% centertocenter = 110;  %%px
% centertocenter = 8;
% initialbox = [-centertocenter/2 centertocenter/2];
% yesno = 0;
colors = [[0.9 0.9 1];[0 0.1 1];[0 0 0];[1 0.1 0];[1 0.9 0.9]];
colors = interp1(1:5,colors,linspace(1,5,64));
count = 1;
numsplinepts = 200;
% pxperm = 1364;
pxperm = 845; %%measured 121817 in imagej from
% {external 1}\scattering\Shag\RAWDATA\prelimData\071116\cal\cal_bed_071116.avi
lengths = zeros(132,1);
lengths(120) = 36.4;
lengths(122) = 39.2;
lengths(123) = 40.1;
lengths(124) = 38.1;
lengths(125) = 36.6;
lengths(128) = 38.0;
lengths(129) = 39.3;
lengths(130) = 37.1;
lengths(131) = 33.1;
lengths(132) = 38.9;
M = struct('snakenum',nan(1,numruns),'pegxing',nan(1,numruns),'vcom',nan(numsplinepts,numruns),'T',nan(numsplinepts,numruns),'time',{{}},'vseg',{{}},'lambdas',{{}},...
    'lambda',{{}},'KappaMax',{{}},'dist',{{}},'kpld',{{}},'ksi',{{}},'thetamax',{{}});
%%
close all
for wz=1:numruns
    clear tempx tempy xpos ypos Curvature x y curvemaxstuff pointd;
    fnamealanlysis = listresults2(wz).name;
    display(fnamealanlysis);
    load(strcat(directoryresults2,'\',fnamealanlysis));
    if abs(atand(angles)) < anglemax && abs(kappas) > kappamin
        xs = nan(numsplinepts,length(xSort));
        ys = nan(numsplinepts,length(xSort));
        % ntimes = length(xSort);
        ntimes = length(xSort);
        for jj=1:ntimes
            %     x = xSort{jj};y = ySort{jj};
            x = xSort{jj};y = ySort{jj};
            if isempty(x) == 0 && length(x)>50
                if x(1)==x(end)
                    x(end) = x(end)+1;
                elseif y(1)==y(end)
                    y(end) = y(end)+1;
                end
                [psx,psy,~] =  cubesplineinterp(x,y,numsplinepts);
                xs(:,jj) = psx;
                ys(:,jj) = psy;
            end
        end
        xss = nan(numsplinepts,ntimes);
        yss = nan(numsplinepts,ntimes);
        for jj=1:numsplinepts
            xss(jj,:) = smooth(xs(jj,:),7);
            yss(jj,:) = smooth(ys(jj,:),7);
            %              plot(xs(jj,:),ys(jj,:),'LineWidth',2);hold on;
            %              plot(xss(jj,:),yss(jj,:),'LineWidth',4);drawnow;pause(0.3);hold off;
        end
        if xss(end,round(ntimes/2))>xss(1,round(ntimes/2))
            x = flipud(xss)./pxperm;
            y = flipud(yss)./pxperm;
        else
            x = xss./pxperm;
            y = yss./pxperm;
        end
        [numpts,numframes] = size(x);
        
        anginc=10;
        [Curvature] = spatialCurvature(x,y,anginc);
        Curvature(1:anginc,:) = [];
%         Curvature = Curvature(anginc+1:end,:);
        cmax = max(abs(Curvature));
        
        [TemporCurvature] = temporalCurvature(x,y,2);
%         TemporCurvature = TemporCurvature(:,isfinite(TemporCurvature(1,:)));
        seginc = 3;
        theta = tangentangle3(x,y,seginc,1);   %%%last 1 means this will subtract the mean theta to correct for snake moving on an angle
        [frequency,vcom] = FrequencyVelocity(x,y,TemporCurvature,Fs,cutoff);
        vseg = sqrt(diff(x,[],2).^2 + diff(y,[],2).^2).*Fs;
        PeriodinFrames = round(Fs/mean(frequency,'omitnan'));
%         [lambdas, lambda,kappamaxfrompeaks,kplds] = KappaLambdas(x,y,Curvature);
        distance = DistancePerUndulation(x,y,PeriodinFrames,TemporCurvature);
        
        pegxing = find(x(1,:) < (mean(pegXY(:,1),'omitnan')-mean(pegR))./pxperm,1,'last');
        
        if strcmp(fnamealanlysis(1),'1')==0
            snnum = str2double(fnamealanlysis(1:2));
            snnum = 100+snnum;
        else
            snnum = str2double(fnamealanlysis(1:3));
        end
        
        
        
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
        t = -pegxing:length(vseg)-pegxing-1;
%         errorbar(t,mean(vseg,'omitnan'),std(vseg,'omitnan'),'o');title([num2str(mean(mean(vseg,'omitnan'))),'\pm',num2str(std(mean(vseg,'omitnan')))]);
%         drawnow;

        M.snakenum(count) = snnum;    %%Snake number
        M.pegxing(count) = pegxing;
        M.vcom(:,count) = vcom;  %%%average com speed
        M.T(:,count) = 1./frequency;   %%%%Temporal period
        M.time{count} = t;
        M.vseg{count} = vseg;
        M.lambdas{count} = {lambdas};  %%%% found average arclength, exlude values outside one std
        M.lambda{count} = {lambda};          %%%% found average wavelength, exlude values outside one std
        M.KappaMax{count} = {kappamaxfrompeaks};
        M.dist{count} = {distance};   %%%distance travelled per undulation--measured as euclidean distance between same point going through zero curvature
        M.kpld{count} = {kplds};
        M.ksi{count} ={ lengths(snnum)./lambdas};  %%%ksi
        M.thetamax{count} = {max(rad2deg(theta))};
        count = count+1;
    end
end
% writeVideo(writerObj,Movie);
% close(writerObj);