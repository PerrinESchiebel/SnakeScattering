function [VarList,SnakeForces] = SnakeForceHistogram_syncRuns_052617(anglemax,kappamin)
%%%%%%%%%%%MAKE A STRUCT WITH ALL RUN NAMES AND MEASURED THETA AND KAPPA%%%%%%%%%%%%%
directoryforcehist = 'F:\Dropbox\SnakeScattering\Figures\snake_codes\snake_mats_for_combinsResults\';
listforcehist = dir([directoryforcehist,'*.mat']);
directoryforcelocations = 'F:\Dropbox\SnakeScattering\Figures\snake_codes\allForces\';

nforceruns = length(listforcehist);
centertocenter = 19;   %%px for shag
% initialbox = [-centertocenter/2 centertocenter/2];
VarList = struct('name',[],'angle',[],'kappa',[],'theta',[],'pegpair',[],'TimePeg1',[],'TimePeg2',[],'tstart',[],'tend',[]);
virtualR = 6.5;   %%%%%From ImageJ measurements: peg diameter ~ 5.5 px. snake diameter ~ 7.5px, so (pegD + snakeD)/2 = 6.5
% pegR = 5.5/2;
%%
anglemax = 15;kappamin = 30;
anglemax = 15;kappamin = 15;
for zz=2:nforceruns
    load([directoryforcehist,listforcehist(zz).name]);
    if abs(atand(angles)) < anglemax && kappas > kappamin
        if strcmp(name(1),'1') == 0
            name = ['1',name];
        end
        VarList(zz).name = name;
        VarList(zz).angle = atand(angles);
        VarList(zz).kappa = kappas;
        pegXY = sortrows(pegXY,2);
        
%         heady = mean(mean(splineY(9:11,pegxing-2:pegxing+2)));
%         [~,I] = sort(abs(pegXY(:,2)-heady));
%         pegpair = I(1:2);
        
        %         hold off;
        %         plot(splineX,splineY,'Color',[0.7 0.7 0.7]);hold on;
        %         viscircles(pegXY,repmat(4,length(pegx),1));
        %         plot(mean(pegXY(:,1)),heady,'h','MarkerFaceColor',[0.7 0.1 0.8],'MarkerSize',12);
        %         plot(pegXY(pegpair,1),pegXY(pegpair,2),'+','MarkerSize',12,'LineWidth',3,'Color',[0.1 0.8 0.8]);drawnow;
        %         response = input('good?','s');
        %
        
        splineX = splineX - mean((pegXY(:,1)));
        splineY = splineY - mean((pegXY(:,2)));
        pegx = pegXY(:,1) - mean((pegXY(:,1)));
        pegy = pegXY(:,2) - mean((pegXY(:,2)));

        %%%%LOOK FOR PAIR OF PEGS SNAKE GOES BETWEEN
         pegxing = find(splineX(10,:) > 0,1,'first');  %%%%%THIS IS TO FIND WHERE THE SNAKE IS DEFINITELY BETWEEN THE PEGS
        heady = mean(mean(splineY(9:11,pegxing-2:pegxing+2)),'omitnan');
        [~,I] = sort(abs(pegy-heady));
        pegpair = I(1:2);
        tstart = find(splineX(1,:) > 0 - virtualR,1,'first');
        tend = find(splineX(end,:) > 0 + virtualR,1,'first');
        %%%%LOOK FOR WHEN SNAKE MAY ACTUALLY TOUCH THESE PEGS
        %%%%PEG1 AND PEG2 REFER TO PEGPAIR(1) AND PEGPAIR(2), NOT
        %%%%NECESSARILY THE PEG NUMBERS
        pegs = [pegx,pegy];
        segdistancesPEG1 = (splineX-pegs(pegpair(1),1)).^2 + (splineY-pegs(pegpair(1),2)).^2;
        segdistancesPEG2 = (splineX-pegs(pegpair(2),1)).^2 + (splineY-pegs(pegpair(2),2)).^2;
        [~,indt1] = find(segdistancesPEG1 <= virtualR^2);
        [~,indt2] = find(segdistancesPEG2 <= virtualR^2);
        
            hold off;
                plot(splineX,splineY,'Color',[0.7 0.7 0.7]);hold on;
                viscircles([pegx,pegy],repmat(5.5,length(pegx),1),'EdgeColor','k');
                viscircles([pegx(pegpair),pegy(pegpair)],repmat(5.5,2,1),'EdgeColor',[0 0.7 0.7]);
                plot(splineX(:,tstart),splineY(:,tstart),'Color',[1 0 0.5],'LineWidth',2);          
                                plot(splineX(:,tend),splineY(:,tend),'Color',[0 1 0.5],'LineWidth',2);      
                drawnow;pause(0.5);
%                 colors = colormap(brewermap(tend-tstart+1,'PuBu'));
%                 for xind = tstart:tend
% %                     plot(x{xind}-mean(pegXY(:,1)),y{xind}-mean(pegXY(:,2)),'Color','k','LineWidth',2.5);
%                     plot(splineX(:,xind),splineY(:,xind),'Color',colors(xind-tstart+1,:),'LineWidth',2);                   
%                 end
%                 drawnow;
%                 plot(mean(pegXY(:,1)),heady,'h','MarkerFaceColor',[0.7 0.1 0.8],'MarkerSize',12);
%                 plot(pegXY(pegpair,1),pegXY(pegpair,2),'+','MarkerSize',12,'LineWidth',3,'Color',[0.1 0.8 0.8]);drawnow;
        %         response = input('good?','s');
        %
        
        
        %         meany = mean(mean(splineY(1:50,1:pegxing)));
        %         count2 = 1;
        %
        %         while (meany > initialbox(1) && meany < initialbox(2)) == 0
        %             if count2 > 15
        %                 display('not in box');
        %                 break
        %             else
        %                 if meany < initialbox(1)
        %                     splineY = splineY+centertocenter;
        %                 else
        %                     splineY = splineY-centertocenter;
        %                 end
        %                 meany = mean(mean(splineY(:,1:pegxing)));
        %                 count2 = count2+1;
        %             end
        %
        %         end
        %         if count2 < 16
        %
        %             [th,rho] = cart2pol(splineX(1:50,pegxing:end)./vot,splineY(1:50,pegxing:end)./vot);
        %             cirr = 6;
        %             deltar = 1;
        %             rs = [cirr,cirr+deltar];
        %             %             ind = find(rho >= rs(1) & rho < rs(2));
%                     VarList(zz).theta = mean(th(rho >= rs(1) & rho < rs(2)));
        VarList(zz).pegpair = pegpair;
        VarList(zz).TimePeg1 = indt1;
        VarList(zz).TimePeg2 = indt2;
        VarList(zz).tstart = tstart;
        VarList(zz).tend = tend;
        %             if response == 'n'
        %                 error('problem');
        %             end
        %         end
    end
    %     clear x y splineX splineY pegXY angles kappas ind th rho indt1 indt2
end
%%
%%%%%%%%%%%%%%%LOADING AND PARSING FORCES%%%%%%%%%%%%%%%%%%%%%%
% colors1 = [linspace(38,0,7); linspace(166,77,7); linspace(154,64,7)]'./255;
% colors2 = [linspace(126, 49,7); linspace(87,27,7); linspace(194,146,7)]'./255;
% maindir = 'G:\scattering\Shag\';
% maindir = 'I:\scattering\Shag\';
% directories = dir(maindir);
% directories = directories((vertcat(directories.isdir)));
% directories = {directories.name};
% directories(ismember(directories,{'.','..'})) = [];
% [num,txt,~] = xlsread('L:\scattering\Shag\vidSync.xlsx');
[num,txt,~] = xlsread('F:\Dropbox (GaTech)\Research\scattering\snake_codes_SHAG\vidSync.xlsx');

syncFileNames = txt(2:end,2);
syncFrameNum = num(:,[3,5]);  %%%%First number is for force data, second is for tracked data
clear num txt

% allangs = [];
% allForces = [];
% allperp = cell(1,1);
% allparr = cell(1,1);
% allperpsum = cell(1,1);
% allparrsum = cell(1,1);
% Fs = 1/200;
truefalse = 0;
% Fmag = cell(3,3);
% forceindices = cell(3,3);
% momperpall = nan(nforceruns,7);
% momparrall = nan(nforceruns,7);
% contactlength = nan(nforceruns,7);
% passdir = nan(nforceruns,7);
% thresh = 0.005;
% count = 1;
% scatterangle = nan(length(directories),50);
% collectperp = nan(1,1);
% collectparr = nan(1,1);
% collecttheta = nan(1,1);
% collectangles = nan(1,1);
% cc = 1;
% ntouch = zeros(nforceruns,55);
SnakeForces = struct('Peg1',[],'Peg2',[]);
% Peg1Forces = cell(1,nforceruns);
% Peg2Forces = cell(1,nforceruns);
NotTouchedForces = cell(1,nforceruns);
%%
count = 1;
for ii = 2:nforceruns
    %     %     for ii = 3
    %     foldername = directories{ii};
    %     if isempty(str2double(foldername))
    %     else
    %         maindir = 'F:\Dropbox\SnakeScattering\Snake_Scattering_matfiles\allForces\';
    %         folderpath = strcat(maindir,'\tracked_Params\forceData\');
    %         vidpath = strcat(maindir,foldername,'\vids\');
    
    %         list = dir([directoryforcelocations,'*.avi']);
%     anglessn = [];
%     Forcessn = [];
%     momperp = nan(length(list),7);
%     momparr = nan(length(list),7);
%     sumperp = nan(length(list),1);
%     sumparr = nan(length(list),1);
    %         for jj=1:length(list)
    name = listforcehist(ii).name;
    name = name(1:end-4);
    for tt = 1:length(VarList)
        if strcmp(VarList(tt).name,name) == 1
            truefalse = 1;
            break
        end
    end
    if strcmp(name,'132_2p3_shag_H_072116') == 1
        truefalse = 0;
    end
    %     if strcmp(name(1:3),'132') == 0
    %         truefalse = 0;
    %     end
    match = 0;
    if truefalse == 1
        if abs(VarList(tt).angle) < anglemax && VarList(tt).kappa > kappamin
%             display(name);
            for aa = 1:length(VarList)
                if strcmp(syncFileNames{aa},[name,'.avi']) == 1
                    display(aa)
                    match = 1;
                    break
                end
            end
            if match == 1
            
            frameAdjust = syncFrameNum(aa,2) - syncFrameNum(aa,1);
            if isfinite(frameAdjust) == 1
%                 display(pegpair);
%             pegpair = VarList(tt).pegpair;
%             indt1 = VarList(tt).TimePeg1 + 1 - frameAdjust;  %%%ADJUST FOR TAKING FRAME DIFFERENCES BY ADDING 1 TO THE FRAME INDEX
%             indt2 = VarList(tt).TimePeg2 + 1 - frameAdjust;
%             
%             indt1(indt1<1) = [];
%             indt2(indt2<1) = [];
            pegpair = VarList(tt).pegpair;
            tstart = VarList(tt).tstart + 1 - frameAdjust;  %%%ADJUST FOR TAKING FRAME DIFFERENCES BY ADDING 1 TO THE FRAME INDEX
            tend = VarList(tt).tend + 1 - frameAdjust;
            
%             indt1(indt1<1) = [];
%             indt2(indt2<1) = [];
            %%%%%%NOW THESE INDICES CORRESPOND TO THE FORCE DATA
            %                     scatterangle(ii,jj) = theta(tt);
            
            fnamex = [directoryforcelocations,name,'_Params_Matrix_ForceX.mat'];
            fnamey = [directoryforcelocations,name,'_Params_Matrix_ForceY.mat'];
            %                     fnamemag = [directoryforcelocations,name,'_Params_Matrix_ForceMags.mat'];
            if exist(fnamex,'file') == 0 || exist(fnamey,'file') == 0
                display(['no forces for ',list(ii).name]);
            else
                load(fnamex);
                load(fnamey);
%                 %                         load(fnamemag);
% %                 F = sqrt(Fx.^2+Fy.^2);
% %                 [npegs,ntimes] = size(Fx);
% %                 clear discontinuityLocs sequencelengths ind longsequencelocsstart longsequencelocsend
% %                 sequencethresh = 7;
% %                 k = strfind(name,'16');
% %                 date = name(k-4:k+1);
%                 %                 if strcmp(date,'081616') || strcmp(date,'082216') || strcmp(date,'090116')
%                 %                     pegpair = VarList(tt).pegpair;
%                 %                 else
%                 %                     pegpair = sort(npegs+1-VarList(tt).pegpair);
%                 %                 end
% %                 hold off;
%                 indt1 = unique(indt1);
%                 indt2 = unique(indt2);
%                 indt1(indt1>=length(Fx)) = [];
%                 indt2(indt2>=length(Fx)) = [];
%                                 Fz1 = Fy(pegpair(1),indt1);
%                 Fx1 = -Fx(pegpair(1),indt1);
%                 Fz2 = Fy(pegpair(2),indt2);
%                 Fx2 = -Fx(pegpair(2),indt2);
timind = tstart:tend;
timind(timind<1) = [];
timind(timind>length(Fx)) = [];
 Fz1 = Fy(pegpair(1),timind);
                Fx1 = -Fx(pegpair(1),timind);
                Fz2 = Fy(pegpair(2),timind);
                Fx2 = -Fx(pegpair(2),timind);

                peglist = 1:size(Fy,1)';
                peglist(pegpair) = [];
                 randind = randi([1 length(peglist)],1);
                Fxn = -Fx(peglist(randind),:);
                Fzn = Fy(peglist(randind),:);

% str = input('Ringing?','s');
% switch str
%     case 'y'
%         hold off
%         plot(1:length(timind),Fx1,1:length(timind),Fz1);hold on;
%         plot(1:length(timind),Fx2,1:length(timind),Fz2);drawnow;
%         [xind,~] = ginput(2);
%         timind(xind(1):xind(2)) = [];
%         Fz1 = Fy(pegpair(1),timind);
%         Fx1 = -Fx(pegpair(1),timind);
%         Fz2 = Fy(pegpair(2),timind);
%         Fx2 = -Fx(pegpair(2),timind);
%         hold off;plot(timind,Fx1,timind,Fz1);
%                 plot(timind,Fx2,timind,Fz2);drawnow;
%                 str = input('better?','s');
%                 if strcmp(str,'n') == 1
%                     error('problem');
%                 end
% end
                SnakeForces(count).Peg1 = [Fx1;Fz1];
                SnakeForces(count).Peg2 = [Fx2;Fz2];
                NotTouchedForces{count} = [Fxn;Fzn];
%                 Peg1Forces{count} = [Fx1(VarList(ii).TimePeg1);Fz1(VarList(ii).TimePeg1)];
%                 Peg2Forces{count} = [Fx2(VarList(ii).TimePeg2);Fz2(VarList(ii).TimePeg2)];
%  hold off;
% plot(Fx','k');hold on;plot(Fy','k');
% plot(timind,Fx1,timind,Fz1);hold on;
% plot(timind,Fx2,timind,Fz2);
% plot(Peg1Forces{count}(1,:));plot(Peg2Forces{count}(1,:));drawnow;
%                    clear indt1 indt2
            end
            end
                %                             plot(-Fx(ll,:),'Color',colors1(ll,:),'LineWidth',4);
                %                             hold on;plot(Fy(ll,:),'Color',colors2(ll,:),'LineWidth',4);
                %                     ind = find(abs(F(ll,:))>thresh);
                %%%%THIS IS GETTING RID OF RINGING
%                     discontinuityLocs = find(diff(ind)>1);
%                     if isempty(discontinuityLocs) == 0
%                         discontinuityLocs = [1,discontinuityLocs,length(ind)];
%                         sequencelengths = diff(discontinuityLocs);
%                         longsequenceinds = find(sequencelengths>sequencethresh);
%                         longsequencelocsstart = discontinuityLocs(longsequenceinds);
%                         longsequencelocsend = discontinuityLocs(longsequenceinds+1);
%                         nseqs = length(longsequencelocsstart);
%                         runs = cell(nseqs,1);
%                         for mm=1:nseqs
%                             runs{mm} = ind(longsequencelocsstart(mm):longsequencelocsend(mm));
%                             Fpmx = sign(Fx(ll,runs{mm}));
%                             Fpmy = sign(Fy(ll,runs{mm}));
%                             if sum(abs(diff(Fpmx))>0) > 5 && sum(abs(diff(Fpmy))>0)>5
%                                 runs{mm} = [];
%                             end
%                         end
%                         ind = [runs{:}];
%                     end
                    %%%%%%%%
%                     if length(ind) < ntimes
%                         z = Fy(ll,ind);
%                         x = -Fx(ll,ind);
%                         mag = F(ll,ind);
%                         if isempty(x) == 0
%                             angleeach{count,1} = tt;
%                             angleeach{count,ll+1} = atan2(x,z);
%                             anglessn = [anglessn,atan2(x,z)];
%                             Forcessn = [Forcessn,mag];
%                             %                             momperp(jj,ll) = sum(Fs.*x);
%                             %                             momparr(jj,ll) = sum(Fs.*z);
%                             Fmag{ii,ll} = F(ll,ind);
%                             forceindices{ii,ll} = ind;
%                             momperpall(ii,ll) = sum(Fs.*x);
%                             momparrall(ii,ll) = sum(Fs.*z);
%                             contactlength(ii,ll) = length(ind);
%                             %                             if ll == pegpair(1)
%                             %                                 pass = 1; %%%%peg on left
%                             %                             elseif ll == pegpair(2)
%                             %                                 pass = 2; %%%%peg on right
%                             %                             elseif ll+1 == pegpair(1)
%                             %                                 pass = 3; %%%%%touched and skipped left hand pair
%                             %                             elseif ll-1 == pegpair(2)
%                             %                                 pass = 4; %%%%%touched and skipped to right hand pair
%                             %                             else
%                             %                                 pass = 5;
%                             %                                 problemnames{cc} = name;
%                             %                                 cc = cc+1;
%                             %                             end
%                             %                             passdir(ii,ll) = pass;
%                             %                             ntouch(ii) = ntouch(ii) + 1;
%                             %                                     plot(ind,x,'o','Color',[0 191 165]./255,'LineWidth',2);
%                             %                                     hold on;plot(ind,z,'o','Color',[140 158 255]./255,'LineWidth',2);
%                         end
%                     end
%                 end
                %                 if sum(isfinite(momperp)) == 1
                %                     pegnum = find(isfinite(momperp(jj,:)));
                %                     singlepass(ii,jj) = momperp(jj,pegnum);
                %                     singledir(ii,jj) = passdir(ii,jj,pegnum);
                %                 end
                
                %                         response = input('good?','s');
                %                         if response == 'n'
                %                             error('problem');
                %                         end
                %                         drawnow;
                
%                                 sumperp = sum(momperp,'omitnan');
%                                 sumparr = sum(momparr,'omitnan');
%                                 collectperp(ii) = sumperp;
%                                 collectparr(ii) = sumparr;
%                                 collecttheta(ii) = VarList(tt).theta;
%                                 collectangles(count) = mean(angles,'omitnan');
                count = count+1;
                hold off;
            end
%             allangs = [allangs,anglessn];
%             allForces = [allForces,Forcessn];
        end
%         truefalse = 0;
    end
    truefalse = 0;
    %     end
end

%     allperp{ii} = momperp;
%     allparr{ii} = momparr;
% allperpsum{ii} = sumperp;
% allparrsum{ii} = sumparr;



% close

%%%%%%%%%%%%%%%%PLOTS%%%%%%%%%%%%%%%
% % histogram(allangs)
% bins = linspace(-180,180,100);
% histogram(rad2deg(allangs),bins,'Normalization','pdf','FaceColor',[0.4 0.5 0.1],'EdgeColor','none')
% % thetalabels{1} = ['-\pi'];
% % thetalabels{2} = ['-\pi','/2'];
% % thetalabels{3} = ['0'];
% % thetalabels{4} = ['\pi','/2'];
% % thetalabels{5} = ['\pi'];
% set(gca,'XTick',[-180 -90 0 90 180],'XLim',[-180 180]);
% set(gca,'YTick',[0 0.7]);
% xlabel('\theta_{Force} [o]');
% ylabel('P(\theta)');
% set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% %%
% parr = cell2mat(allparr');
% parr = parr(isfinite(parr));
% histogram(parr,100,'Normalization','pdf','FaceColor',[0.5 0.3 0.7])
% xlabel('\DeltaP_{||} (Ns)')
% ylabel('Probability Density')
% set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% set(gca,'XTick',[-0.04 0 0.04],'YTick',[0 200])
% figure,
% parr = cell2mat(allparrsum');
% histogram(parr,100,'Normalization','pdf','FaceColor',[0.5 0.3 0.7])
% xlabel('\DeltaP_{||} (Ns)')
% ylabel('Probability Density')
% set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% set(gca,'XTick',[-0.04 0 0.04],'YTick',[0 200])
% %%
% parr = cell2mat(allperp');
% parr = parr(isfinite(parr));
% histogram(parr,100,'Normalization','pdf','FaceColor',[0.2 0.4 0.6])
% xlabel('\DeltaP_{\perp,pegs} (Ns)')
% ylabel('Probability Density')
% set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% set(gca,'XTick',[-0.04 0 0.04],'YTick',[0 100])
% figure,
% parr = cell2mat(allperpsum');
% histogram(parr,100,'Normalization','pdf','FaceColor',[0.2 0.4 0.6])
% xlabel('\DeltaP_{\perp} (Ns)')
% ylabel('Probability Density')
% set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% set(gca,'XTick',[-0.04 0 0.04],'YTick',[0 100])
% %%
% figure
% plot(reshape(momperpall,1,numel(momperpall)),...
%     reshape(contactlength,1,numel(momperpall)).*Fs,'o','MarkerFaceColor',[0.4 0.3 0.3],'MarkerEdgeColor','k','MarkerSize',12)
% set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% set(gca,'XTick',[-0.04 0 0.04],'YTick',[0 1 2 3])
% ylabel('Time [s]');xlabel('\DeltaP_{\perp} [Ns]');
% figure
% plot(abs(reshape(momperpall,1,numel(momperpall))),...
%     reshape(contactlength,1,numel(momperpall)).*Fs,'o','MarkerFaceColor',[0.4 0.3 0.3],'MarkerEdgeColor','k','MarkerSize',12)
% set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% set(gca,'XTick',[-0.04 0 0.04],'YTick',[0 1 2 3])
% ylabel('Contact Time [s]');xlabel('abs(\DeltaP_{\perp}) [Ns]');
% % figure
% % plot(abs(reshape(momperpall,1,numel(momperpall)))./reshape(contactlength,1,numel(momperpall)),...
% %     reshape(contactlength,1,numel(momperpall)).*Fs,'o','MarkerFaceColor',[0.4 0.3 0.3],'MarkerEdgeColor','k','MarkerSize',12)
% % set(gca,'linewidth',4,'fontsize',36,'fontname','helvetica','color',[1 1 1])
% % set(gca,'XTick',[-0.04 0 0.04],'YTick',[0 1 2 3])
% % ylabel('Contact Time [s]');xlabel('abs(\DeltaP_{\perp}) [Ns]');
% %%
% % collectperp = reshape(collectperp,numel(collectperp),1);
% % collecttheta = reshape(collecttheta,numel(collecttheta),1);
% % ind = find(collectperp == 0);
% % collectperp(ind) = [];
% % collecttheta(ind) = [];
% plot(collecttheta,collectperp,'o')
% hold on;
% for ii=1:23
%     for jj=1:45
%         if singlepass(ii,jj) == 0
%         else
%             plot(singlepass(ii,jj),collecttheta(ii,jj),'or');hold on;
%         end
%     end
% end