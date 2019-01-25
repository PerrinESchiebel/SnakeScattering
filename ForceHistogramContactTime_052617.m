% function [SnakeForces,RingInd] = SnakeForceHistogram_syncRuns_clickRinging(anglemax,kappamin,clickringing)
%%
%%%%%%%%%%%MAKE A STRUCT WITH ALL RUN NAMES AND MEASURED THETA AND KAPPA%%%%%%%%%%%%%
%%%%%USED anglemax = 15;kappamin = 30;
directoryforcehist = 'F:\Dropbox\SnakeScattering\Figures\snake_codes\snake_mats_for_combinsResults\';  %%%% TRACKED SNAKES AFTER SNAKE RUN SORTING
listforcehist = dir([directoryforcehist,'*.mat']);
directoryforcelocations = 'F:\Dropbox\SnakeScattering\Figures\snake_codes\allForces\';   %%%%FORCES FROM TRACKED PEGS
% [num,txt,~] = xlsread('L:\scattering\Shag\vidSync.xlsx');      %%%%BY HAND, CHOSE MATCHING FRAMES IN THE TWO VIDEOS
[num,txt,~] = xlsread('vidSync.xlsx');      %%%%BY HAND, CHOSE MATCHING FRAMES IN THE TWO VIDEOS
syncFileNames = txt(2:end,2);
syncFrameNum = num(:,[3,5]);  %%%%First number is for force data, second is for tracked data
clear num txt
nforceruns = length(listforcehist);
SnakeForces = struct('Peg1',[],'Peg2',[]);
virtualR = 6.5 + 10;   %%%%%From ImageJ measurements: peg diameter ~ 5.5 px. snake diameter ~ 7.5px, so (pegD + snakeD)/2 = 6.5
%%%%% 012317, I noticed sometimes
%%%%% the time is a little late
%%%%% because of the tracking not
%%%%% getting all the way to the
%%%%% head so I added 3 px
count = 1;
dmat = [{'081616'},{'082216'},{'090116'},{'082316'}];
RingInd = nan(100,2);
totaltouch = nan(100,1);
forcetouch = nan(100,1);
savestuff = nan(10000,100);
theta = nan(100,1);
PegTouch = nan(100,1);
% anglemax = 15;kappamin = 30;
anglemax = 15;kappamin = 15;
thresh = 0.003;

vot = 17.5; %%%avg for all snakes IN CM
vot = vot*5*0.25*2.54;   %%%%CONVERT TO PX
rs = [6*vot,7*vot];


%%%%  THE FIRST MATRIX IS BROKEN, START AT 2
for mq=1:nforceruns
    display(mq)
    load([directoryforcehist,listforcehist(mq).name]);
    truefalse = 0;
    if abs(atand(angles)) < anglemax && kappas > kappamin
        %%%%SOME FILES WERE SAVED WITHOUT THE '1' E.G. SNAKE 20 VS 120
        if strcmp(name(1),'1') == 0
            name = ['1',name];
        end
        if strcmp(name,'132_2p3_shag_H_072116') == 0
            for aa = 1:length(syncFileNames)
                if strcmp(syncFileNames{aa},[name,'.avi']) == 1
                    %                     display(aa)
                    truefalse = 1;
                    break
                end
            end
            if truefalse ==0;
                display([name,'    is missing sync values']);
            else
                pegXY = sortrows(pegXY,2);
                splineX = splineX - mean((pegXY(:,1)));
                splineY = splineY - mean((pegXY(:,2)));
                pegx = pegXY(:,1) - mean((pegXY(:,1)));
                pegy = pegXY(:,2) - mean((pegXY(:,2)));
                
                %%%%LOOK FOR PAIR OF PEGS SNAKE GOES BETWEEN
                pegxing = find(splineX(10,:) > 0,1,'first');  %%%%%THIS IS TO FIND WHERE THE SNAKE IS DEFINITELY BETWEEN THE PEGS
                heady = mean(mean(splineY(9:11,pegxing-2:pegxing+2)),'omitnan');
                [~,I] = sort(abs(pegy-heady));
                pegpair = I(1:2);
                 [th,rho] = cart2pol(splineX(:,pegxing:end),splineY(:,pegxing:end));
                  ind = find(rho >= rs(1) & rho < rs(2));
                  savestuff(1:numel(ind),count) = th(ind);
                %%%%THE PEG DATA FROM THE FIRST DAY TRACKED THE PEGS IN REVERSE
                %%%%ORDER, SO CHECK FOR THIS AND CORRECT
                
                %%%%FIND WHEN HEAD FIRST CROSSES LINE OF PEGS AND TAIL LAST
                %%%%CROSSES
                tstart = find(splineX(1,:) > 0 - virtualR,1,'first');
                tend = find(splineX(end,:) > 0 + virtualR,1,'first');
                if isempty(tstart) || isempty(tend)
                    display([name,'   not enough snake']);
                else
                    
%                     subplot(2,1,1);
%                     plot(splineX,splineY,'Color',[0.7 0.7 0.7]);hold on;
%                     viscircles([pegx,pegy],repmat(5.5,length(pegx),1),'EdgeColor','k');
%                     viscircles([pegx(pegpair),pegy(pegpair)],repmat(5.5,2,1),'EdgeColor',[0 0.7 0.7]);
%                     plot(splineX(:,tstart),splineY(:,tstart),'Color',[1 0 0.5],'LineWidth',2);
%                     plot(splineX(:,tend),splineY(:,tend),'Color',[0 1 0.5],'LineWidth',2);      drawnow;hold off;
                    
                    
                    
                    
                    frameAdjust = syncFrameNum(aa,2) - syncFrameNum(aa,1);  %%%%DIFFERENCE IN TIME BETWEEN FORCE AND KINEMATICS
                    if isfinite(frameAdjust) == 0
                        display([name,'   has issue with sync numbers']);
                    else
                        tstartforce = tstart + 1 - frameAdjust;  %%%ADJUST FOR TAKING FRAME DIFFERENCES BY ADDING 1 TO THE FRAME INDEX
                        tendforce = tend + 1 - frameAdjust;
                        
                        fnamex = [directoryforcelocations,name,'_Params_Matrix_ForceX.mat'];
                        fnamey = [directoryforcelocations,name,'_Params_Matrix_ForceY.mat'];
                        if exist(fnamex,'file') == 0 || exist(fnamey,'file') == 0
                            display(['no forces for ',list(ii).name]);
                        else
                            load(fnamex);
                            load(fnamey);
                            [npegs,ntimes] = size(Fx);
                            k = strfind(name,'16');
                            k = k(end);
                            date = name(k-4:k+1);
                            ppair = pegpair;
                            if any(strcmp(date,dmat))
                            else
                                %                         pegpair = sort(npegs+1-pegpair);
                                pegpair = 8 - pegpair;
                            end
                            timind = tstartforce:tendforce;
                            %%%%SNAKE PEG CROSSING ONLY ACCURATE WITHIN A COUPLE
                            %%%%FRAMES, SO IT IS POSSIBLE FOR tstartforce TO BE <1
                            timind(timind<1) = [];
                            timind(timind>ntimes) = [];
                            %%%% ON 01/23/17 I DOUBLED CHECKED WITH VIDEOS. BOTH FX
                            %%%% AND FY NEED TO BE FLIPPED TO GO WITH CONVENTION.
                            %%%% SNAKE PUSHING FORWARD IS +Z. SNAKE PUSHING RIGHT
                            %%%% IS +X
                            
                            Fz1 = -Fy(pegpair(1),timind);
                            Fx1 = -Fx(pegpair(1),timind);
                            Fz2 = -Fy(pegpair(2),timind);
                            Fx2 = -Fx(pegpair(2),timind);
                            
%                             subplot(2,1,2);
%                             plot(Fx','k');hold on;plot(Fy','Color',[0.5 0.5 0.5]);
                            Fx = [Fx1;Fx2];
                            Fz = [Fz1;Fz2];
                            
%                             plot(Fx','.-','Color',[0.8 0.8 0],'LineWidth',2);hold on;plot(Fz','.-','Color',[0 0.8 0.8],'LineWidth',2);
%                             line([1,length(Fx)],[thresh thresh],'LineStyle',':','Color','k');
%                             line([1,length(Fx)],-[thresh thresh],'LineStyle',':','Color','k');hold off;
%                             drawnow;
%                             str = input('pegtouched?','s');
%                             if strcmp(str,'y')
%                                 PegTouch(count) = 1;
%                             else
%                                 PegTouch(count) = 0;
%                             end
                            %                             line([timind(1),timind(end)],[thresh, thresh])
                            %                             line([timind(1),timind(end)],-[thresh, thresh]);
                            %                             drawnow;hold off;pause(1);
                            
                            Fmag = sqrt(Fx.^2 + Fz.^2);
                            touchtime = find(Fmag(1,:)>thresh | Fmag(2,:)>thresh);
                            
                            
%                             if exist('clickringing','var') == 1&& clickringing == 1
%                                 hold off;
%                                 subplot(2,1,1);
%                                 plot(splineX,splineY,'Color',[0.7 0.7 0.7]);hold on;
%                                 viscircles([pegx,pegy],repmat(5.5,length(pegx),1),'EdgeColor','k');
%                                 viscircles([pegx(ppair),pegy(ppair)],repmat(5.5,2,1),'EdgeColor',[0 0.7 0.7]);
%                                 plot(pegx(pegpair),pegy(pegpair),'+','LineWidth',2,'Color',[0.6 0 0.8]);
%                                 plot(splineX(:,tstart),splineY(:,tstart),'Color',[1 0 0.5],'LineWidth',2);
%                                 plot(splineX(:,tend),splineY(:,tend),'Color',[0 1 0.5],'LineWidth',2);
%                                 hold off;
%                                 subplot(2,1,2);
%                                 plot(-Fx','k');hold on;plot(-Fy','k');
%                                 plot(timind,Fx1,'LineWidth',3)
%                                 plot(timind,Fz1,'LineWidth',3);hold on;
%                                 plot(timind,Fx2,'LineWidth',3)
%                                 plot(timind,Fz2,'LineWidth',3);drawnow;
%                                 str = input('Ringing?','s');
%                                 switch str
%                                     case 'y'
%                                         hold off
%                                         plot(1:length(timind),Fx1,1:length(timind),Fz1);hold on;
%                                         plot(1:length(timind),Fx2,1:length(timind),Fz2);drawnow;
%                                         [xind,~] = ginput(2);
%                                         timind(xind(1):xind(2)) = [];
%                                         Fz1 = -Fy(pegpair(1),timind);
%                                         Fx1 = -Fx(pegpair(1),timind);
%                                         Fz2 = -Fy(pegpair(2),timind);
%                                         Fx2 = -Fx(pegpair(2),timind);
%                                         hold off;plot(timind,Fx1,'.',timind,Fz1,'.');
%                                         plot(timind,Fx2,'.',timind,Fz2,'.');drawnow;
%                                         str = input('better?','s');
%                                         RingInd(count,:) = xind;
%                                         if strcmp(str,'n') == 1
%                                             keyboard
%                                         end
%                                 end
%                             end
                            SnakeForces(count).Peg1 = [Fx1;Fz1];
                            SnakeForces(count).Peg2 = [Fx2;Fz2];
                            totaltouch(count) = length(timind);
                            forcetouch(count) = length(touchtime);
                            count = count+1;
                        end
                        clear x y splineX splineY pegXY angles kappas timind
                    end
                end
            end
        end
    end
end