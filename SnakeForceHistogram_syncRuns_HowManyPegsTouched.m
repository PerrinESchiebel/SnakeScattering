%%%Created from SnakeForceHistogram_syncRuns on 
%%%092617
%%%added option to not judge accuracty on 022718
function [SnakeForces,RingInd] = SnakeForceHistogram_syncRuns_HowManyPegsTouched(anglemax,kappamin,clickringing,judgeruns)
%%
%%%%%%%%%%%MAKE A STRUCT WITH ALL RUN NAMES AND MEASURED THETA AND KAPPA%%%%%%%%%%%%%
%%%%%USED anglemax = 15;kappamin = 30;
directoryforcehist = 'F:\Dropbox\SnakeScattering\Figures\snake_codes\snake_mats_for_combinsResults\';  %%%% TRACKED SNAKES AFTER SNAKE RUN SORTING
listforcehist = dir([directoryforcehist,'*.mat']);
directoryforcelocations = 'F:\Dropbox\SnakeScattering\Figures\snake_codes\allForces\';   %%%%FORCES FROM TRACKED PEGS
[num,txt,~] = xlsread('F:\Dropbox\Spreadsheets\vidSync.xlsx');      %%%%BY HAND, CHOSE MATCHING FRAMES IN THE TWO VIDEOS
syncFileNames = txt(2:end,2);
syncFrameNum = num(:,[3,5]);  %%%%First number is for force data, second is for tracked data
clear num txt
nforceruns = length(listforcehist);
SnakeForces = struct('Peg1',[],'Peg2',[],'pegstouched',[],'x',[],'y',[]);
virtualR = 6.5 + 10;   %%%%%From ImageJ measurements: peg diameter ~ 5.5 px. snake diameter ~ 7.5px, so (pegD + snakeD)/2 = 6.5
%%%%% 012317, I noticed sometimes
%%%%% the time is a little late
%%%%% because of the tracking not
%%%%% getting all the way to the
%%%%% head so I added 3 px
threshold = 0.005;
centertocenter = 19;   %%px for shag
initialbox = [-centertocenter/2 centertocenter/2];
count = 1;
dmat = [{'081616'},{'082216'},{'090116'},{'082316'}];
RingInd = nan(100,2);
anglemax = 15;kappamin = 30;
%%%%  THE FIRST MATRIX IS BROKEN, START AT 2
for mq=1:nforceruns
%     display(mq)
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
            if truefalse ==0
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
                %%%%THE PEG DATA FROM THE FIRST DAY TRACKED THE PEGS IN REVERSE
                %%%%ORDER, SO CHECK FOR THIS AND CORRECT
                
                %%%%FIND WHEN HEAD FIRST CROSSES LINE OF PEGS AND TAIL LAST
                %%%%CROSSES
                tstart = find(splineX(1,:) > 0 - virtualR,1,'first');
                tend = find(splineX(end,:) > 0 + virtualR,1,'first');
                if isempty(tstart) || isempty(tend)
                    display([name,'   not enough snake']);
                else
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
                            if exist('clickringing','var') == 1&& clickringing == 1
                                hold off;
                                subplot(2,1,1);
                                plot(splineX,splineY,'Color',[0.7 0.7 0.7]);hold on;
                                viscircles([pegx,pegy],repmat(5.5,length(pegx),1),'EdgeColor','k');
                                viscircles([pegx(ppair),pegy(ppair)],repmat(5.5,2,1),'EdgeColor',[0 0.7 0.7]);
                                plot(pegx(pegpair),pegy(pegpair),'+','LineWidth',2,'Color',[0.6 0 0.8]);
                                plot(splineX(:,tstart),splineY(:,tstart),'Color',[1 0 0.5],'LineWidth',2);
                                plot(splineX(:,tend),splineY(:,tend),'Color',[0 1 0.5],'LineWidth',2);
                                hold off;
                                subplot(2,1,2);
                                plot(-Fx','k');hold on;plot(-Fy','k');
                                plot(timind,Fx1,'LineWidth',3)
                                plot(timind,Fz1,'LineWidth',3);hold on;
                                plot(timind,Fx2,'LineWidth',3)
                                plot(timind,Fz2,'LineWidth',3);drawnow;
                                str = input('Ringing?','s');
                                switch str
                                    case 'y'
                                        hold off
                                        plot(1:length(timind),Fx1,1:length(timind),Fz1);hold on;
                                        plot(1:length(timind),Fx2,1:length(timind),Fz2);drawnow;
                                        [xind,~] = ginput(2);
                                        timind(xind(1):xind(2)) = [];
                                        Fz1 = -Fy(pegpair(1),timind);
                                        Fx1 = -Fx(pegpair(1),timind);
                                        Fz2 = -Fy(pegpair(2),timind);
                                        Fx2 = -Fx(pegpair(2),timind);
                                        hold off;plot(timind,Fx1,'.',timind,Fz1,'.');
                                        plot(timind,Fx2,'.',timind,Fz2,'.');drawnow;
                                        str = input('better?','s');
                                        RingInd(count,:) = xind;
                                        if strcmp(str,'n') == 1
                                            keyboard
                                        end
                                end
                            end
                            Fmag1 = sqrt(Fx1.^2+Fz1.^2);
                            Fmag2 = sqrt(Fx2.^2+Fz2.^2);
                            ind1 = find(Fmag1>threshold);
                            ind2 = find(Fmag2>threshold);
                            num1 = sum(Fmag1>threshold);
                            num2 = sum(Fmag2>threshold);
                            if num1 > 5 && num2 < 5 || num1 < 5 && num2 >5
                                numpegstouched = 1;
                            elseif num1>10 && num2>10
                                numpegstouched = 2;
                            else
                                numpegstouched = 0;
                            end
                            
                            if exist('judgeruns','var') == 1&& judgeruns == 1
                                colors  = colormap(jet(size(Fx,1)));
                                subplot(3,1,3);
                                plot(Fmag1,'LineWidth',2,'Color',colors(pegpair(1),:));hold on;plot(Fmag2,'LineWidth',2,'Color',colors(pegpair(2),:))
                                plot(ind1,Fmag1(ind1),'o');plot(ind2,Fmag2(ind2),'o');
                                hold off
                                subplot(3,1,1);
                                plot(splineX,splineY,'Color',[0.7 0.7 0.7]);hold on;
                                viscircles([pegx,pegy],repmat(5.5,length(pegx),1),'EdgeColor','k');
                                viscircles([pegx(ppair),pegy(ppair)],repmat(5.5,2,1),'EdgeColor',[0 0.7 0.7]);axis equal tight;hold off;
                                title(numpegstouched,'FontSize',32);
                                subplot(3,1,2);
                                for stx = 1:length(colors)
                                    plot(Fx(stx,:),'Color',colors(stx,:),'LineWidth',2);hold on;
                                    plot(Fy(stx,:),'--','Color',colors(stx,:),'LineWidth',2);
                                end
                                hold off;
                                drawnow;
                                str = input('OK? type n if no','s');
                                switch str
                                    case 'n'
                                        numpegstouched = input('How many pegs?');
                                end
                            end
                            
                            SnakeForces(count).Peg1 = [Fx1;Fz1];
                            SnakeForces(count).Peg2 = [Fx2;Fz2];
                            SnakeForces(count).pegstouched = numpegstouched;
                            %%
                            %%
                            %                 display('yes')
                            %                 pegR = 4;
                            %                 diameterInches = 0.25;
                            %                 rs = mean(pegR(pegR>0));
                            %                 conversion = diameterInches/2*2.54/mean(rs);
                            %                 pegXY = sortrows(pegXY,2);
                            %                 splineX = splineX - mean((pegXY(:,1)));
                            %                 splineY = splineY - mean((pegXY(:,2)));
                            %                 pegx = pegXY(:,1) - mean((pegXY(:,1)));
                            %                 pegy = pegXY(:,2) - mean((pegXY(:,2)));
                            pegxing = find(splineX(1,:) < 0,1,'last');
                            meany = mean(mean(splineY(1:50,1:pegxing)));
                            count2 = 1;
                            while (meany > initialbox(1) && meany < initialbox(2)) == 0
                                if count2 > 15
                                    display('not in box');
                                    break
                                else
                                    if meany < initialbox(1)
                                        splineY = splineY+centertocenter;
                                    else
                                        splineY = splineY-centertocenter;
                                    end
                                    meany = mean(mean(splineY(:,1:pegxing)));
                                    count2 = count2+1;
                                end
                                
                            end
                            %                 display(count2)
                            if count2 < 16
                                SnakeForces(count).x = splineX;
                                SnakeForces(count).y = splineY;
                                display(mq);
                                %                     display(count)
                                %                     xMatsstf{count} = splineX;
                                %                     zMatsstf{count} = splineY;
                                %                     pegxstf{count} = pegx;
                                %                     pegystf{count} = pegy;
                                %                     [theta{count},rho{count}] = cart2pol(splineX,splineY);
                                %                     if ischar(reverse)
                                %                         if strcmp(reverse,'false')
                                %                             reverse = 0;
                                %                         else
                                %                             reverse = 1;
                                %                         end
                                %                     end
                                %                     reverses(count) = reverse;
                                %                     count = count + 1;
                            end
                            
                        end
                        count = count+1;
                    end
                    clear x y splineX splineY pegXY angles kappas timind
                end
            end
        end
    end
end