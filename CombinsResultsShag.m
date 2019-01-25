%%last update 062617
%%%This is nearly identical to CombinsResultsShag_v2. Only minor
%%%improvements for neatness
%%%091118 added a snake number variable
function [xtail,ztail,xb,zb,xhead,zhead,xcm,zcm,xcm0,zcm0,xMatsstf,zMatsstf,pegxstf,pegystf,theta,reverses,SnakeNumbers] = CombinsResultsShag(anglemax,kappamin,snakenumber)
curdir = pwd;
% curdir = curdir(1:end-16);
curdir = curdir(1:40);
directoryresults = [curdir,'Data\snake_mats_for_combinsResults\'];
listresults = dir([directoryresults,'*.mat']);
load([directoryresults,listresults(1).name]);
if exist('ww','var')
    error('You may have an issue with the indexing variable already existing in the loaded files');
end
count = 1;
% centertocenter = 110;  %%px
% centertocenter = 8;
centertocenter = 19;   %%px for shag
initialbox = [-centertocenter/2 centertocenter/2];
xMatsstf = cell(1,1);
zMatsstf = cell(1,1);
pegxstf = cell(1,1);
pegystf = cell(1,1);
theta = cell(1,1);
rho = cell(1,1);
reverses = nan(100,1);
SnakeNumbers = nan(100,1);
yesno = 1;
reverse = 'false';
% pointsinavg = 1:50; %%I tried referencing this to find meany e.g.
% splineX(pointsinavg,1:pegxing) but it is much slower than just typing
% 1:50 or : or whatever directly in
for ww = 1:length(listresults)
    fprintf(1,'Trial number %i \n', ww);
    load([directoryresults,listresults(ww).name]);
    if exist('snakenumber','var')
        if strcmp(name(1),'1')==0
            snnum = str2double(name(1:2));
            snnum = 100+snnum;
        else
            snnum = str2double(name(1:3));
        end
        yesno = any(snnum == snakenumber);
    end
    if yesno == 1
        %         display('yes')
        if exist('angles','var') && exist('kappas','var')
            %             figure(9);plot(ww,abs(rad2deg(angles)),'o');hold on;drawnow
            %             figure(10);plot(ww,kappas,'o');hold on;drawnow
            if abs(atand(angles)) < anglemax && kappas > kappamin
                    SnakeNumbers(ww) = str2double(listresults(ww).name(1:3));
                %                 display('yes')
%                 pegR = 4;
%                 diameterInches = 0.25;
%                 rs = mean(pegR(pegR>0));
%                 conversion = diameterInches/2*2.54/mean(rs);
                pegXY = sortrows(pegXY,2);
%                 splineX = splineX - mean((pegXY(:,1)));
%                 splineY = splineY - mean((pegXY(:,2)));
%                 pegx = pegXY(:,1) - mean((pegXY(:,1)));
%                 pegy = pegXY(:,2) - mean((pegXY(:,2)));
                splineX = splineX - ((pegXY(4,1)));
                splineY = splineY - ((pegXY(4,2)));
                pegx = pegXY(:,1) - ((pegXY(4,1)));
                pegy = pegXY(:,2) - ((pegXY(4,2)));
                pegxing = find(splineX(1,:) < 0,1,'last');
%                 t1 = max([1,pegxing-100]);
                meany = mean(mean(splineY(1:50,1:pegxing)));
                count2 = 1;
                while (meany > initialbox(1) && meany < initialbox(2)) == 0
                    if count2 > 15
                        fprintf(1,'not in box');
                        break
                    else
                        if meany < initialbox(1)
                            splineY = splineY+centertocenter;
                        else
                            splineY = splineY-centertocenter;
                        end
                        %                         meany = mean(mean(splineY(1:50,1:pegxing)));
                        meany = mean(mean(splineY(1:50,1:pegxing)));
                        count2 = count2+1;
                    end
                    
                end
                %                 display(count2)
                if count2 < 16
                    %                     display(count)
                    xMatsstf{count} = splineX;
                    zMatsstf{count} = splineY;
                    pegxstf{count} = pegx;
                    pegystf{count} = pegy;
                    [theta{count},rho{count}] = cart2pol(splineX,splineY);
                    if ischar(reverse)
                        if strcmp(reverse,'false')
                            reverse = 0;
                        else
                            reverse = 1;
                        end
                    end
                    reverses(count) = reverse;
                    count = count + 1;
                end
                
                
            end
            clear x y splineX splineY pegXY angles kappas
        end
    end
end

xMatsstf(cellfun('isempty', xMatsstf)) = [];
zMatsstf(cellfun('isempty', zMatsstf)) = [];
xcm = cell(size(xMatsstf));
zcm = cell(size(xMatsstf));
xb = cell(size(xMatsstf));
zb = cell(size(xMatsstf));
xcm0 = zeros(size(xMatsstf));
zcm0 = zeros(size(xMatsstf));
xhead = cell(size(xMatsstf));
zhead = cell(size(xMatsstf));
xtail = cell(size(xMatsstf));
ztail = cell(size(xMatsstf));

for ww = 1:length(xMatsstf)
    xcm{ww} = mean(xMatsstf{ww},1);
    zcm{ww} = mean(zMatsstf{ww},1);
    xhead{ww} = xMatsstf{ww}(1:10,:);
    zhead{ww} = zMatsstf{ww}(1:10,:);
    xb{ww} = xMatsstf{ww}(:);
    zb{ww} = zMatsstf{ww}(:);
    xcm0(ww) = xcm{ww}(1);
    zcm0(ww) = zcm{ww}(1);
    xtail{ww} = xMatsstf{ww}(1:100,:)';
    ztail{ww} = zMatsstf{ww}(1:100,:)';
end
% diameterInches = 2;
% rs = mean(pegRadii(pegRadii>0));
% conversion = diameterInches/2*2.54/mean(rs);