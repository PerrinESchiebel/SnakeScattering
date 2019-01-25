figure
nboots = 10000;
bins = -62.5:5:62.5;
bb = bins(1:end-1)+diff(bins)/2;
colors = [[1 0 0];[1 0.3 0];[1 1 0];[0 1 0];[0 0 1]];
AllTrue = nan(1,nboots);
for jj=1:nboots
% for jj=9:nboots
    [~,peaklocs] = findpeaks(snakeboots(jj,:),'MinPeakProminence',0.004); % from real dist, the smallest prominence is 0.0052
    [~,valleylocs] = findpeaks(-snakeboots(jj,:),'MinPeakProminence',0.004);
%     plot(bb,snakeboots(jj,:),'LineWidth',3,'Color','k');hold on;
    if length(peaklocs) >= 3 && length(valleylocs)>=2
%         centerpeak = peaklocs(bb(peaklocs)>=-5 & bb(peaklocs)<7);
        
        centerpeak = peaklocs(bb(peaklocs)>-10 & bb(peaklocs)<15);
        if length(centerpeak)>1
            %             centerpeak = round(mean(centerpeak));
            centerpeak = centerpeak(bb(centerpeak)==max(bb(centerpeak)));
        end
        leftpeak = peaklocs(bb(peaklocs)>=-30 & bb(peaklocs)<-10);
        if length(leftpeak)>1
            %             leftpeak = round(mean(leftpeak));
            leftpeak = leftpeak(bb(leftpeak)==max(bb(leftpeak)));
        end
%         rightpeak = peaklocs(bb(peaklocs)>15 & bb(peaklocs)<=30);
                rightpeak = peaklocs(bb(peaklocs)>15 & bb(peaklocs)<=35);
        if length(rightpeak)>1
            %             rightpeak = round(mean(rightpeak));
            rightpeak = rightpeak(bb(rightpeak)==max(bb(rightpeak)));
        end
        if isempty(centerpeak) || isempty(leftpeak) || isempty(rightpeak)
            AllTrue(jj) = 0;
%             fprintf(1,'All peaks not found run %i \n',jj);
            title('all peaks not found');
        else
            middlevalley = bb(valleylocs)>-2 & bb(valleylocs)<10;
            if sum(middlevalley) > 0
%                 fprintf(1,'Middle Valley in run %i \n',jj);
                title('middle valley found');
            else
                %             peaklocs = peaklocs(centerpeak-1:centerpeak+1);
                
                %      centervalleys = find(bb(valleylocs)>-20 & bb(valleylocs)<20);
                leftvalley = valleylocs(valleylocs>leftpeak & valleylocs<centerpeak);
%  leftvalley = valleylocs(bb(valleylocs)>-20 & bb(valleylocs)<5);        
                if length(leftvalley)>1
                    %                     leftvalley = round(mean(leftvalley));
                               
                    leftvalley = leftvalley(bb(leftvalley)==min(bb(leftvalley)));
                end
%                 rightvalley = valleylocs(bb(valleylocs)>5 & bb(valleylocs)<20);
                rightvalley = valleylocs(valleylocs>centerpeak & valleylocs<rightpeak);
                if length(rightvalley)>1
                    %                     rightvalley = round(mean(rightvalley));
                    
                    rightvalley = rightvalley(bb(rightvalley)==min(bb(rightvalley)));
                end
                if isempty(leftvalley) || isempty(rightvalley)
                    AllTrue(jj) = 0;
%                     fprintf(1,'All valleys not found run %i \n', jj);
                    title('all valleys not found');
                else
                    %         valleylocs = valleylocs(centervalleys);
                    %                 valleylocs = centervalleys;
                    p1 = snakeboots(jj,leftpeak);
                    p3 = snakeboots(jj,centerpeak);
                    p5 = snakeboots(jj,rightpeak);
                    p2 = snakeboots(jj,leftvalley);
                    p4 = snakeboots(jj,rightvalley);
                    
%                     plot(bb,snakeboots(jj,:));hold on;
%                     indices = [leftpeak, leftvalley, centerpeak, rightvalley, rightpeak];
%                     scatter(bb(indices),snakeboots(jj,indices),100,1:5,'o','filled');colormap(colors);
                    %          plot(bb(peaklocs),snakeboots(jj,peaklocs),'*');
                    %          plot(bb(valleylocs),snakeboots(jj,valleylocs),'s');
                    %          text(bb(peaklocs(1)),p1,'p1');
                    %          text(bb(peaklocs(2)),p3,'p3')
                    %          text(bb(peaklocs(3)),p5,'p5')
                    %          text(bb(valleylocs(1)),p2,'p2')
                    %          text(bb(valleylocs(2)),p4,'p4')
                    %                     drawnow;pause(0.5);hold off;
                    if p1>p2 && p2<p3 && p3>p4 && p4<p5
                        AllTrue(jj) = 1;
                    else
                        AllTrue(jj) = 0;
                    end
                end
            end
        end
    end
%     drawnow;hold off;
end
 sum(AllTrue,'omitnan')./length(AllTrue)