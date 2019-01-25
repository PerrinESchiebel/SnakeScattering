function AllCurve = CreateMatrix4PCA_SHAG_v2(DDdirectory,CurveOrTan,mperpx,plotyesno)
nsplinepts = 100;
if nargin < 3
    %     px2cm = 12;
    CurveOrTan = 'Curve';
end
anginc = 5;  %%%%%This is ideal for finding curvature 100 pts on Chionactis
seginc = 3;  %%%%%This is ideal for finding tangent angle 100 pts on Chionactis
listlist = dir([DDdirectory,'*.mat']);
numfiles = size(listlist,1);
AllCurve = cell(numfiles,1);
segin = 1;
segout = 100;
for qq = 1:numfiles
    if qq == 11 || qq == 23 || qq==34 ||qq==41 || qq==52
    else
        filename = listlist(qq).name;
        display(filename);
        %     clear px2cm;
        %     px2cm = conversions(strcmp(filename,{conversions.name})).conversionFactor;
        load(strcat(DDdirectory,'',filename));
        if qq == 11 || qq == 23 || qq==34 ||qq==41 || qq==52
        else
            xs = nan(nsplinepts,length(snakeX));
            ys = nan(nsplinepts,length(snakeY));
            % xs = nan(numsplinepts,length(snakeX));
            % ys = nan(numsplinepts,length(snakeX));
            ntimes = length(snakeX);
            % ntimes = length(snakeX);
            for jj=1:ntimes
                x = snakeX{jj};
                y = snakeY{jj};
                %     x = snakeX{jj};y = snakeY{jj};
                if isempty(x) == 0 && length(x)>50
                    if x(1)==x(end)
                        x(end) = x(end)+1;
                    elseif y(1)==y(end)
                        y(end) = y(end)+1;
                    end
                    [psx,psy,~] =  cubesplineinterp(x,y,nsplinepts,'space',0,1);
                    xs(:,jj) = psx;
                    ys(:,jj) = psy;
                end
            end
            % xs = xs(5:end-4,30:end);
            % ys = ys(5:end-4,30:end);
            [nsplinepts,ntimes] = size(xs);
            xss = nan(nsplinepts,ntimes);
            yss = nan(nsplinepts,ntimes);
            for jj=1:nsplinepts
                xss(jj,:) = smooth(xs(jj,:),7);
                yss(jj,:) = smooth(ys(jj,:),7);
            end
            xx = flipud(xss);
            yy = flipud(yss);
            if qq == 1
                fct = 70;
            elseif qq == 2
                fct = 1;
            elseif qq == 3
                fct = 100;
            elseif qq == 6
                fct = 50;
            elseif qq == 7 || qq == 17 || qq == 12 || qq == 15 || qq ==18 || qq==26 || qq==27 || ...
                    qq==29 || qq==33 || qq==36 || qq==37 || qq==38 || qq==40 ||qq==47
                fct = 30;
            elseif qq == 20
                fct = 100;
            elseif qq ==45 || qq==46
                fct = 50;
                %     fct = 70;
            else
                %     fct = round(nframes/10);
                fct = 1;
            end
            if qq==44
%                  xout = splineX(segin:segout,fct:end-20);
%                 yout = splineY(segin:segout,fct:end-20);
%                 xout(xout<100) = NaN;
%                 yout(yout<100) = NaN;
                xout = xx(segin:segout,fct:end-20);
                yout = yy(segin:segout,fct:end-20);
                xout(xout<100) = NaN;
                yout(yout<100) = NaN;
%                     y = y(x>100);
%                     x = x(x>100);
                
            else
                xout = xx(segin:segout,fct:end);
                yout = yy(segin:segout,fct:end);
            end
            
            
            
        end
 if plotyesno == 1
            subplot(2,1,1);
                plot(cell2mat(snakeX),cell2mat(snakeY),'.k');hold on;plot(xout,yout,'.','Color',[0.3 0.3 0.8]);axis tight;hold off;
 end

        x = xout;
        y = yout;
        
        %%
        [ns,~] = size(x);
        if ns == 500
            x = x(1:5:500,:);
            y = y(1:5:500,:);
        end
        switch CurveOrTan
            case 'Curve'
                clear Curvature;
                Curvature = spatialCurvature(x,y,anginc); %%%%NOTE THIS IS RETURNED IN PIXELS
                Curvature(Curvature(:,1)==0,:)=[];
                if length(mperpx) == 1
                    Curvature = Curvature./mperpx;
                else
                    Curvature = Curvature./mperpx(qq);
                end
                Curvature(abs(Curvature)>200) = NaN;
                Curvature(isnan(Curvature)==1) = max(max(abs(Curvature)));
                AllCurve{qq} = Curvature;
            case 'Tan'
                clear Tangent;
                theta = tangentangle3(x,y,seginc);
                theta(theta(:,1)==0,:)=[];
                theta(isnan(theta)==1) = max(max(abs(theta)));
                AllCurve{qq} = theta;
        end
        if plotyesno == 1
                subplot(2,1,2);pcolor(Curvature);shading flat;colormap(redblue);colorbar;axis square;drawnow;
        end

    end
    %     pcolor(rad2deg(AllCurve{q}));
    %     shading flat
    %     drawnow
end