%%%%%%%Using x,y lab frame snake kinematics, calculate the curvature along
%%%%%%%the body at each frame IN PIXELS!!!!!!!!!!
%%%anginc is how many points to skip forward and back to calculate
%%%for most videos anginc = 1 will do

%%%%%12-11-15 I checked what happens if you do 100 vs 500 points with and
%%%%%without the smoothingspline procedure. If you do 500 points you can do
%%%%%with or without. 100 points seems to capture curvature find. However,
%%%%%the smoothing spline doesn't do a great job faithfully following 100
%%%%%points.

%%%%%02-12-18 added a line to cut the nans off the front before returning
%%%dim(x) -> numpts x numframes

function Curvature = spatialCurvature(x,y,anginc,smoothparam,plotyesno)
% tic
if nargin < 5
    plotyesno = 0;
end
    [numpts,numframes] = size(x); 
    if numpts ~= 100 && numpts ~= 500
        x = x';
        y = y';
        [numpts,numframes] = size(x);
    end
    if numpts ~= 100 && numpts ~= 500
        fprintf(1,'error: number of tracked points');
        return;
    end


    Curvature = nan(numpts-anginc,numframes);
    Angles = nan(numpts-anginc,numframes);
    for mm=1:numframes
        clear pasttopresentvector presenttonextvector dotprod crossprod;
        tempx = x(:,mm);
        tempy = y(:,mm);


        ttrk=1:numpts;
        if exist('smoothparam','var')
            ft = fittype( 'smoothingspline');
            opts = fitoptions( ft );
            opts.SmoothingParam = smoothparam;  %%%%original value = 0.001
            [xData, yData] = prepareCurveData( [], tempx );
            Ftx = fit( xData, yData, ft, opts );
            [xData, yData] = prepareCurveData( [], tempy );
            Fty = fit( xData, yData, ft, opts );
            tempx=Ftx(ttrk);
            tempy=Fty(ttrk);
            if plotyesno == 1
                figure(3);plot(x(:,mm),y(:,mm));hold on;plot(tempx,tempy,'--k');axis tight;drawnow;hold off;
            end
        end

        %angle, radius and curvature all along the track
        %Calculate using formula for a fit circle
%         spcurve = zeros(1,numpts-anginc);
        SplineAngles = zeros(1,numpts-anginc)';
        SplineRadius = zeros(1,numpts-anginc)';
        SplineCurvature = zeros(1,numpts-anginc)';
        for ii=1+anginc:numpts-anginc
            pasttopresentvector= [tempx(ii)-tempx(ii-anginc),tempy(ii)-tempy(ii-anginc)];
            presenttonextvector= [tempx(ii+anginc)-tempx(ii),tempy(ii+anginc)-tempy(ii)];
%             dtan = presenttonextvector-pasttopresentvector;
%             spcurve(ii) = sqrt(dtan(1)^2+dtan(2)^2);
%             dotprod=pasttopresentvector(1)*presenttonextvector(1)+pasttopresentvector(2)*presenttonextvector(2);
            dotprod=dot(pasttopresentvector, presenttonextvector);
            dotprod=dotprod/(norm(pasttopresentvector)*norm(presenttonextvector));
            SplineAngles(ii) = acos(dotprod); 
            crossprod=pasttopresentvector(1)*presenttonextvector(2)-pasttopresentvector(2)*presenttonextvector(1);

            pastpresslope=(tempy(ii)-tempy(ii-anginc))/(tempx(ii)-tempx(ii-anginc));
            presfutslope=(tempy(ii+anginc)-tempy(ii))/(tempx(ii+anginc)-tempx(ii));
            circx=(pastpresslope*presfutslope*(tempy(ii-anginc)-tempy(ii+anginc))+presfutslope*(tempx(ii-anginc)+...
                tempx(ii))-pastpresslope*(tempx(ii)+tempx(ii+anginc)))/(2*(presfutslope-pastpresslope));
            circy=(-1*(circx-(tempx(ii-anginc)+tempx(ii))/2)/pastpresslope)+(tempy(ii-anginc)+tempy(ii))/2;
            SplineRadius(ii)=((tempx(ii)-circx)^2+(tempy(ii)-circy)^2)^0.5;
            SplineCurvature(ii)=1/SplineRadius(ii);
            if crossprod<0
                 SplineRadius(ii)=SplineRadius(ii)*-1;
                 SplineCurvature(ii)=SplineCurvature(ii)*-1;
                 SplineAngles(ii)=SplineAngles(ii)*-1;
            end
        end           
    Curvature(:,mm) = SplineCurvature;
    Angles(:,mm) = SplineAngles;
    end  
%     toc
Curvature(1:anginc,:) = [];
end