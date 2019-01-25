function [snakeX,snakeY,splineX,splineY,pegCenters,pegRadii,keep] = trackSnakeInVideoOnShagNewTFB(videoFile,backgroundThreshold,seSize,minThresh,maxThresh,radiusRange,plotOn)

cropLeft = true;
if ~cropLeft
    cropLeftPix = 1;
else
    cropLeftPix = 200;
end

if nargin < 2 || isempty(backgroundThreshold)
    backgroundThreshold = 8;
end
if nargin < 3 || isempty(seSize)
    seSize = 10;
end
if nargin < 4 || isempty(minThresh)
    minThresh = 100;
end
if nargin < 5 || isempty(maxThresh)
    maxThresh = 500;
end
if nargin < 6 || isempty(radiusRange)
    radiusRange = [11 14];
end
if nargin < 7 || isempty(plotOn)
    plotOn = true;
    h = figure(10);clf;
    saveVideo = false;
end

dThresh = 10;
sortDir = '+';
savePath = 'D:\Jennifer\shag\';
if ~exist(savePath,'dir')
mkdir(savePath);
end

vid = VideoReader(videoFile);
frame1 = double(rgb2gray(read(vid,1)));

frame1 = frame1(:,cropLeftPix:end);

[xx,yy] = meshgrid(1:size(frame1,2),1:size(frame1,1));

[pegCenters,pegRadii] = imfindcircles(frame1,radiusRange,'ObjectPolarity','bright','method','twostage');
if ~isempty(pegCenters)
    pegImage = zeros(size(frame1));
    
    thetas = linspace(0,2*pi,1000);
    for i = 1:length(pegRadii)
        pegImage = pegImage + inpolygon(xx,yy,pegCenters(i,1)+1.25*pegRadii(i)*cos(thetas),pegCenters(i,2)+1.25*pegRadii(i)*sin(thetas));
    end
else
    disp('No pegs identified, try increasing radius range.')
    keyboard
    %     test = imfill(imgradient(frame1));
    %     test(test>75) = 75;
    %     test(test<75) = 0;
    %     test(test~=0) = 1;
    %     pegImage = imdilate(imopen(imfill(test),strel('disk',10)),strel('disk',4));
    %     rp = regionprops(logical(pegImage),'Centroid','MajorAxisLength','MinorAxisLength');
    %     for i = 1:length(rp)
    %         pegCenters(i,:) = rp(i).Centroid;
    %         pegRadii(i) = mean([rp(i).MajorAxisLength rp(i).MinorAxisLength])/2;
    %     end
end






numFrames = floor(vid.Duration*vid.FrameRate);
snakeXY = cell(numFrames-1,1);
endpoints = cell(numFrames-1,1);
snakeX = cell(numFrames-1,1);
snakeY = cell(numFrames-1,1);
keep = zeros(numFrames-1,1);

[ll ww]  = size(frame1);
snakes = zeros(ll,ww,numFrames-1);

for j = 1:(numFrames-1)
    
    frame2 = double(rgb2gray(read(vid,j+1)));
    
    frame2 = frame2(:,cropLeftPix:end);
    
    %make image of pegs to subtract from images
    pegImage = zeros(size(frame1));
    
    thetas = linspace(0,2*pi,1000);
    for i = 1:length(pegRadii)
        pegImage = pegImage + inpolygon(xx,yy,pegCenters(i,1)+1.5*pegRadii(i)*cos(thetas),pegCenters(i,2)+1.5*pegRadii(i)*sin(thetas));
    end
    
    if mod(j,ceil(numFrames/4)) == 0  %mod(j,100)==0
        disp(['Frame ' num2str(j) ' / ' num2str(numFrames-1)]);
%         break
    end
    
%     [snakeXY{j},~,~,endpoints{j}] = isolateSnakeAndPegsOnShagTFB(frame1,frame2,pegImage,backgroundThreshold,seSize,radiusRange,minThresh,maxThresh);
        [test,blah,snake] = isolateTFB(frame1,frame2,pegImage,backgroundThreshold,seSize,radiusRange,minThresh,maxThresh);

        snakes(:,:,j) = blah;
        
         frame1 = frame2;
        
end

[x,y,endpoints] = TrackForwardandBackward(snakes);
for j = 1:(numFrames-1)
    try
   snakeXY{j}(:,1) = x{j};
   snakeXY{j}(:,2) = y{j};
    catch me
    end
end

for j = 1:(numFrames-1)

    if ~isempty(snakeXY{j})  
        try 
            testttt = (endpoints{j}(1,:));
%     if ~isempty(snakeXY{j})
            
        keep(j) = 1;
        [xSort,ySort] = sortSnakeDataPoints_PES(snakeXY{j}(:,1),snakeXY{j}(:,2),dThresh,endpoints{j}(1,:));
%                 [xSort,ySort] = sortSnakeDataPoints_PES(snakeXY{j}(:,1),snakeXY{j}(:,2),dThresh,sortDir,endpoints{j}(1,:));
        catch me2
            disp([me2.message  ' j = ' (num2str(j))]);
        end
%         [xSort,ySort] = sortSnakeDataPoints_PES(x,y,dThresh,sortDir,endpoints{j}(1,:));
        
        snakeX{j} = xSort;
        snakeY{j} = ySort;
        
        
        
        if plotOn %&& mod(j,100) == 0
            imagesc(frame1),title(num2str(j)),colormap(gray),hold on
            plot(snakeXY{j}(:,1),snakeXY{j}(:,2),'o','color','b','markersize',2,'markerfacecolor','b')
            %plot(splineX(:,j),splineY(:,j),'color','b','linewidth',2)
            viscircles(pegCenters,pegRadii,'edgecolor','r','linewidth',2,'DrawBackgroundCircle',false);
            axis equal tight;
            set(gca,'linewidth',2,'fontsize',20,'fontname','helvetica','Position',[0.1300 0.1100 0.7750 0.8150])
            if saveVideo
                set(h,'PaperPositionMode','auto')
                print([savePath num2str(j,'%03i')],'-dpng','-r0')
            end
            pause(0.02);
        end
        hold off
    else  snakeX{j} = [];
        snakeY{j} = [];
        
    end
    
    
   
end



snakeX = snakeX(keep==1);
snakeY = snakeY(keep==1);
[splineX,splineY] = MMS_bspline_skel(snakeX,snakeY,200);


