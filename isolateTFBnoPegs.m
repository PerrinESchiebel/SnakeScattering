function [test,blah,snake] = isolateTFBnoPegs(frame1,frame2,pegImage,backgroundThreshold,seSize,radiusRange,minThresh,maxThresh)
% % % % % % % % % % % % % % % isolateTFB(frame1,frame2,pegImage,backgroundThreshold,seSize,radiusRange,minThresh,maxThresh,iter)

imageDiff = frame2-frame1;

% find peg outlines if they are not already known
if nargin < 6 || isempty(radiusRange)
    %radiusRange = [30 45]; % 2 in pegs
    radiusRange = [11 14]; % 1 in pegs
end

if nargin < 7 || isempty(minThresh)
    minThresh = 100;
end
if nargin < 8 || isempty(maxThresh)
    maxThresh = 600;
end

if nargin < 3 || isempty(pegImage)
    
       
%     [xx,yy] = meshgrid(1:size(imageDiff,2),1:size(imageDiff,1));
%     
%     [pegCenters,pegRadii] = imfindcircles(frame1,radiusRange,'ObjectPolarity','bright','method','twostage');
%     pegImage = zeros(size(imageDiff));
%     
%     thetas = linspace(0,2*pi,1000);
%     for i = 1:length(pegRadii)
%         pegImage = pegImage + inpolygon(xx,yy,pegCenters(i,1)+pegRadii(i)*cos(thetas),pegCenters(i,2)+pegRadii(i)*sin(thetas));
%     end
%     
% else [pegCenters,pegRadii] = imfindcircles(pegImage,radiusRange,'ObjectPolarity','bright');
    pegImage = 0;
end

if nargin < 4 || isempty(backgroundThreshold)
    backgroundThreshold = 15;
end
if nargin < 5 || isempty(seSize)
    seSize = 25;
end



% test = -(imageDiff - double(pegImage).*imageDiff);
test = -(imageDiff - 0);
test(test <= backgroundThreshold) = 0;

% starting oct 3, poorly behaved pixel at (x,y) = (220,576)?  median filter
% seems to fix it
test = medfilt2(test);
blah = imdilate(test,strel('disk',seSize));
blah(blah > 0) = 1;
% bw = imclose(blah,strel('disk',15));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
snake = bwmorph(blah,'thin',inf);
% % % 
% % % snakes = test;
% % % [x,y,endpoints] = TrackForwardandBackward(snakes);
% % % 
% % % 
% % % % throw away regions that are too small
% % % % rp = regionprops(snake,'PixelIdxList');
% % % % 
% % % % for i = 1:length(rp)
% % % %     if length(rp(i).PixelIdxList) < minThresh
% % % %         snake(rp(i).PixelIdxList) = 0;
% % % %     end
% % % %     
% % % %     if length(rp(i).PixelIdxList) > maxThresh
% % % %         snake(rp(i).PixelIdxList) = 0;
% % % %     end
% % % % end
% % % % 
% % % % 
% % % % [snakeY,snakeX] = find(snake == 1);
% % % snakeX = x; snakeY = y;
% % % 
% % % snakeXY = [snakeX snakeY];
% % % 
% % % [rr,cc] = find(bwmorph(snake == 1,'endpoints'));
% % % endPoints = [cc,rr];
% % % 
% % % % if isempty(snakeX) && iter > 100
% % % %     imagesc(imageDiff)
% % % %     keyboard
% % % % end
% % %  