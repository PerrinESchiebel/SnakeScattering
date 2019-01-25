function [x,y,endpoints] = TrackForwardandBackward(snakes)
square = [270 800 280 475];
[~,~,numFrames] = size(snakes);
endpoints = cell(1,numFrames);
middleframe = round(numFrames/2);
test = snakes(:,:,middleframe);
test(1:square(3),:) = 0;
test(square(4):end,:) = 0;
test(:,1:square(1)) = 0;
test(:,square(2):end) = 0;
CC = bwconncomp(test,4);
if CC.NumObjects > 1
    numPixels = cellfun(@numel,CC.PixelIdxList);
    idx = find(any([numPixels > 6000;numPixels < 1000]));
    for jj = 1:length(idx)
        test(CC.PixelIdxList{idx(jj)}) = 0;
    end
end
CC = bwconncomp(test,4);
XC = regionprops(CC,'centroid');
XC = cat(1,XC.Centroid);
 


% try again
test = snakes(:,:,middleframe);
    CC = bwconncomp(test,4);
XC = regionprops(CC,'centroid');
XC = cat(1,XC.Centroid);

%%%%%TRack backwards and forwards
XCmid = XC;
clear x y
x = cell(1,1);y = cell(1,1);
for jj = fliplr(1:middleframe)
    test = snakes(:,:,jj);
    CC = bwconncomp(test,4);
    nobj = CC.NumObjects;
    if nobj == 0
    elseif nobj > 1
        XC2 = regionprops(CC,'centroid');
        XC2 = cat(1,XC2.Centroid);
        dists = sqrt((XC2(:,1)-repmat(XC(1),nobj,1)).^2 + (XC2(:,2)-repmat(XC(2),nobj,1)).^2);
        ind = find(dists==min(dists));
        list = 1:nobj;
        list(ind) = [];
        for ll = list
            test(CC.PixelIdxList{ll}) = 0;
        end
        XC = XC2(ind,:);
        clear XC2 CC L idx
    elseif nobj == 1
        XC = regionprops(CC,'centroid');
        XC = cat(1,XC.Centroid);
        
    end
    skels = bwmorph(test,'thin',Inf);
    [rr,cc] = find(bwmorph(skels == 1,'endpoints'));
    endpoints{jj} = [cc,rr];
    [y{jj},x{jj}] = find(skels==1);
end
XC = XCmid;
for jj = middleframe+1:numFrames
    test = snakes(:,:,jj);
    CC = bwconncomp(test,4);
    nobj = CC.NumObjects;
    if nobj == 0
    elseif CC.NumObjects > 1
        XC2 = regionprops(CC,'centroid');
        XC2 = cat(1,XC2.Centroid);
        dists = sqrt((XC2(:,1)-repmat(XC(1),nobj,1)).^2 + (XC2(:,2)-repmat(XC(2),nobj,1)).^2);
        ind = find(dists==min(dists));
        list = 1:nobj;
        list(ind) = [];
        for ll = list
            test(CC.PixelIdxList{ll}) = 0;
        end
        XC = XC2(ind,:);
        clear XC2 CC L idx
    elseif nobj == 1
        XC = regionprops(CC,'centroid');
        XC = cat(1,XC.Centroid);
        
    end
    skels = bwmorph(test,'thin',Inf);
    [cc,rr] = find(bwmorph(skels == 1,'endpoints'));
    endpoints{jj} = [rr,cc];
    [y{jj},x{jj}] = find(skels==1);
end
end