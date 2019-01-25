function [xSort,ySort] = sortSnakeDataPoints_PES(x,y,dThresh,startPoint)


if nargin < 3 || isempty(dThresh)
    
    dThresh = 10;
end




% sorting pixels based on closest distance

if size(x,1) < size(x,2)
    x = x';
    y = y';
end


% find distances b etween all particle pairs
xx = repmat(x,1,length(x));
yy = repmat(y,1,length(y));

 D = sqrt((xx-xx').^2+(yy-yy').^2);
% A = sparse(D < 2);
% Z = graphallshortestpaths(A);
% [ii,jj] = find(Z == max(Z(~isinf(Z))));
% 


% if qqq == 98
%     keyboard
%     end

% set diagonal terms to inf
D(1:length(x)+1:length(x)*length(x)) = inf;
indices = zeros(size(x));
% xones = cellfun(@(L) L(1),x);
% xends = cellfun(@(L) L(end),x);

indices(1) = find(x == startPoint(1) & y == startPoint(2));

% indices(1) = find(min(x-startPoint(1)) & min(y-startPoint(2)));
% if strcmp(sortDir,'+')
%     [~,indices(1)] = min(ii);
% elseif strcmp(sortDir,'-')
%     [~,indices(1)] = max(ii);    
% end
% if strcmp(sortDir,'+y')
%     [~,indices(1)] = min(y);
% elseif strcmp(sortDir,'+x')
%     [~,indices(1)] = min(x);
% elseif strcmp(sortDir,'-y')
%     [~,indices(1)] = max(y);
% elseif strcmp(sortDir,'-x')
%     [~,indices(1)] = max(x);
% end

D(:,indices(1)) = inf;
count = 1;

while count < length(x)
    % find min distance to next pixel
    [minD,ind] = min(D(indices(count),:));
    if minD > dThresh
        break;
    end
    indices(count+1) = ind;
    
   % exclude matches from being considered in the future
    D(:,indices(count+1)) = inf;
    count = count+1;
end

indices = indices(indices > 0);
xSort = x(indices);
ySort = y(indices);



if isempty(xSort) || isnan(sum(xSort))
    keyboard
end