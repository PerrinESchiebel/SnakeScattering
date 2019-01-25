function [splinedSnakesX,splinedSnakesY] = trackAllShagVideosTFB(directory)

if nargin<1
   directory =  'E:\scattering\ShagTalbot\123016\vids\';
end

splineDir = 'E:\scattering\ShagTalbot\matFiles\';
% splineDir = 'E:\scattering\ShagTalbot\matFiles1115\';

files = dir([directory '*.avi']);

splinedSnakesX = cell(size(files));
splinedSnakesY = cell(size(files));
pegCenters = cell(size(files));
pegRadii = cell(size(files));

% mkdir('splinedSnakesOnShagNew\');

% num = length(dir(splineDir));
ff = dir(splineDir);
num = str2num(ff(end).name(1:end-4));

if isempty(num)
    num = 0;
end 

% % % onlyDoThisToday = true;
% % % if onlyDoThisToday
% % %     if strcmp(directory,'E:\scattering\ShagTalbot\111416\vids\')
% % %         startIndex = 24;
% % %     else 
% % %         startIndex = 1;
% % %     end
% % % else
% % %         startIndex = 1;
% % % end


startIndex = 1;
for i = startIndex:length(files)
    path = [directory files(i).name];
    savenum = num+i-startIndex+1;    %minus 2 if using numfiles
    if startIndex > 1
        savenum = savenum;   %skip broken one   if commented: continue including broken one
    end
         disp([num2str(i) ' / ' num2str(length(files)) ' (num ' num2str(savenum) ')']);

    [snakeX,snakeY,splineX,splineY,pegXY,pegR,keep] = trackSnakeInVideoOnShagNewTFB(path); %#ok<*ASGLU>
%     save(['splinedSnakesOnShagNew\' num2str(i,'%03i') '.mat'],'path','pegXY','pegR','snakeX','snakeY','splineX','splineY','keep')
   save([splineDir num2str(savenum,'%03i') '.mat'],'path','pegXY','pegR','snakeX','snakeY','splineX','splineY','keep')
   
    splinedSnakesX{i} = splineX;
    splinedSnakesY{i} = splineY;
    pegCenters{i} = pegXY;
    pegRadii{i} = pegR;
    
end
disp(['Done - ' directory(end-12:end-5)]);