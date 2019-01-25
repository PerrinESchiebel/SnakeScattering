
function [psx,psy] = MMS_bspline_skel(x,y,numsplinepts,square)

[numframes] =length(x);
if nargin == 2
    numsplinepts = 100;
end

psx = nan(numsplinepts, numframes);
psy = nan(numsplinepts, numframes);

currentDirectory = pwd;

addpath([currentDirectory '\GUI_V4\']);
addpath([currentDirectory '\GUI_V4\helper_functions\'])

track_Params.b_Spline.pred_Limit    = .01;
track_Params.b_Spline.reg_Term      = 1e-3;
track_Params.b_Spline.num_Basis     = 10;
track_Params.b_Spline.order         = 3;
track_Params.b_Spline.res           = 11;
track_Params.b_Spline.extra         = 0;
track_Params.b_Spline.type          = 'a-periodic-approx';
gui_States.b_Spline_Forward_Obj = b_Splines_Class(track_Params.b_Spline, 0);
% switch format
%     case 'Data'
%         for jj = 1:numframes
%             % b.states = fliplr(Data.Forward.pose(logical((Data.type==0).*(~isnan(Data.Forward.pose(:,1,jj)))'), :, (Data.frame_Index == jj))');
%             b.states = fliplr([x(:,jj),y(:,jj)]');
%             gui_States.b_Spline_Forward_Obj.process_Input(b);
%             vec = gui_States.b_Spline_Forward_Obj.eval_X(linspace(0,1,numsplinepts))*gui_States.b_Spline_Forward_Obj.basis_Coeff';
%         %         plot(x(:,jj),y(:,jj),'o');hold on;
%         %         plot(vec(:,1),vec(:,2));axis equal tight;drawnow;pause(0.1);hold off;
%             psx(:,jj) = vec(:,1);
%             psy(:,jj) = vec(:,2);
%         end

%     case 'other'
for jj = 1:numframes
    % b.states = fliplr(Data.Forward.pose(logical((Data.type==0).*(~isnan(Data.Forward.pose(:,1,jj)))'), :, (Data.frame_Index == jj))');
    b.states = [x{jj},y{jj}]';
    
    
    if isempty(b.states)
    else
        if nargin == 4
            b.states(:,b.states(1,:)<square(1)) = [];
            b.states(:,b.states(1,:)>square(2)) = [];
            b.states(:,b.states(2,:)<square(3)) = [];
            b.states(:,b.states(2,:)>square(4)) = [];
        end
        gui_States.b_Spline_Forward_Obj.process_Input(b);
        vec = gui_States.b_Spline_Forward_Obj.eval_X(linspace(0,1,numsplinepts))*gui_States.b_Spline_Forward_Obj.basis_Coeff';
        %         plot(x(:,jj),y(:,jj),'o');hold on;
        %         plot(vec(:,1),vec(:,2));axis equal tight;drawnow;pause(0.1);hold off;
        psx(:,jj) = vec(:,1);
        psy(:,jj) = vec(:,2);
    end
end
% end