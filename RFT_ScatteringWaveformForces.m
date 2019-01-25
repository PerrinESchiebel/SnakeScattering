%%last edited 06/15/2017
%%fixed the radius 05/08/17
addpath F:\Dropbox\Research\RFT\surface_RFT\calc_jointPower
tic;
dir = pwd;
ss = strfind(dir,'Dropbox');
path = dir(1:ss-1);
addpath([path,'Dropbox\Research\RFT\surface_RFT'])
z = 8;
load(strcat(path,'Dropbox\Research\RFT\surface_RFT\','f',num2str(z),'para.mat'));
L = 0.373;
% numsegs = 100;
dsPost = 0.0254*0.25/2;
numsegs = round(L/dsPost);
numtimes = 70;
mu = 0.1;
snakemass = 0.019;
% snakemass = 0.005;
bellydrag = snakemass*9.81*mu;
% snake = [-1,0.390,0.0035];  %%%%[wavespeed, length, radius]

snake = [-1,0.390,0.004];  %%%%[wavespeed, length, radius]
% L = 0.390;

frequency = 4;
Po = 0.003;

liftfact = 0.8;
liftfact = 1;
muAl = 0.2;  %%%friction coefficient between Al and glass particles
fFac = [1,mu/muAl];   %%% [scaling factor for Fperp, scaling factor for Fpara];
                                                                                         %%FcnstW = scaling factor for "walls" of snake body
Fcnst = (z*snake(2)/numsegs)/(z*0.03);  %%%Because depth dep. is tied into the fit functions, do not need it explicitly. 
                                                                                            %%%drag plate width = 0.03cm
                                                                                              %%%%FcnstV = scaling factor for bottom of snake
                                                                                              %%%if no belly drag set to zero
                                                                                              

xis = 2;
kmax = 25.2;
kplds = kmax*L/xis;

%%%%%%%%INITIALIZE MATRICES
nxis = length(xis);
nkplds = length(kplds);

vels = nan(nxis,nkplds);
torques = nan(nxis,nkplds);
allexits = nan(numtimes,nkplds,nxis);
vwmatall = nan(numtimes,3,nkplds);
powers = nan(nxis,nkplds);
internalTauksis = nan(numtimes,numsegs,nxis);
internalTaukap = nan(numtimes,numsegs,nkplds);
Pjoint = nan(numtimes,numsegs,nxis,nkplds);
Pdiss = nan(numtimes,numsegs,nxis,nkplds);
slips = nan(nxis,nkplds);
fforPo = nan(nxis,nkplds);
allPjoint = nan(nxis,nkplds);
sumPint = nan(nxis,nkplds);
sumPext = nan(nxis,nkplds);
freqs = nan(nxis,nkplds);
Fx = nan(numsegs,numtimes,nxis,nkplds);
Fy = nan(numsegs,numtimes,nxis,nkplds);
%%%%FOR SOLVER
vnot = [0.1,0,0];
opts = optimset('Diagnostics','off', 'Display','Off');
AllForces = nan(numsegs,numtimes,3);
for ll = 1:nxis
    xi = xis(ll);
    curvetau = nan(numtimes,length(kplds));
%     frequency = abs((xi*snake(1))/snake(2));
    for jj = 1:nkplds
        thm = kplds(jj)/(2*pi);       
        [xCoG,yCoG,vxCoG,vyCoG,th,dzetadt,lift] = CoGPos2(thm,xi,L,frequency,numsegs,numtimes);
        vwmat = nan(numtimes,3);
        taus = nan(numtimes,1);
        Powert = nan(numtimes,1);
        slip_angle = nan(numtimes,1);
        for ii=1:numtimes-1           
            x = xCoG(ii,:)';
            y = yCoG(ii,:)';
            theta = th(ii,:)';
            dzdt = dzetadt(ii,:)';
            liftind = find(abs(lift(ii,:))>liftfact);
            vCoG = [vxCoG(ii,:)',vyCoG(ii,:)'];
            searchfunc2D = @(V) sum(RFT_reportJointPower(V,x,y,theta,dzdt,liftind,bellydrag,vCoG,Fcnst,fFac,z,f), 1); %%%square snake!            
            [vout, ~, exitflag] = fsolve(searchfunc2D, vnot, opts);
            allexits(ii,jj,ll) = exitflag;
            if exitflag == 1
                if abs(vout(1))>5
                else
                    vwmat(ii,:) = [mean(vCoG) + vout(1:2),vout(3)];
                    [FandT, cosT, Tint, Pext,P] = RFT_reportJointPower(vout,x,y,theta,dzdt,liftind,bellydrag,vCoG,Fcnst,fFac,z,f);   %%square snake!
                    vnot = vout;
                    internalTaukap(ii,:,jj) = Tint;
                    internalTauksis(ii,:,ll) = Tint;
                    taus(ii) = max(Tint);
                    Powert(ii) = sum(abs(Pext));
                    Pjoint(ii,:,ll,jj) = P;
                    Pdiss(ii,:,ll,jj) = Pext;
                    slip_angle(ii) = mean(acosd(abs(cosT)));
                    AllForces(:,ii,:) = FandT;
%                     allslips(:,ii) = acosd(abs(cosT));
                end
            end
            curvetau(ii,jj) = max(taus);
        end
        allPjoint(ll,jj) = max(squeeze(max(Pjoint(:,:,ll,jj))));
        sumPint(ll,jj) = sum(sum(Pjoint(:,:,ll,jj),'omitnan'));
        sumPext(ll,jj) = sum(sum(Pdiss(:,:,ll,jj),'omitnan'));
        PowerhereNoF = max(squeeze(max(Pjoint(:,:,ll,jj))))/frequency;
        fforPo(ll,jj) = Po/PowerhereNoF;
%     allPjoint(ll,jj) = max(squeeze(max(Pjoint(:,:,ll,jj))));
%         sumPint(ll,jj) = sum(sum(Pjoint(:,:,ll,jj),'omitnan'));
%         sumPext(ll,jj) = sum(sum(Pdiss(:,:,ll,jj),'omitnan'));
%         PowerhereNoF = max(squeeze(max(Pdiss(:,:,ll,jj))))/frequency;
%         fforPo(ll,jj) = Po/PowerhereNoF;
%          PowerhereNoF = intP(ll,jj)/frequency;
%         fforintPo(ll,jj) = intPo/PowerhereNoF;
        vwmatall(:,:,ll) = vwmat;
        vels(ll,jj) = mean(vwmat(isnan(vwmat(:,1))==0,1));
        torques(ll,jj) = (max(curvetau(:,jj)));
        powers(ll,jj) = sum(Powert(isnan(Powert)==0));
        slips(ll,jj) = mean(slip_angle(isnan(slip_angle)==0));
        freqs(ll,jj) = frequency;    
        Fx(:,:,ll,jj) = AllForces(:,:,1);
        Fy(:,:,ll,jj) = AllForces(:,:,2);
    end
end
vscaled = vels.*fforPo./freqs;
% vscaled2 = vels.*fforintPo./freqs;
Fmag = sqrt(AllForces(:,:,1).^2+AllForces(:,:,2).^2);
toc
%%
thresh = 0.0037; % standard deviation of forces taken from a peg at random which wasn't contacted=0.0037;
load('F:\Dropbox (GaTech)\Research\scattering\Data\ForceMatsOnlyTouched\AllForcesTouch_hand_trimmed_ringing.mat')
nforcetrials = length(SnakeForces);
mag1 = cell(1,nforcetrials);
mag2 = cell(1,nforcetrials);
for ii=1:length(SnakeForces)
    m1 = 0;m2 = 0;
p1 = SnakeForces(ii).Peg1;
p2 = SnakeForces(ii).Peg2;
if isempty(p1) == 0
m1 = sqrt(p1(1,:).^2+p1(2,:).^2);
m1 =  m1(m1>=thresh);
mag1{ii} = m1;
end
if isempty(p2)==0
m2= sqrt(p2(1,:).^2+p2(2,:).^2);
m2 = m2(m2>=thresh);
mag2{ii} = m2;
end
end
% histogram([cell2mat(mag1),cell2mat(mag2)],'Normalization','pdf');
%%
figure

Fposts = [cell2mat(mag1),cell2mat(mag2)].*1000;
FRFT = reshape(Fmag.*1000,1,numel(Fmag));

[vals,bins] = histcounts(FRFT,'Normalization','pdf');
plot([bins(1),bins(1:end)+(bins(2)-bins(1))/2],[0,vals,0],'LineWidth',3,'Color','b')
hold on;
bins = 0:2:60;
[vals,bins] = histcounts(Fposts,bins,'Normalization','pdf');
plot([bins(1),bins(1:end)+(bins(2)-bins(1))/2],[0,vals,0],'LineWidth',3,'Color','r')
legend('Force magnitudes from posts','Force magnitudes from RFT')
xlabel('Force (mN)')
ylabel('Probability Density');
set(gca,'FontSize',24,'LineWidth',2)
intruderWidth = 30; %mm
dsPost = 25.4*0.25; %mm
    % SandAvg = 179/intruderWidth*dsPost;
CarpetAvg = 880/intruderWidth*dsPost;
SandAvg = 179;
CarpetAvg = 880;



ratMedians = median(Fposts,'omitnan')/median(FRFT,'omitnan');
RatAvgDrag = CarpetAvg/SandAvg;
maxRat = 1550/SandAvg;
minRat = 440/SandAvg;
mPost = median(Fposts,'omitnan');

sigmaPost = std(Fposts,'omitnan');


annotation('textbox',[0.5 0.5, 0.5 0.5],'String',strcat('median of Fpost over median of Fsand = ',num2str(ratMedians)),'LineStyle','none','FontSize',18)
annotation('textbox',[0.6 0.5, 0.6 0.5],'String',strcat('mean of average drag ratio = ',num2str(RatAvgDrag)),'LineStyle','none','FontSize',18)
annotation('textbox',[0.6 0.5, 0.6 0.5],'String',strcat('range of average drag ratio = ',num2str(minRat),' to ',num2str(maxRat)),'LineStyle','none','FontSize',18)