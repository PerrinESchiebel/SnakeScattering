function [V,lambda,amp] = PCAsnakes(AllCurve,plotyesno)
if nargin == 1
    plotyesno = 1;
end
y = cell2mat(AllCurve');
y = y(:,isfinite(y(1,:)))';
[nt,ns] = size(y);
if ns > nt
    error('Matrix may not be oriented correctly');
end
C = cov(y);
[V,lambda]=eig(C);  %%%%V is matrix columns of which are ortho eigenvectors
lambda = diag(lambda);
amp = y*V;
if plotyesno == 1
    figure,plot(linspace((100-ns)/(2*100),1-(100-ns)/(2*100),ns),rad2deg(V(:,end)),'--','LineWidth',4,'Color',[95,95,93]./255);hold on;
    plot(linspace((100-ns)/(2*100),1-(100-ns)/(2*100),ns),rad2deg(V(:,end-1)),'LineWidth',4,'Color',[145,144,146]./255);hold on;
    %     legend('Mode 1','Mode 2');set(gca,'FontSize',36,'FontWeight','bold','LineWidth',4);
    xlabel('Fraction of arclength');ylabel('Degrees');
    
    figure,plot(flipud(lambda),'o','MarkerSize',12,'MarkerFaceColor',[95,95,93]./255);hold on;
    xlabel('Mode');ylabel('Eigenvalue');set(gca,'FontSize',36,'FontWeight','bold','LineWidth',4);
    sigsum = nan(ns,1);
    for ii=1:ns
        sigsum(ii) = sum(lambda(end-ii+1:end));
    end
    figure,line([-2 20],[1 1],'Color',[0.3,0.3,0.3],'LineWidth',4,'LineStyle','--');hold on;
    plot(sigsum/sum(lambda),'o','MarkerSize',12,'MarkerFaceColor',[95,95,93]./255);
    xlabel('Mode');ylabel('\sigma_K^2');set(gca,'FontSize',36,'FontWeight','bold','LineWidth',4);
    xticklabs = cell(20,1);xticklabs{1} = '1';xticklabs{2} = '2';xticklabs{10} = '10';xticklabs{20} = '20';
    set(gca,'YTick',[0 1],'XTick',1:20,'XTickLabel',xticklabs);
    ylim([0 1.3]);xlim([-2 20]);
    figure
    npts = size([AllCurve{1}],1);
    [hh,xx] = hist3([amp(:,end),amp(:,end-1)],[npts,npts]);  %%%%%%%%%
    pcolor(xx{1},xx{2},hh);shading flat,
    %     load blackbodyColormap,colormap(cmap),colorbar
    set(gca,'FontSize',36,'FontWeight','bold','LineWidth',4);
    xlabel('\alpha_1 (currently: arb units)');ylabel('\alpha_2');
    
end
end