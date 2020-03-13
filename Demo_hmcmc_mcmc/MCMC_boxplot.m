clear

updatePlotFlag = 0;     % 0/1 update plot during MCMC

%% - target distribution (multivariate normal)
mu = [5,8];
Sigma = [2 0.8; 0.8 .9];
% Sigma = [1 0.999; 0.999 1];
prob = @(X)mvnpdf(X,mu,Sigma);


%% ground truth
% x_true = mvnrnd(mu,Sigma,2000);     % sample from true distribution
xlist = linspace(mu(1)-4,mu(1)+4,190);
ylist = linspace(mu(2)-4,mu(2)+4,200);
[X,Y] = meshgrid(xlist,ylist);
Z= reshape(prob([X(:),Y(:)]),size(X));

%% start values
xp_init = [2,12];

%% proposal distribution (assume normal)
nrep = 100;
epsilon_list = repmat([0.05,0.1,0.5,1,5,10,50,100],nrep,1);
ESS_all = nan(size(epsilon_list));

for ie = 1:numel(epsilon_list)
    epsilon = epsilon_list(ie)
    Sigma_prop = epsilon*eye(length(xp_init));
    
    % --- asign start values
    xp = xp_init;
    
    % --- start MCMC
    nmcmc = 500;
    param_mcmc = nan(nmcmc,length(xp));
    counter = 0;
    
    tic
    for ii=1:nmcmc
        
        % propose a move from proposal joint-distribution
        xp_tmp = mvnrnd(xp,Sigma_prop);
        
        % calculate log ratio
        logr = log(mvnpdf(xp_tmp,mu,Sigma)) - log(mvnpdf(xp,mu,Sigma));
        
        % take it or not
        if log(rand)<logr
            xp = xp_tmp;
            counter = counter+1;
        end
        
        % --- collect
        param_mcmc(ii,:) = xp;
        
        % --- plot
        if updatePlotFlag == 1
            if ii==1
                hf = figure;
                subplot(9,1,1:5);
                hold on; box on;
                [~,hc] = contour(X,Y,Z,10,'b','LineWidth',1);
                plot(xp_init(1),xp_init(2),'o');
                hplot1a = plot([xp_init(1);param_mcmc(:,1)],[xp_init(2);param_mcmc(:,2)],'r.-','linewidth',0.05);
                hplot1b = plot(param_mcmc(ii,1),param_mcmc(ii,2),'co','MarkerFaceColor','c');
                xlabel('x1');
                ylabel('x2');
                subplot(9,1,7:9);
                hold on; box on;
                hplot2 = plot(param_mcmc);
                set(gca,'xlim',[1,nmcmc]);
                xlabel('iteration');
                legend({'x1','x2'},'Location','best');
                htitle = sgtitle(sprintf('%d of %d',ii,nmcmc));
            else
                set(hplot1a,'XData',[xp_init(1);param_mcmc(:,1)],'YData',[xp_init(2);param_mcmc(:,2)]);
                set(hplot1b,'XData',param_mcmc(ii,1),'YData',param_mcmc(ii,2));
                set(htitle,'String',sprintf('%d of %d',ii,nmcmc));
                hplot2(1).YData = param_mcmc(:,1);
                hplot2(2).YData = param_mcmc(:,2);
                uistack(hc,'top');
            end
            pause(0.001);
        end
    end
    total_time = toc;
    
    % if updatePlotFlag == 1
    %     return;
    % end
    % --- burn-in
    if epsilon<=0.2
        nburnin=30;
    elseif epsilon<1
        nburn = 10;
    else
        nburnin = 0;
    end
    nburnin = 0;
    param_postburnin = param_mcmc(nburnin+1:end,:);
    
    % --- MCMC standard error by effective sample size
    [ESS,~] = multiESS(param_postburnin);    % book assume no correlation between thetas
    ESS = round(ESS);
    
    ESS_all(ie) = ESS;
    
end
return

% save msplot_MCMC_singlemode_boxplot.mat;


%% plot
load('msplot_MCMC_singlemode_boxplot.mat');
figure
boxplot(ESS_all,epsilon_list(1,:),'Notch','on','OutlierSize',.1);
set(gca,'XMinorTick','off','YMinorTick','on');
set(gca,'TickLength',[0.02,0.02]);
set(gca,'ylim',[20,110]);
%set(gca,'xtick',[1:2:size(epsilon_list,2)],'XTickLabel',cellstr(num2str(epsilon_list(1,1:2:end)')));

xlabel('\Delta');
ylabel('Effective sample size');