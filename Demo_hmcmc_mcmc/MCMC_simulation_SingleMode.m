clear

updatePlotFlag = 1;     % 0/1 update plot during MCMC

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
epsilon = 5;
Sigma_prop = epsilon*eye(length(xp_init));


%% HMCMC
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
    
    % calculate log ratio (Metropolis)
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

%% analyze result
fprintf('Accepted rate is %f%% (%d out of %d)\n',counter/nmcmc*100,counter,nmcmc);
figure
plot(param_mcmc);
legend({'x1','x2'},'Location','best');
xlabel('iteration');

% --- burn-in
nburnin=0;
param_postburnin = param_mcmc(nburnin+1:end,:);

% --- mean and variance of parameters
E_param = mean(param_postburnin);
V_param = var(param_postburnin);

% --- MCMC standard error by effective sample size
[ESS,Sigma] = multiESS(param_postburnin);    % book assume no correlation between thetas
ESS = round(ESS);
% MCMCSE = sqrt(diag(cov(param_mcmc))'/ESS);    % useing covariance of mcmc
MCMCSE = sqrt(diag(Sigma)'/ESS);    % useing covariance from mutliESS

% --- display result
array2table([E_param; sqrt(V_param); MCMCSE],'VariableNames',{'x1','x2',},'RowNames',{'Mean','SE','MCMCSE'})
fprintf('Effective sample size is %d\n',ESS);


%% plot
figure
subplot(9,1,1:5);
hold on; box on;
[~,hc] = contour(X,Y,Z,5,'b','LineWidth',1);
plot(xp_init(1),xp_init(2),'o');
plot([xp_init(1);param_mcmc(:,1)],[xp_init(2);param_mcmc(:,2)],'r-')
xlabel('x1');
ylabel('x2');
uistack(hc,'top');  
subplot(9,1,7:9);
plot(param_mcmc);
xlabel('iteration');
legend({'x1','x2'},'Location','best');
