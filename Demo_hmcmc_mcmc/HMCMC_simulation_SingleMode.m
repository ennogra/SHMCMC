clear

updatePlotFlag = 1;     % 0/1 update plot during MCMC
dualAvgFlag = 0;        % 0/1 use dual averaging
nutsFlag = 0;           % 0/1 use NUTS
updateMassFlag = 0;     % 0/1. better not update mass in presence of multiple modes

%% - target distribution (multivariate normal)
mu = [5,8];
Sigma = [2 0.8; 0.8 .9];
% Sigma = [1 0.99; 0.99 1];


%% define U, dU, K, dK
% --- potential energy (directly define)
U1 = @(X)(X-mu)*inv(Sigma)*(X-mu)'/2 - log(1/sqrt((2*pi)^size(X,2)*det(Sigma)));
% % dU = @(x)(x-mu)*inv(Sigma);
% dU1 = @(X)dPotE(X,mu,Sigma,U1);
 dU1 = @(X)deal((X-mu)*inv(Sigma),U1(X));
% % dU = @(x)gradest(U,x);

% --- potential energy (define using external function)
param.mu = mu;
param.Sigma = Sigma;
U2 = @(X)potentialE_SingleMode(X,param);
opts = admOptions('i', 1);  % only differentiate the 1st input argument (not param)
% dU2 = @(X)admDiffFor(@potentialE_SingleMode,1,X,param,opts);    % Forward mode
%dU2 = @(X)admDiffVFor(@potentialE_SingleMode,1,X,param,opts);   % Vectorized forward mode
% dU2 = @(X)admDiffRev(@potentialE_SingleMode,1,X,param,opts);    % Reverse mode
dU2 = @(X)admDiffFD(@potentialE_SingleMode,1,X,param,opts);     % Numerical finite difference method (not an auto diff method)

% --- choose format of potential energy
% - use analytical expression
U = U1;
dU = dU1;

% - use AdiMAT
% U = U2;
% dU = dU2;

% --- kninetic energy
%M = inv(diag(diag(Sigma)));     %  mass inverse can be 
M = eye(2);
K = @(p)p*inv(M)*p'/2+log(det(M));
dK = @(p)p*inv(M);

%% ground truth
% x_true = mvnrnd(mu,Sigma,2000);     % sample from true distribution 
xlist = linspace(mu(1)-4,mu(1)+4,190);
ylist = linspace(mu(2)-4,mu(2)+4,200);
[X,Y] = meshgrid(xlist,ylist);
Z = reshape(mvnpdf([X(:),Y(:)],mu,Sigma),size(X));
%Z= exp(-reshape(U([X(:),Y(:)]),size(X)));

%% start values
xp_init = [2,12];
[dPotE_init,PotE_init] = dU(xp_init);
logp_xp_init = -PotE_init;

%% find good initial value for epsilon
pp_init = mvnrnd(zeros(1,length(xp_init)),M);
logp_pp_init = -K(pp_init);
logp_init = logp_xp_init + logp_pp_init;
epsilon_init = 1;   % large start value is preferred than small start value
[epsilon,counter] = heuristic4epsilon(xp_init,pp_init,dPotE_init,epsilon_init,U,K,dU,dK,logp_init);
fprintf('Inialized epsilon=%f, after %d iterations of heuristic\n',epsilon,counter);

% --- heuristic method works the best when only one mode exists.
% if it got stuck duing dual averging, do not use the heuristic method
% for the starting value; instead, specify a small value of inital epsilon
epsilon = .3;  % force epsilon to be this value

%% dual averaging settings
if dualAvgFlag == 1
    delta = 0.65;   % optimal acceptance probability
    gamma = 0.15;
    t0 = 10;
    kappa = 0.55;
    Madap = Inf;    % Inf is preferred especially in presence of multiple modes in the probablility distribution
    H_init = 0;
    Hm = H_init;
    logepsilon_init = log(10*epsilon);
    logepsilon_mean = log(epsilon);
end

%% NUTS or simple HMCMC settings
if nutsFlag == 1
    delta_max = 500;   % stopping rule (equation 3), large enough so not to interfere 
else
    lambda_target = .1;    % target total number simulation length
end

%% Mass update settings
if updateMassFlag == 1
    cov_start_idx = 10*length(xp_init);
end

%% HMCMC
% --- asign start values
xp = xp_init;
logp_xp = logp_xp_init;
dPotE = dPotE_init;

% --- start MCMC
nmcmc = 500;
param_mcmc = nan(nmcmc,2);
epsilon_nmcmc = nan(nmcmc,1);

tic
for ii=1:nmcmc

    % --- update mass
    if updateMassFlag == 1 && ii > cov_start_idx
        M = inv(diag(diag(nancov(param_mcmc(floor(cov_start_idx/2):ii-1,:)))));
     %   diag(M)'
    end    

    % --- move one step in HMCMC
    trial_hmcmc = 0;
    done_hmcmc = 0;
    while done_hmcmc == 0  % to prevent ill random pp values
        try
            if nutsFlag == 0    % use simple HMCMC
                [xp,logp_xp,dPotE,logr] = hmcmcwalk(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,lambda_target);                
            elseif nutsFlag == 1 && dualAvgFlag == 1    % use NUTS with dual averaging for epsilon
                [xp,logp_xp,dPotE,alpha,nalpha] = nuts(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,delta_max);
            elseif nutsFlag == 1 && dualAvgFlag == 0    % use NUTS without dual averaging
                [xp,logp_xp,dPotE] = nuts_noDA(xp,logp_xp,dPotE,U,K,dU,dK,M,epsilon,delta_max);
            end
            done_hmcmc = 1;
        catch
            trial_hmcmc = trial_hmcmc + 1;                        
            fprintf('For the %dth iteration: re-try %d time(s)\n', ii,trial_hmcmc);
        end
    end
    
    % --- adapt epsilon with dual averaging
    if dualAvgFlag == 1
        if ii<=Madap
            if nutsFlag == 1    % update Hm for nuts
                Hm = (1-1/(ii+t0))*Hm + 1/(ii+t0)*(delta-alpha/nalpha);
            else                % update Hm for simple HMCMC
                Hm = (1-1/(ii+t0))*Hm + 1/(ii+t0)*(delta-min(1,exp(logr)));
            end 
            logepsilon = logepsilon_init - sqrt(ii)/gamma*Hm;
            logepsilon_mean = ii^(-kappa)*logepsilon +(1-ii^(-kappa))*logepsilon_mean;  % moving average
            epsilon = exp(logepsilon);
        else
            epsilon = exp(logepsilon_mean);
        end
    end
    
    % --- collect
    param_mcmc(ii,:) = xp;
    epsilon_nmcmc(ii) = epsilon;    % for dual averaging
    
    % --- plot
    if updatePlotFlag == 1
        if ii==1
            hf = figure;
            ha1 = subplot(9,1,1:4);
            hold on; box on;
            [~,hc] = contour(X,Y,Z,10,'b','LineWidth',1);
            plot(xp_init(1),xp_init(2),'o');
            hplot1a = plot([xp_init(1);param_mcmc(:,1)],[xp_init(2);param_mcmc(:,2)],'r.-','linewidth',0.05);
            hplot1b = plot(param_mcmc(ii,1),param_mcmc(ii,2),'co','MarkerFaceColor','c');
            xlabel('x1');
            ylabel('x2');
                        
            ha2 = subplot(9,1,6:7);
            hplot2 = plot(param_mcmc);
            set(gca,'xlim',[1,nmcmc]);
            set(gca,'xticklabel',[]);
            legend({'x1','x2'},'Location','best');            
            
            ha3 = subplot(9,1,8:9);            
            hplot3 = plot(epsilon_nmcmc);
            set(gca,'yscale','log');
            ylabel('\epsilon');
            xlabel('iteration');
            set(gca,'xlim',[1,nmcmc]);

            htitle = sgtitle(sprintf('%d of %d',ii,nmcmc));
        else
            set(hplot1a,'XData',[xp_init(1);param_mcmc(:,1)],'YData',[xp_init(2);param_mcmc(:,2)]);
            set(hplot1b,'XData',param_mcmc(ii,1),'YData',param_mcmc(ii,2));
            uistack(hc,'top');
            hplot2(1).YData = param_mcmc(:,1);
            hplot2(2).YData = param_mcmc(:,2);            
            set(hplot3,'YData',epsilon_nmcmc);
            set(htitle,'String',sprintf('%d of %d',ii,nmcmc));                                    
        end
        pause(0.001);
    end
end
total_time = toc;
% if updatePlotFlag == 1
%     return;
% end

% save msplot_HMCMC_single

%% analyze result
figure
plot(param_mcmc);
legend({'x1','x2'},'Location','best');
xlabel('iteration');

% --- burn-in
nburnin=5;
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
subplot(9,1,1:4);
hold on; box on;
[~,hc] = contour(X,Y,Z,5,'b','LineWidth',1);
plot(xp_init(1),xp_init(2),'o');
plot([xp_init(1);param_mcmc(:,1)],[xp_init(2);param_mcmc(:,2)],'r-')
xlabel('x1');
ylabel('x2');
uistack(hc,'top');  

subplot(9,1,6:7);
plot(param_mcmc);
set(gca,'xlim',[1,nmcmc]);
set(gca,'xticklabel',[]);
legend({'x1','x2'},'Location','best');

subplot(9,1,8:9);
plot(epsilon_nmcmc);
set(gca,'yscale','log');
ylabel('\epsilon');
xlabel('iteration');


