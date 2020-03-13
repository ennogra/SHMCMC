% script to run Subspace HMCMC

clear

global fitFlag x2
global lbub

%% some flags ans settings
% --- about data
weightFlag = 1;       % 0/1 use weight for intensity

% --- about HMCMC
updatePlotFlag = 1;    % 0/1 update plot during MCMC
dualAvgFlag = 0;        % 0/1 use dual averaging
nutsFlag = 1;         % 0/1 use NUTS (must use unit Mass and cannot update Mass)

%% Load data
s = load('mrgSAXS_1percent_150nm.dat');

% reduce # of data points in log space
% [~,fitIndex] = ismember(...
%     interp1(s(:,1),s(:,1), logspace(log10(s(1,1)),log10(s(end,1)),133),'nearest'),...
%     s(:,1));
fitIndex = 1:5:size(s,1);
qdata =     s(fitIndex,1);
Idata = s(fitIndex,2);
Idata_err = s(fitIndex,3);

% normalize 
Idata_max = max(Idata);
Idata = Idata/Idata_max;
Idata_err = Idata_err/Idata_max;

if weightFlag == 0
    weight = ones(size(qdata));
elseif weightFlag == 1
    weight = Idata./Idata_err;
    weight = weight/sum(weight)*length(qdata);
end
% 
% figure
% errorbar(qdata,Idata,Idata_err,'o-');
% set(gca,'yscale','log');
% set(gca,'xscale','log');
% return
%% initialize and prepare parameters
Inorm = 1;
Ibk = 1e-5;
R = 600;
sigma_R = 20;
sigma_q = 1e-4;

% --- [start, lb, ub, fitflag]
fitparam = [...
    Inorm           0.1                 10                  1
    Ibk             0                   1e-3                1     
    R               300                 1000                1
    sigma_R         1                   100                 1
    sigma_q         0                   1e-3                1
    ];  
x       = fitparam(:,1);
lb      = fitparam(:,2);
ub      = fitparam(:,3);
fitFlag = fitparam(:,4);

% --- transform parameter so it become unbounded
% for any unbounded parameters, better rescale its initial value to 1. This
% needs to be done manually.
if any(lb-ub>=0)  % ub must be larger than lb
    error('incorrect bound');
end
if any(x<=lb | x>=ub)
    error('start value cannot be at or outside the bounds')
end
lbub = [lb,ub];
lbub(:,3) = NaN;  % Last column indicate the transformation method; 
lbub(isfinite(lb) & isfinite(ub),3)     = 1;    % (a, b)
lbub(isfinite(lb) & ~isfinite(ub),3)    = 2;    % (a, Inf)
lbub(~isfinite(lb) & isfinite(ub),3)    = 3;    % (-Inf,b);
lbub(~isfinite(lb) & ~isfinite(ub),3)   = 4;    % (-Inf,Inf);
lbub(fitFlag==0,3) = 5;                         % no transformation
xtrans = transformx(x,lbub);

% --- assign to-be and not-to-be fitted parameters
x1 = xtrans(fitFlag==1);
x2 = xtrans(fitFlag==0);

%% preview data
tic
I_cal = fcn_saxs_sphere_hmcmc(x1,qdata);     % preview initial conditions before fitting
toc

% --- plot
figure
hold on;
errorbar(qdata,Idata,Idata_err,'o');
plot(qdata,I_cal,'-','linewidth',1.5);
hold off; box on;
set(gca,'xscale','log');
set(gca,'yscale','log');

%% Initialize parameters 
% ---  initialize variance of residual to guessed values
Sigmavi = var(log10(I_cal) - log10(Idata));

% --- initial U and dU values
xp_init = [x1(:)',log(Sigmavi)];

%% Define potential, kinetic energies and their differentials for entire parameter space
sspFlag = true(1,length(xp_init));
XFix = [];

% --- define potential and differential
% U = @(X)potentialE(X,qdata,Idata,weight,sspFlag,XFix);
opts = admOptions('i', 1);  % only differentiate the 1st input argument (not param)
% dU = @(X)admDiffFor(@potentialE,1,X,qdata,Idata,weight,opts);    % Forward mode
% dU = @(X)admDiffVFor(@potentialE,1,X,qdata,Idata,weight,opts);     % Vectorized forward mode
% dU = @(X)admDiffRev(@potentialE,1,X,qdata,Idata,weight,opts);    % Reverse mode
dU = @(X)admDiffFD(@potentialE,1,X,qdata,Idata,weight,sspFlag,XFix,opts);     % Numerical finite difference method (not an auto diff method)


%% initialize dU
tic
[dPotE_init,PotE_init] = dU(xp_init);
fprintf('One differetial takes %.5f sec\n',toc);

%% --- determine epsilon values
increment_target = 0.5; % product of epsilon and dPotE for the momentum's increment
% --- auto determine subspaces
nssp = 2;
[~,~,sspBin] = histcounts(log10(abs(dPotE_init)+eps),nssp);
sspFlag = false(nssp,length(dPotE_init));
epsilon = nan(1,nssp);
for ii=1:nssp
    sspFlag(ii,sspBin==ii) = true;
    epsilon(ii) = increment_target/max(abs(dPotE_init(sspFlag(ii,:))));     % use the maximum abs(dPotE) to define the epsilon value
end

% % --- manual assign substpace
% nssp = 2;
% sspFlag = logical([
%     0   1   0   1   1   1
%     1   0   1   0   0   0
%     ]);
% epsilon = nan(1,nssp);
% for ii=1:nssp
%     epsilon(ii) = increment_target/max(abs(dPotE_init(sspFlag(ii,:))));     % use the maximum abs(dPotE) to define the epsilon value
% end

% --- force epsilon values
%epsilon = [0.052,0.0066];
epsilon

%% --- dual averaging settings
if dualAvgFlag == 1
    delta = 0.8;   % optimal acceptance probability
    gamma = 0.05;
    t0 = 10;
    kappa = 0.75;
    Madap = Inf;    % Inf is preferred especially in presence of multiple modes in the probablility distribution
    H_init = 0;
    Hm = ones(1,nssp)*H_init;
    logepsilon_init = log(10*epsilon);
    logepsilon_mean = log(epsilon);
end

%% --- NUTS or simple HMCMC settings
if nutsFlag == 1  
    delta_max = 500;        % stopping rule (equation 3), large enough so not to interfere
else
    nstep_target = 10;
    lambda_target = nstep_target*epsilon;    % target total number simulation length
end


%% start HMCMC
% --- assign initial values for each subspace
xp = xp_init;
xp_ssp      = cell(1,nssp);
dPotE_ssp   = cell(1,nssp);
PotE_ssp    = cell(1,nssp);
logp_xp_ssp = cell(1,nssp);
for issp=1:nssp
    xp_ssp{issp} = xp(sspFlag(issp,:));
    idU = @(X)admDiffFD(@potentialE,1,X,qdata,Idata,weight,sspFlag(issp,:),xp(~sspFlag(issp,:)),opts);     % Numerical finite difference method (not an auto diff method)
    [dPotE_ssp{issp},tmp_PotE] = idU(xp_ssp{issp});
    logp_xp_ssp{issp} = -tmp_PotE;
end

% --- start MCMC
nmcmc = 5000;
param_mcmc = nan(nmcmc,length(xp));
epsilon_mcmc = nan(nmcmc,nssp);
dPotE_mcmc = nan(nmcmc,length(xp));

param_mcmc0 = xp;
epsilon_mcmc0 = epsilon;
dPotE_mcmc0 = dPotE_init;

%%
warning('off','MATLAB:Figure:FigureSavedToMATFile');
tic
for ii=1:nmcmc
    
    fprintf('%d of %d iterations:\n', ii,nmcmc);
    for issp  = 1:nssp
        % --- define U, dU, K, dK
        iU  = @(X)potentialE(X,qdata,Idata,weight,sspFlag(issp,:),xp(~sspFlag(issp,:)));
        idU = @(X)admDiffFD(@potentialE,1,X,qdata,Idata,weight,sspFlag(issp,:),xp(~sspFlag(issp,:)),opts);     
        iM  = eye(nnz( sspFlag(issp,:) ));  % default
        iK  = @(p)p*inv(iM)*p'/2 + log(det(iM));    % Euclidean-Gaussian kinetic energy
        idK = @(p)p*inv(iM);
        
        % --- move one step in HMCMC
        trial_hmcmc = 0;
        done_hmcmc = 0;
        while done_hmcmc == 0  % to prevent ill random pp values
            try
                if nutsFlag == 0    % use simple HMCMC
                    [xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp},logr(issp)] = hmcmcwalk(xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp},iU,iK,idU,idK,iM,epsilon(issp),lambda_target(issp));
             %       [xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp},logr(issp)] = hmcmcwalk(xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp},iU,iK,idU,idK,iM,epsilon(issp),nstep_target*epsilon(issp));                    
                elseif nutsFlag == 1 && dualAvgFlag == 1    % use NUTS with dual averaging for epsilon
                    [xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp},alpha(issp),nalpha(issp)] = nuts(xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp},iU,iK,idU,idK,iM,epsilon(issp),delta_max);
                elseif nutsFlag == 1 && dualAvgFlag == 0    % use NUTS without dual averaging
                    [xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp}] = nuts_noDA(xp_ssp{issp},logp_xp_ssp{issp},dPotE_ssp{issp},iU,iK,idU,idK,iM,epsilon(issp),delta_max);
                end
                done_hmcmc = 1;
                fprintf('      %d of %d subspace (%.1f sec) with epsilon=%.5f, xp: ',issp,nssp,toc,epsilon(issp));
                fprintf(['[',repmat('%.10f ',1,length(xp_ssp{issp})),'\b]\n'],xp_ssp{issp});
                fprintf(['                                                      dPotE: [',repmat('%.10f ',1,length(dPotE_ssp{issp})),'\b]\n'],dPotE_ssp{issp});                                
                
            catch
                trial_hmcmc = trial_hmcmc + 1;
                fprintf('For the %dth iteration: re-try %d time(s)\n', ii,trial_hmcmc);
            end
        end
        
        % --- assemble xp
        xp = nan(1,length(xp_init));
        dPotE = xp;
        for jssp = 1:nssp
            xp(sspFlag(jssp,:)) = xp_ssp{jssp};
            dPotE(sspFlag(jssp,:)) = dPotE_ssp{jssp};
        end        
    end
    
    % --- adapt epsilon with dual averaging
    if dualAvgFlag == 1
        for issp = 1:nssp
            if ii<=Madap
                if nutsFlag == 1    % update Hm for nuts
                    Hm(issp) = (1-1/(ii+t0))*Hm(issp) + 1/(ii+t0)*(delta-alpha(issp)/nalpha(issp));
                else                % update Hm for simple HMCMC
                    Hm(issp) = (1-1/(ii+t0))*Hm(issp) + 1/(ii+t0)*(delta-min(1,exp(logr(issp))));
                end
                logepsilon(issp) = logepsilon_init(issp) - sqrt(ii)/gamma*Hm(issp);
                logepsilon_mean(issp) = ii^(-kappa)*logepsilon(issp) +(1-ii^(-kappa))*logepsilon_mean(issp);  % moving average
                epsilon(issp) = exp(logepsilon(issp));
            else
                epsilon(issp) = exp(logepsilon_mean(issp));
            end
        end
%     else
%         for issp = 1:nssp
%             epsilon(issp) = increment_target/max(abs(dPotE(sspFlag(issp,:))));     % use the maximum abs(dPotE) to define the epsilon value
%         end
    end
    
    % --- collect result
    param_mcmc(ii,:) = xp;
    epsilon_mcmc(ii,:) = epsilon;    
    dPotE_mcmc(ii,:) = dPotE;
    
    % --- update plot
    if updatePlotFlag == 1
        I_cal = fcn_saxs_sphere_hmcmc(xp(1:end-1),qdata);     % preview initial conditions before fitting
        if ii==1
            hfig = figure;
            subplot(3,1,1:2);
            hold on;
            errorbar(qdata,Idata,Idata_err,'o');
            hline = plot(qdata,I_cal,'-','linewidth',1.5);
            hold off; box on;
            set(gca,'xscale','log');
            set(gca,'yscale','log');
            subplot(3,1,3);
            hline2 = plot(epsilon_mcmc);
            set(gca,'xlim',[0,nmcmc]);
            xlabel('iteration');
            ylabel('\epsilon');
            set(gca,'yscale','log');
            htitle = sgtitle(sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc));            
        end
        hline.YData = I_cal;
        htitle.String = sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc);
        for ihline2 = 1:length(hline2)
            set(hline2(ihline2),'XData',1:nmcmc,'YData',epsilon_mcmc(:,ihline2));
        end
        pause(0.01);        
    end
    
    save tmp.mat        % iteratively save result
end
total_time =toc;


% --- save result
if nutsFlag == 1 && dualAvgFlag == 1
    save fresult_subspace_NUTS_DA_v2.mat
elseif nutsFlag == 1 && dualAvgFlag == 0
    save fresult_subspace_NUTS_NoDA_v2.mat    
elseif nutsFlag == 0 && dualAvgFlag == 1    
    save fresult_subspace_Simple_DA_v2.mat        
elseif nutsFlag == 0 && dualAvgFlag == 0
    save fresult_subspace_Simple_NoDA_v2.mat            
end
