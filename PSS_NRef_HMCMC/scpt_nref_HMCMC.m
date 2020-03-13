
% lines or blockes commented with @@@@@ can be changed

clear all

%% define global parameters 
global fitFlag x2 
global lbub
global nsets   % data sets for multiple curve fittings
global lambda  conv_method % nuetron wavelength
global interface_profile plotEDPFlag       
global binedp_delta_cutoff edp_0 edp
global nSlice nlayer

%% some flags ans settings
% --- about beam
lambda = 4.75;              % neutron wavelength

% --- about density profile
binedp_delta_cutoff = 1e-12;
interface_profile   = 2;        % 0/1/2: error / hyperbolic tangent/cumulative skewed norm
nSlice = 2000;              % number of slices to slice the entire film. 
            % The larger the slower; keep this value smaller than the total
            % film thickness, i.e. use ~1A as the smallest slice thickness.
plotEDPFlag = 0;                % 0/1: plot edp in the call function
            
% --- about data
weightFlag = 1;       % 0/1 use weight for intensity
conv_method = 1;      % Only takes 1 for now. 1/2/3: vectorization/forloop/conv: 3 is the fastest (only 1 is implemented for now)


% --- about HMCMC
updatePlotFlag = 1;    % 0/1 update plot during MCMC
dualAvgFlag = 0;        % 0/1 use dual averaging
nutsFlag = 1;         % 0/1 use NUTS (must use unit Mass and cannot update Mass)

% --- plot settings for multiple curves
plotshift = 0.1;    % shift factor for multiple data sets when plotting


%% load data
% @@@@@ --- load data; you may need to change data format to columns
fpath = pwd;
flist = {
   'PSS25nm_I6mM_NaY9_RAW.dat'; 
    };

nsets = length(flist); 
qset = cell(1,nsets);  % collect data sets
Iset = cell(1,nsets); 
Ierr_set = cell(1,nsets); 
Iweight_set = cell(1,nsets);

% --- points to delete
ind_bad = {
     []         % uncomment if no bad point to remove for this set
%     [30 34 49 54]
   };
% ind_bad = cell(length(flist),1);

for iset=1:nsets
    % --- for DAT format
    if strcmpi(flist{iset}(end-2:end),'ABS')
        s = importdata(fullfile(fpath,flist{iset}),' ',5);
        s = s.data;
    elseif strcmpi(flist{iset}(end-2:end),'DAT')
        s = load(fullfile(fpath,flist{iset}));     % for DAT format
    end

    % -- remove bad sigma q
    s(ind_bad{iset},:) = [];

    [~,ia,~] = unique(s(:,1),'rows','first');       % remove points of identical qz
    s = s(ia,:);
    ind = s(:,2)<=0 | s(:,1)>0.08;   % % remove negative I and too noise data
    s(ind,:) = [];
    
    % --- reduce data density
    fitRange = 1:1:size(s,1);
    s = s(fitRange,:);
    
    q00 = s(:,1);     % q values (A^-1); a list column
    I00 = s(:,2);        % intensity
    I00_err = s(:,3);       % weight
    
    % --- normalize intensity
    I00_max = max(I00);
    I00 = I00/I00_max;
    I00_err = I00_err/I00_max;
    
    % --- add weight
    if weightFlag == 0
        I00_weight = ones(length(q00),1);
    elseif weightFlag == 1
        I00_weight = I00./I00_err;
        I00_weight = I00_weight/sum(I00_weight)*length(q00);     % normalize so the sum is the number of points
    end
        
    % --- collect
    qset{iset} = q00;
    Iset{iset} = I00;
    Ierr_set{iset} = I00_err;
    Iweight_set{iset} = I00_weight;
end

% --- view raw data
figure
hold on; box on;
for iset=1:nsets
    errorbar(qset{iset},Iset{iset}.*qset{iset}.^4*plotshift^(iset-1),Ierr_set{iset}.*qset{iset}.^4*plotshift^(iset-1),'o-');
end
set(gca,'yscale','log');
legend(flist,'Location','best','Interpreter','none');

%% initialize paramters
% some constant SLDs
sld_Si = 2.07e-6;       % silicon SLD in unit 1/A^2
sld_D2O = 6.36e-6;      % D2O SLD
sld_SiOx = 3.475e-6;    % SiOx SLD
sld_Air = 0;            % Air SLD

Inorm = 1.06;
Ibk = 2.5e-9;

resdawave = 0.055;       % relative dispersion of wavelength
resconstant = 0.0;      % constant resoltuion for qz;

alphai_offset = 0.0032;       % angle offset

% --- sld of subphase and environment
sld_sub = sld_D2O;      % subphase SLD unit 1/A^2
sld_env = sld_Si;       % environment SLD unit 1/A^2 

% --- roughness on the subphase
sigma_sub = 93;         % it is film/D2O interface here

% --- give layer profile from top to bottom (excluding env and sub)
% [thickness, SLD, roughness of the layer-top interface (toward env), interface asymmetry factor alpha]
layeredp = [
          18.4193079599663                  sld_SiOx          3.48101757750399                         0    % determined from reference reflectivity 
         0.658126327922431      5.08604356617186e-06          12.7563701598632          12.1083736630674
          206.094858111351      1.19498225032736e-09          278.516742811243          11.3491353717867
          47.6067164820199      3.34848579245317e-06           488.22952715644          9.83595599482467
           183.54682191343      1.48426800033844e-06          37.0933107952107          7.57328000969597
          127.833567042829      6.19857144198796e-11          51.0558108245798          1.82105294763546
          ];
nlayer = size(layeredp,1);  % number of layers
fitFlag_layeredp = [
    0           0           0           0
    1           1           1           1
    1           1           1           1 
    1           1           1           1
    1           1           1           1
    1           1           1           1
    ];
% lb_layeredp = eps*ones(nlayer,4);   % lower boundary 
lb_layeredp = [ ...
    0*ones(nlayer,1)      0*ones(nlayer,1)          0*ones(nlayer,1)      [0;0;0;0;-Inf;-Inf]
    ];     % lower boundary 
% ub_layeredp = Inf*ones(nlayer,4);     % upper boundary
ub_layeredp = [ ...
    Inf*ones(nlayer,1)      Inf*ones(nlayer,1)          Inf*ones(nlayer,1)      [Inf;Inf;Inf;Inf;Inf;Inf]
    ];     % upper boundary

%% construct fitting parameters
% @@@@@ use 0/1 for the last coloumn for fixed/fitted paramters
fitparam = [...
    Inorm           0.8                 1.2                 1
    Ibk             0                   1e-4                1     
    resdawave       0                   Inf                 1
    resconstant     -1                  10000               0
    alphai_offset   -0.1                0.1                 1
    
	sld_env         1.9e-6              2.2e-6              0
    sld_sub         5.0e-6              7.0e-6              0
    sigma_sub       0                   Inf                 1
    
    layeredp(:)     lb_layeredp(:)      ub_layeredp(:)      fitFlag_layeredp(:)
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
% if any(x==lb | x==ub)
%     error('start value cannot be at the bounds')
% end
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


%% prevew data
tic
[I_cal_set,qset_true] = fcn_neutron_ref_hmcmc(x1,qset);     % preview initial conditions before fitting
toc

normFlag = 1;   % 1/2: q^4/none

% --- plot
figure
hold on; box on;
% plot data
for iset = 1:nsets  
    if normFlag == 1
        errorbar(qset{iset},Iset{iset}.*qset_true{iset}.^4*plotshift^(iset-1),Ierr_set{iset}.*qset_true{iset}.^4*plotshift^(iset-1),'o');
    elseif normFlag == 2
        errorbar(qset{iset},Iset{iset}*plotshift^(iset-1),Ierr_set{iset}*plotshift^(iset-1),'o');           
    end
end
% plot calculation
for iset = 1:nsets
    if normFlag == 1
        plot(qset{iset},I_cal_set{iset}.*qset_true{iset}.^4*plotshift^(iset-1),'k-','linewidth',1.5);
    elseif normFlag == 2
        plot(qset{iset},I_cal_set{iset}*plotshift^(iset-1),'k-','linewidth',1.5);
    end
end
set(gca,'xscale','linear');
set(gca,'yscale','log');
% legend(strcat('data set',cellstr(num2str((1:nsets)'))),'Location','best')
legend(flist,'Location','best','Interpreter','none')

if plotEDPFlag == 1
    return 
end


%% Define potential, kinetic energies and their differentials
% --- define potention  and differential
U = @(X)potentialE(X,nsets,qset,Iset,Iweight_set);
opts = admOptions('i', 1);  % only differentiate the 1st input argument (not param)
% dU = @(X)admDiffFor(@potentialE,1,X,nsets,qset,Iset,Iweight_set,opts);   % Forward mode
% dU = @(X)admDiffVFor(@potentialE,1,X,nsets,qset,Iset,Iweight_set,opts);     % Vectorized forward mode
% dU = @(X)admDiffRev(@potentialE,1,X,nsets,qset,Iset,Iweight_set,opts);    % Reverse mode
dU = @(X)admDiffFD(@potentialE,1,X,nsets,qset,Iset,Iweight_set,opts);     % Numerical finite difference method (not an auto diff method)

% --- define mass matrix
M = eye(length(x1)+nsets);  % default
if nutsFlag == 1 % when NUTS is choosen, force M to be unit
    M = eye(length(x1)+nsets); 
end

% --- define kinetic energy and differential
% K = @(p)p*inv(M)*p'/2 + log(det(M))/2; % Riemannian-Gaussian kinetic energy
K = @(p)p*inv(M)*p'/2 + log(det(M));    % Euclidean-Gaussian kinetic energy
dK = @(p)p*inv(M);          % differation of kinetic energy w.r.t momentum

%% Initialize parameters 
plotEDPFlag = 0;    % disable online plot of EDP

% ---  initialize variance of residual to guessed values
Sigmavi = nan(1,nsets);
for iset=1:nsets   
    Sigmavi(iset) = var(log10(I_cal_set{iset}) - log10(Iset{iset}));
end

% --- initial U and dU values
xp_init = [x1(:)',log(Sigmavi)];
[dPotE_init,PotE_init] = dU(xp_init);

%% --- find good initial value for epsilon (may skip this part)
pp_init = mvnrnd(zeros(size(xp_init)),M);
logp_init = -PotE_init - K(pp_init);
epsilon_init = 0.1;   % large start value is preferred than small start value
%tic
%[epsilon,counter] = heuristic4epsilon(xp_init,pp_init,dPotE_init,epsilon_init,U,K,dU,dK,logp_init);
%fprintf('Inialized epsilon=%f, after %d iterations of heuristic done in %f sec\n',epsilon,counter,toc);

% --- force epsilon if a different inital value is desired
epsilon = 0.0004;

%% --- dual averaging settings
if dualAvgFlag == 1
    delta = 0.8;   % optimal acceptance probability
    gamma = 0.05;
    t0 = 10;
    kappa = 0.75;
    Madap = Inf;    % Inf is preferred especially in presence of multiple modes in the probablility distribution
    H_init = 0;
    Hm = H_init;
    logepsilon_init = log(10*epsilon);
    logepsilon_mean = log(epsilon);
end

%% --- NUTS or simple HMCMC settings
if nutsFlag == 1  
    delta_max = 500;        % stopping rule (equation 3), large enough so not to interfere
else
    lambda_target = 10*epsilon;    % target total number simulation length
end


%% start HMCMC
% --- assign initial values
xp = xp_init;
dPotE = dPotE_init;
logp_xp = -PotE_init;

% --- start MCMC
nmcmc = 5000;   % # of iterations
param_mcmc = nan(nmcmc,length(xp));    
epsilon_mcmc = nan(nmcmc,1);

% --- the 0th iteration to store the intial values
param_mcmc0 = xp_init;
epsilon_mcmc0 = epsilon;

warning('off','MATLAB:Figure:FigureSavedToMATFile');
tic
for ii=1:nmcmc
   
    
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
            fprintf('%d of %d iterations (%.1f sec) with epsilon=%.5f: ',ii,nmcmc,toc,epsilon);
%            fprintf(['[',repmat('%f ',1,6),'\b]\n'],xp([1:5,end]));     % only display 6 parameters
            fprintf(['[',repmat('%f ',1,length(xp)),'\b]\n'],xp);    
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
    
    % --- collect result
    param_mcmc(ii,:) = xp;
    epsilon_mcmc(ii) = epsilon;
    
    % --- update plot
    if updatePlotFlag == 1
        [I_cal_set,qset_true] = fcn_neutron_ref_hmcmc(xp(1:end-nsets),qset);     
        if ii==1
            hfig = figure('position',[200,200,600,700]);
            % --- plot curve
            hax1 = subplot(11,1,1:4);
            hold on; box on;
            for iset = 1:nsets
                errorbar(qset{iset},Iset{iset}.*qset_true{iset}.^4*plotshift^(iset-1),Ierr_set{iset}.*qset_true{iset}.^4*plotshift^(iset-1),'o');
            end
            % plot calculation
           %  hline1 = [];
            for iset = 1:nsets
                hline1(iset) = plot(qset{iset},I_cal_set{iset}.*qset_true{iset}.^4*plotshift^(iset-1),'k-','linewidth',1.5);
            end
            legend(flist,'Location','best','Interpreter','none')
            set(gca,'yscale','log');
            set(gca,'XMinorTick','on','YMinorTick','on');            
            xlabel('q');
            ylabel('R\timesq^4');
            % --- plot SLD
            hax2 = subplot(11,1,6:9);            
            hold on; box on;
            hline2 = plot(edp_0(:,1),edp_0(:,2),'r-','linewidth',1);
            ylabel('SLD');      
            xlabel('height (A)');            
            set(gca,'XMinorTick','on','YMinorTick','on');
            % --- plot epsilon
            hax4 =subplot(11,1,11);
            hline4 = plot(epsilon_mcmc);
            set(gca,'xlim',[0,nmcmc]);
            xlabel('iteration');
            ylabel('\epsilon');
            set(gca,'yscale','log');
            htitle = sgtitle(sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc));
        else
            for iset = 1:nsets
                hline1(iset).YData = I_cal_set{iset}.*qset_true{iset}.^4*plotshift^(iset-1);
            end
            set(hline2,'Color','b','linewidth',0.5); % change previous color to blue and linewidth to thin            
            hline2 = plot(edp_0(:,1),edp_0(:,2),'r-','linewidth',1,'Parent',hax2);  % add new line
            set(hline4,'XData',1:nmcmc,'YData',epsilon_mcmc);             
            htitle.String = sprintf('%d of %d: %.1f sec elapsed',ii,nmcmc,toc);
        end
        pause(0.01);        
    end
    
    save tmp.mat        % iteratively save result
end
total_time =toc;

% --- add inital values
param_mcmc = [param_mcmc0; param_mcmc];
epsilon_mcmc = [epsilon_mcmc0; epsilon_mcmc];

if nutsFlag == 1 && dualAvgFlag == 1
    save fresult_NUTS_DA_v1.mat
elseif nutsFlag == 1 && dualAvgFlag == 0
    save fresult_NUTS_NoDA_v1.mat    
elseif nutsFlag == 0 && dualAvgFlag == 1    
    save fresult_Simple_DA_v1.mat        
elseif nutsFlag == 0 && dualAvgFlag == 0
    save fresult_Simple_NoDA_v1.mat            
end

