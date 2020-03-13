%% ========================== post HMCMC analysis ==========================
load fresult_hmcmc_NUTS_NoDA_v1_5000mcmc.mat

%% trace plot
% % --- clear NaN (in case of incomplete iteration loops);
nmcmc = nnz(all(~isnan(param_mcmc),2));
param_mcmc(nmcmc+1:end,:) = [];
epsilon_mcmc(nmcmc+1:end) = [];

% --- inverse transform to normal scale
x1_mcmc = inversetransformx(param_mcmc(:,1:end-1)',lbub(logical(fitFlag),:))';
Sigmavi_mcmc = exp(param_mcmc(:,end));
param_mcmc_inverted = [x1_mcmc,Sigmavi_mcmc];

% --- trace plot of parameters
idx_plot = 1:size(param_mcmc,2);


% idx_plot = [1:8,size(param_mcmc,2)];
% ylabel_str = {'I_0','z_{offset}','d_{ptba}','\sigma_{dtba}','d_{ps}','\sigma_{ps}','d_{au}','f_e','\sigma^2'};


nsubplots = length(idx_plot)+1;
figure
for kk=1:length(idx_plot)
    subplot(nsubplots,1,kk);
    % --- plot transformed parameters
    % plot(param_mcmc(:,idx_plot(kk)));
    % --- plot inverse transformed (in true scale)
    plot(param_mcmc_inverted(:,idx_plot(kk)));
%    ylabel(ylabel_str{kk});
end
% --- plot epsilon
subplot(nsubplots,1,kk+1);
plot(epsilon_mcmc)
set(gca,'yscale','log')
ylabel('\epsilon');
xlabel('iteration');
linkaxes(findall(gcf,'type','axes'),'x');

%% Plot variance trace
figure
plot(Sigmavi_mcmc*1e3);
xlabel('Iteration');
set(gca,'XMinorTick','on','YMinorTick','on');
%set(gca,'xtick',[200:200:1500]);
set(gca,'xlim',[1,5000]);
ylabel('\sigma^2 (\times 10^{-3})');
set(gca,'TickLength',[0.02,0.02]);

%% collect result
I_cal_tmp = cell(nmcmc,nsets);
q_true_tmp = cell(nmcmc,nsets);
edp_all = cell(nmcmc,1);
for ii=1:nmcmc
    disp(ii);
    [I_cal_tmp(ii,:),q_true_tmp(ii,:)] = fcn_neutron_ref_hmcmc(param_mcmc(ii,1:end-nsets),qset);   
    edp_all{ii} = edp_0;
end
% --- convert to mat
I_cal_all = cell(1,nsets);
q_true_all = cell(1,nsets);
for iset = 1:nsets
    I_cal_all{iset} = cell2mat(I_cal_tmp(:,iset)');
    q_true_all{iset} = cell2mat(q_true_tmp(:,iset)');
end

save msplot_hmcmc.mat

%% statistics of parameters
load msplot_hmcmc.mat

%%
% --- burn-in
nburnin= 4000; 
nstable = nmcmc-nburnin;
param_postburnin = param_mcmc(nburnin+1:end,:);
param_mcmc_inverted_postburnin = param_mcmc_inverted(nburnin+1:end,:);

% --- mean and variance of parameters
E_param = mean(param_mcmc_inverted_postburnin);
V_param = var(param_mcmc_inverted_postburnin);

% --- MCMC standard error by effective sample size
[ESS,Sigma] = multiESS(param_mcmc_inverted_postburnin);    % book assume no correlation between thetas
%[ESS,Sigma] = multiESS(param_postburnin);    % book assume no correlation between thetas
ESS = round(ESS);
% MCMCSE = sqrt(diag(cov(param_mcmc))'/ESS);    % useing covariance of mcmc
MCMCSE = sqrt(diag(Sigma)'/ESS);    % useing covariance from mutliESS

array2table([E_param; sqrt(V_param);MCMCSE],'RowNames',{'Mean','SE','MCMCSE'})
fprintf('Effective sample size is %d\n',ESS);


%% matrix plot (slow for high dimension)
figure
[S,AX,BigAx,H,HAx] = plotmatrix(param_mcmc_inverted_postburnin); % plot parameters in normal scale
% [S,AX,BigAx,H,HAx] = plotmatrix(param_postburnin); % plot parameters in normal scale
set(AX,'xminortick','on','yminortick','on','ticklength',[0.04,0.04]);

list_label = {'I_0','I_b','\Delta\lambda','\sigma_6','\Delta\alpha',...
    'h_1','h_2','h_3','h_4','h_5',...
    'SLD_1','SLD_2','SLD_3','SLD_4','SLD_5',...
    '\sigma_1','\sigma_2','\sigma_3','\sigma_4','\sigma_5',...
    '\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5',...
    '\sigma^2'};
for ii=1:size(param_mcmc_inverted_postburnin,2)
%     ylabel(AX(ii,1),list_label{ii},'Rotation',90);
    ylabel(AX(ii,1),{list_label{ii};' '});
    xlabel(AX(end,ii),{' ';' ';list_label{ii}});    
end
set(H,'EdgeColor','none')


%% statistics analysis of parameter
x=zeros(length(fitFlag),nmcmc);
x(fitFlag==1,:) = param_mcmc(:,1:end-1)';
x(fitFlag==0,:) = repmat(x2,1,nmcmc);
x = inversetransformx(x,lbub);  % inverse transform is done for all parameters (including non-fitted)
x1_mcmc = x(fitFlag==1,:)';

x = x(:,nburnin+1:end);

Inorm           = x(1,:);
Ibk             = x(2,:);
resdawave       = x(3,:);
resconstant     = x(4,:);
alphai_offset   = x(5,:);
sld_env         = x(6,:);
sld_sub         = x(7,:);
sigma_sub       = x(8,:);
layeredp        = reshape(x(9:end,:),nlayer,[],nstable); % 3D matrix


param_other=[Inorm',Ibk',resdawave',alphai_offset',param_mcmc_inverted_postburnin(:,end)];
E_param = mean(param_other);
V_param = var(param_other);
array2table([E_param; sqrt(V_param); E_param./sqrt(V_param); (1-tcdf(E_param./sqrt(V_param),nmcmc))*2],'RowNames',{'Mean','SE','t-stat','p-value'},'VariableNames',{'I0','Ibk','res','alphai_offset','Sigmav'})



layeredp_a = squeeze(layeredp(2:end,1,:))';     % for height
layeredp_b = squeeze(layeredp(2:end,2,:))';     % for SLD
layeredp_c = [squeeze(layeredp(2:end,3,:))',sigma_sub(:)];     % for roughness
layeredp_d = squeeze(layeredp(2:end,4,:))';     % for skew
dof = length(xp)-1;     % degreee of freedom

stats = [];

E_param = mean(layeredp_a);
V_param = var(layeredp_a);
array2table([E_param; sqrt(V_param); E_param./sqrt(V_param); (1-tcdf(E_param./sqrt(V_param),nmcmc))*2],'RowNames',{'Mean','SE','t-stat','p-value'},'VariableNames',{'h1','h2','h3','h4','h5'})

E_param = mean(layeredp_b);
V_param = var(layeredp_b);
array2table([E_param; sqrt(V_param); E_param./sqrt(V_param); (1-tcdf(E_param./sqrt(V_param),nmcmc))*2],'RowNames',{'Mean','SE','t-stat','p-value'},'VariableNames',{'SLD1','SLD2','SLD3','SLD4','SLD5'})


E_param = mean(layeredp_c);
V_param = var(layeredp_c);
array2table([E_param; sqrt(V_param); E_param./sqrt(V_param); (1-tcdf(E_param./sqrt(V_param),nmcmc))*2],'RowNames',{'Mean','SE','t-stat','p-value'},'VariableNames',{'sigma1','sigma2','sigma3','sigma4','sigma5','sigma6'})


E_param = mean(layeredp_d);
V_param = var(layeredp_d);
array2table([E_param; sqrt(V_param); E_param./sqrt(V_param); (1-tcdf(E_param./sqrt(V_param),nmcmc))*2],'RowNames',{'Mean','SE','t-stat','p-value'},'VariableNames',{'alpha1','alpha2','alpha3','alpha4','alpha5'})


%% --- plot reflectivity
% --- get post-burin
edp_post = edp_all(nburnin+1:end);
I_cal_post =  cell(size(I_cal_all));
q_true_post = cell(size(q_true_all));
for iset=1:nsets
    I_cal_post{iset} = I_cal_all{iset}(:,nburnin+1:end);    
    q_true_post{iset} = q_true_all{iset}(:,nburnin+1:end);
end

normFlag = 1;   % 1/2: q^4/log

figure('pos',[573,10,500,700]);
subplot(3,1,1);
plot(Sigmavi_mcmc*1e3);
xlabel('Iteration');
set(gca,'XMinorTick','on','YMinorTick','on');
%set(gca,'xtick',[200:200:1500]);
set(gca,'xlim',[1,5000]);
ylabel('\sigma^2 (\times 10^{-3})');
set(gca,'TickLength',[0.02,0.02]);

subplot(3,1,2);
hold on; box on;
for iset = 1:nsets
    qz = qset{iset};
    qdata = q_true_post{iset};   
    
    Idata = Iset{iset};
    Ierr = Ierr_set{iset};
    I_cal = I_cal_post{iset};
    
    if normFlag == 1 % scale with qz^4
        errorbar(qz,Idata.*qz.^4*plotshift^(iset-1) * 1e7,Ierr.*qz.^4*plotshift^(iset-1)* 1e7,'o','markersize',5,'CapSize',2,'LineWidth',0.5);  % data
        I_cal_by_q4 = I_cal.*qdata.^4;
        lb_I_cal = quantile(I_cal_by_q4,0.025,2);   % 2.5% quantile
        ub_I_cal = quantile(I_cal_by_q4,0.975,2);   % 97.5% quantile
        patch_x = [qz; flip(qz)];
        patch_y = [lb_I_cal;flip(ub_I_cal)];
        patch(patch_x,patch_y*plotshift^(iset-1)*1e7,'m','EdgeColor','m','FaceAlpha',.5);
    elseif normFlag == 2  % --- no scale
        errorbar(qz,Idata*plotshift^(iset-1),Ierr*plotshift^(iset-1),'o','markersize',4,'CapSize',2,'LineWidth',0.5);  % data
        lb_I_cal = quantile(I_cal,0.025,2);   % 2.5% quantile
        ub_I_cal = quantile(I_cal,0.975,2);   % 97.5% quantile
        patch_x = [qz; flip(qz)];
        patch_y = [lb_I_cal;flip(ub_I_cal)];
        patch(patch_x,patch_y*plotshift^(iset-1),'m','EdgeColor','m','FaceAlpha',.5);        
    end    
end
if normFlag == 1    
    set(gca,'yscale','log');
    ylabel(['R\timesq_z^4 (a.u.\times',char(197),'^{-4})']);
elseif normFlag == 2
    set(gca,'yscale','log');
    ylabel('I');        
end
set(gca,'xminortick','on','YMinorTick','on');
xlabel(['q_z (',char(197),'^{-1})']);
set(gca,'TickLength',[0.025,0.025]);
set(gca,'xlim',[0,0.085]);
legend({'data','95% confidence'},'box','off');
%hlines = get(gca,'Children');   
%legend(hlines(end:-2:1),flist,'location','best','interpreter','none')


% --- plot SLD
z_mean = linspace(min(cellfun(@(x)x(1,1)',edp_post)),max(cellfun(@(x)x(end,1)',edp_post)),size(edp_post{1},1));
z_mean = z_mean(:);

subplot(3,1,3);
hold on; box on;
sld_interp = cell2mat(cellfun(@(x)interp1(x(:,1),x(:,2),z_mean)',edp_post,'UniformOutput',false))';
lb_sld = quantile(sld_interp,0.025,2);   % 2.5% quantile
ub_sld = quantile(sld_interp,0.975,2);   % 97.5% quantile
sld_median = nanmedian(sld_interp,2);
sld_mean = nanmean(sld_interp,2);
sld_mode = mode(sld_interp,2);
patch_x = [z_mean; flip(z_mean)];
patch_y = [lb_sld;flip(ub_sld)];
patch(patch_x,patch_y*1e6,'m','EdgeColor','none','FaceAlpha',.5);
plot(z_mean,sld_median*1e6,'-','LineWidth',1);
% plot(z_mean,sld_mean,'-','LineWidth',1.5);
% plot(z_mean,sld_mode,'-','LineWidth',1.5);
xlabel(['Height above Si (',char(197),')']);
ylabel(['SLD (\times10^{-6} ',char(197),'^{-2})']);
set(gca,'TickLength',[0.02,0.02]);
set(gca,'xlim',[-200,1500]);
set(gca,'xminortick','on','YMinorTick','on');
legend({'95% confidence','median'},'box','off','location','best');