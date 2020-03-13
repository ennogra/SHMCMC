%% ========================== post MCMC analysis ==========================
load fresult_shmcmc_NUTS_NoDA_v1_5000mcmc.mat

%% add initial values
dPotE_mcmc = [dPotE_init; dPotE_mcmc];

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
figure('position',[100,100,450,300]);
ha = tight_subplot(size(param_mcmc,2),1,[0.01 0],[.1 .05],[.125 .025]);
axes(ha(1));
plot([0:nmcmc-1],param_mcmc_inverted(:,1));
set(gca,'ylim',[1.1,1.5]);
%ylabel('I_0');
legend('I_0','Location','best','box','off');

axes(ha(2));
plot([0:nmcmc-1],param_mcmc_inverted(:,2)*1e5);
set(gca,'ylim',[0,3]);
%ylabel({'I_b','(10^{-5})'});
legend('I_b (10^{-5})','Location','best','box','off');

axes(ha(3));
plot(0:nmcmc-1,param_mcmc_inverted(:,3));
set(gca,'ylim',[730,790]);
%ylabel(['R (',char(197),')']);
legend(['R (',char(197),')'],'Location','best','box','off');

axes(ha(4));
plot([0:nmcmc-1],param_mcmc_inverted(:,4));
set(gca,'ylim',[46,78]);
%ylabel(['\sigma_R (',char(197),')']);
legend(['\sigma_R (',char(197),')'],'Location','best','box','off');

axes(ha(5));
plot([0:nmcmc-1],param_mcmc_inverted(:,5)*1e5);
set(gca,'ylim',[0,30]);
%ylabel({'\sigma_q', ['(10^{-5} ',char(197),'^{-1})']});
legend(['\sigma_q (10^{-5} ',char(197),'^{-1})'],'Location','best','box','off');

axes(ha(6));
plot([0:nmcmc-1],param_mcmc_inverted(:,6)*1e3);
set(gca,'ylim',[0,6]);
%ylabel({'\sigma^2', '(10^{-3})'});
legend('\sigma^2 (10^{-3})','Location','best','box','off');
 
xlabel('iteration');
set(ha,'XMinorTick','on');
set(ha,'YMinorTick','on'); 
linkaxes(ha,'x');
set(ha,'XTick',[0:200:800]); 
set(ha,'xlim',[0,1000]);
set(ha,'ticklength',[0.015,0.015]);
set(ha(1:end-1),'XTickLabel','');


%% --- plot dPotE
figure('pos',[500,100,450,200]);
imagesc(abs(dPotE_mcmc)',[0,100]);
colormap('jet');
xlabel('Iteration');
%ylabel('Index of parameter');
set(gca,'TickLength',[0.02,0.02]);
set(gca,'XMinorTick','on','YMinorTick','off');
set(gca,'TickDir','out');
%set(gca,'xtick',[100:100:900])
set(gca,'ytick',1:6);
set(gca,'yticklabel',{'I_0','I_b','R','\sigma_R','\sigma_q','\sigma^2'});
text(3800,4.8,'|\nablaU(x)|','Interpreter','tex','Color','w','FontSize',9);
hbar = colorbar('south');
set(hbar,'Color','w','Position',[0.65,0.31,0.22,0.05]);

%% --- burn-in
nburnin=51;
nstable = nmcmc-nburnin;
param_postburnin = param_mcmc(nburnin+1:end,:);
param_mcmc_inverted_postburnin = param_mcmc_inverted(nburnin+1:end,:);

% --- mean and variance of parameters
E_param = mean(param_mcmc_inverted_postburnin);
V_param = var(param_mcmc_inverted_postburnin);

array2table([E_param; sqrt(V_param); E_param./sqrt(V_param); (1-tcdf(E_param./sqrt(V_param),nmcmc))*2],'RowNames',{'Mean','SE','t-stat','p-value'},'VariableNames',{'I0','bk','R','sigma_R','sigma_q','variance'})


% --- MCMC standard error by effective sample size
[ESS,Sigma] = multiESS(param_mcmc_inverted_postburnin);    % book assume no correlation between thetas
% [ESS,Sigma] = multiESS(param_postburnin);    % book assume no correlation between thetas

ESS = round(ESS);
% MCMCSE = sqrt(diag(cov(param_mcmc))'/ESS);    % useing covariance of mcmc
MCMCSE = sqrt(diag(Sigma)'/ESS);    % useing covariance from mutliESS

array2table([E_param; sqrt(V_param);MCMCSE],'VariableNames',{'I0','bk','R','sigma_R','sigma_q','variance'},'RowNames',{'Mean','SE','MCMCSE'})
fprintf('Effective sample size is %d\n',ESS);

% --- matrix plot
param_mcmc_inverted_postburnin_scaled = param_mcmc_inverted_postburnin;
param_mcmc_inverted_postburnin_scaled(:,2) = param_mcmc_inverted_postburnin_scaled(:,2)*1e5;
param_mcmc_inverted_postburnin_scaled(:,5) = param_mcmc_inverted_postburnin_scaled(:,5)*1e5;
param_mcmc_inverted_postburnin_scaled(:,6) = param_mcmc_inverted_postburnin_scaled(:,6)*1e3;

figure('position',[200,125,450,400]);
%[S,AX,BigAx,H,HAx] = plotmatrix(param_postburnin);    % plot transformed parameters
[S,AX,BigAx,H,HAx] = plotmatrix(param_mcmc_inverted_postburnin_scaled); % plot parameters in normal scale
set(AX,'xminortick','on','yminortick','on','ticklength',[0.05,0.05]);
ylabel(AX(1,1),'I_0');
ylabel(AX(2,1),'I_b (10^{-5})');
ylabel(AX(3,1),['R (',char(197),')']);
ylabel(AX(4,1),['\sigma_R (',char(197),')']);
ylabel(AX(5,1),{'\sigma_q',['(10^{-5}',char(197),'^{-1})']});
ylabel(AX(6,1),'\sigma^2 (10^{-3})');
xlabel(AX(end,1),'I_0');
xlabel(AX(end,2),'I_b (10^{-5})');
xlabel(AX(end,3),['R (',char(197),')']);
xlabel(AX(end,4),['\sigma_R (',char(197),')']);
xlabel(AX(end,5),['\sigma_q (10^{-5}',char(197),'^{-1})']);
xlabel(AX(end,6),'\sigma^2 (10^{-3})');
set(H,'EdgeColor','none')
set(AX,'ticklength',[0.05,0.05]);


%% load result
I_cal_all = nan(length(qdata),nmcmc);
for ii=1:nmcmc
    disp(ii) 
    I_cal_all(:,ii) = fcn_saxs_sphere_hmcmc(param_mcmc(ii,1:end-1),qdata);     % preview initial conditions before fitting
end
save msplot_shmcmc.mat

%% plot saxs
load msplot_shmcmc.mat

% plot selected saxs trace
idx_plot = [1,4,11,21];
figure('position',[100,225,400,350]);
hold on; box on;
idxdata=[1:1:70,71:2:length(qdata)];
errorbar(qdata(idxdata),Idata(idxdata),Idata_err(idxdata),'o','MarkerSize',6);
% plot(-z_mean,delta_median,'-','LineWidth',1.5);
% plot(-z_mean,delta_mean,'-','LineWidth',1.5);
%newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
colororder('default')
plot(qdata,I_cal_all(:,idx_plot),'-','LineWidth',1.5);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'XMinorTick','on','YMinorTick','on');
xlabel(['q (',char(197),'^{-1})']);
ylabel('I (a.u.)');
set(gca,'ylim',[1e-5,2]);
set(gca,'ticklength',[0.025,0.025]);
%set(gca,'YAxisLocation','right','XAxisLocation','top');
legend(['data';'initial'; cellstr(strcat('#',num2str(idx_plot(2:end)'-1)))],'Location','best','Box','off');


%% plot stationary
I_cal = I_cal_all(:,nburnin+1:end); % (stable ones)
lb_I_cal = quantile(I_cal,0.025,2);   % 2.5% quantile
ub_I_cal = quantile(I_cal,0.975,2);   % 97.5% quantile
I_cal_median = nanmedian(I_cal,2);
I_cal_mean = nanmean(I_cal,2);
I_cal_mode = mode(I_cal,2);
patch_x = [qdata; flip(qdata)];
patch_y = [lb_I_cal;flip(ub_I_cal)];

% --- plot all stabled
figure('position',[573,125,450,350]);
hold on; box on;
errorbar(qdata,Idata,Idata_err,'o','MarkerSize',6);
patch(patch_x,patch_y,'m','EdgeColor','none','FaceAlpha',1);
% plot(-z_mean,delta_median,'-','LineWidth',1.5);
% plot(-z_mean,delta_mean,'-','LineWidth',1.5);
%plot(qdata,I_cal_mode,'r-','LineWidth',1.5);
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'XMinorTick','on','YMinorTick','on');
xlabel(['q (',char(197),'^{-1})']);
ylabel('I (a.u.)');
set(gca,'ticklength',[0.025,0.025]);
set(gca,'YAxisLocation','right','XAxisLocation','top');


%% plot saxs trace and stationary in one figure
% --- plot all stabled
figure('position',[573,125,1000,200]);
%ha = tight_subplot(2,1,[0 0],[.1 .01],[.125 .025]);
idxdata=[1:1:70,71:2:length(qdata)];

ha = [];
ha(1)=subplot(1,2,1);
hold on; box on;
errorbar(qdata(idxdata),Idata(idxdata),Idata_err(idxdata),'o','MarkerSize',6);
% plot(-z_mean,delta_median,'-','LineWidth',1.5);
% plot(-z_mean,delta_mean,'-','LineWidth',1.5);
%newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
colororder('default')
plot(qdata,I_cal_all(:,idx_plot),'-','LineWidth',1.5);
set(gca,'xscale','log');
set(gca,'yscale','log');
ylabel('I (a.u.)');
%set(gca,'ylim',[1e-5,2]);
xlabel(['q (',char(197),'^{-1})']);
%set(gca,'YAxisLocation','right','XAxisLocation','top');
legend(['data';'initial'; cellstr(strcat('#',num2str(idx_plot(2:end)'-1)))],'Location','best','Box','off');
set(gca,'ticklength',[0.025,0.025]);
set(gca,'XMinorTick','on','YMinorTick','on');


ha(1)=subplot(1,2,2);
hold on; box on;
errorbar(qdata(idxdata),Idata(idxdata),Idata_err(idxdata),'o','MarkerSize',6);
patch(patch_x,patch_y,'m','EdgeColor','m','FaceAlpha',1);
% plot(-z_mean,delta_median,'-','LineWidth',1.5);
% plot(-z_mean,delta_mean,'-','LineWidth',1.5);
%plot(qdata,I_cal_mode,'r-','LineWidth',1.5);
set(gca,'xscale','log');
set(gca,'yscale','log');
%set(gca,'YAxisLocation','right','XAxisLocation','top');
xlabel(['q (',char(197),'^{-1})']);
set(gca,'yticklabel','');
%set(gca,'ylim',[1e-5,2]);
set(gca,'XMinorTick','on','YMinorTick','on');

%linkaxes(ha,'xy');
set(gca,'ticklength',[0.025,0.025]);
%set(ha,'ytick',[1e-4,1e-2,1])
legend({'data','95% confidence'},'Location','southwest','box','off');

