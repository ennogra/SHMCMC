% --- needed for the computation of refraction index
addpath(fullfile(pwd,'xrayrefraction'),fullfile(pwd,'xrayrefraction','AtomicScatteringFactor'));


%% ========================== post MCMC analysis ==========================
load fresult_hmcmc_Simple_NoDA_v1_2000mcmc.mat

%% add initial values
param_mcmc = [param_mcmc0; param_mcmc];
epsilon_mcmc = [epslion_mcmc0; epsilon_mcmc];
dPotE_mcmc = [dPotE_init; dPotE_mcmc];

%% trace plot
% % --- clear NaN (in case of incomplete iteration loops);
nmcmc = nnz(all(~isnan(param_mcmc),2));
param_mcmc(nmcmc+1:end,:) = [];
epsilon_mcmc(nmcmc+1:end) = [];
dPotE_mcmc(nmcmc+1:end,:) = [];

% --- inverse transform to normal scale
x1_mcmc = inversetransformx(param_mcmc(:,1:end-1)',lbub(logical(fitFlag),:))';
Sigmavi_mcmc = exp(param_mcmc(:,end));
param_mcmc_inverted = [x1_mcmc,Sigmavi_mcmc];

% --- trace plot of parameters
%idx_plot = 1:size(param_mcmc,2);


idx_plot = [1:8,size(param_mcmc,2)];
ylabel_str = {'I_0','z_{offset}','d_{ptba}','\sigma_{dtba}','d_{ps}','\sigma_{ps}','d_{au}','f_e','\sigma^2'};


nsubplots = length(idx_plot)+1;
figure
for kk=1:length(idx_plot)
    subplot(nsubplots,1,kk);
    % --- plot transformed parameters
    % plot(param_mcmc(:,idx_plot(kk)));
    % --- plot inverse transformed (in true scale)
    plot(param_mcmc_inverted(:,idx_plot(kk)));
    ylabel(ylabel_str{kk});
end
% --- plot epsilon
subplot(nsubplots,1,kk+1);
plot(epsilon_mcmc)
set(gca,'yscale','log')
ylabel('\epsilon');
xlabel('iteration');
linkaxes(findall(gcf,'type','axes'),'x');

%% --- plot dPotE
figure('pos',[500,100,430,230]);
imagesc(abs(dPotE_mcmc)',[0,200]);
colormap('jet');
xlabel('Iteration');
ylabel('Index of parameter');
set(gca,'TickLength',[0.02,0.02]);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'TickDir','out');
set(gca,'xtick',[100:200:1900])
text(500,22,'|\nablaU(x)|','Interpreter','tex','Color','w','FontSize',9);
hbar = colorbar('south');
set(hbar,'Color','w','Position',[0.65,0.4,0.22,0.05]);

%% collect x-ray results
Is_cal_all = nan(length(Is),nmcmc);
Y1_all = nan(length(Is),nmcmc);
Y2_all = nan(length(Is),nmcmc);
Y_elastic_all = nan(length(Is),nmcmc);
rho_all = cell(1,nmcmc);
nSlice = 800;  
for ii=1:nmcmc
    fprintf('Calculating %d of %d\n',ii,nmcmc);
    [Is_cal_all(:,ii),Y1_all(:,ii),Y2_all(:,ii),Y_elastic_all(:,ii)] = fcn_flu_spline_hmcmc(param_mcmc(ii,1:end-1),vzs); 
    rho_all{ii} = rho;
end
save msplot_hmcmc.mat

%% load resutl back and plot
load msplot_hmcmc.mat
 
vz0_mean = vz0;
Inorm_mean = Inorm;
alphaf_list = atand((vzs-vz0_mean)/sdd);

zlist = linspace(min(cellfun(@(x)x(end,1),rho_all)),max(cellfun(@(x)x(1,1),rho_all)),size(rho,1));
zlist = zlist(:);
edp_interp = cell2mat(cellfun(@(x)interp1(x(:,1),x(:,2),zlist),rho_all,'UniformOutput',false));
phi_interp = cell2mat(cellfun(@(x)interp1(x(:,1),x(:,3),zlist),rho_all,'UniformOutput',false));


%idx_plot = [1,6,11,21,41,81,161];%,101,201];   % index to be plotted
idx_plot = 1:nmcmc;   % index to be plotted
nidx = length(idx_plot);

yplot_scale = (Is(end)-Ibk)/Inorm_mean;


figure('pos',[500,200,400,450]);
subplot(3,1,1)  % plot sigma2 trace
plot(Sigmavi_mcmc);
set(gca,'xminortick','on','Yminortick','on');
set(gca,'yscale','log');
ylabel('log(\sigma^2)');

set(gca,'TickLength',[0.025,0.025]);
xlabel('Iteration');
ylabel('\sigma^2');

set(gca,'xlim',[0,2000]);

subplot(3,1,2);
hold on; box on;
color_list =  jet(nidx);
color_list(:,4) = 0.05; 
plot(alphaf_list,(Is-Ibk)/Inorm_mean/yplot_scale,'o','markersize',4,'color',[0 0.4470 0.7410,1]);
for ii=1:nidx
    plot(alphaf_list,(Is_cal_all(:,idx_plot(ii))-Ibk)/Inorm_mean/yplot_scale,'linewidth',.5,'color',color_list(ii,:));
end

set(gca,'XTickLabelMode','auto')
set(gca,'YTickLabelMode','auto')
set(gca,'TickLength',[0.025,0.025]);
set(gca,'xlim',[0,0.8]);
set(gca,'ylim',[0,7]);
set(gca,'xminortick','on')
set(gca,'yminortick','on')
ylabel('XWFH (a.u.)');
xlabel('Exit Angle (deg)');

colormap('jet');
hcbar = colorbar('north');
hcbar.Ticks=[0,0.5,1];
hcbar.TickLabels = [0,round((nmcmc-1)/2),nmcmc-1];


subplot(3,1,3);
hold on; box on;
for ii=1:nidx
    plot(zlist,phi_interp(:,idx_plot(ii)),'linewidth',.5,'Color',color_list(ii,:));
end
set(gca,'xminortick','on')
set(gca,'yminortick','on')
set(gca,'TickLength',[0.025,0.025]);
set(gca,'xlim',[0,800]);
%set(gca,'ylim',[0,0.025]);
xlabel(['Height above Pd mirror (',char(197),')']);
% set(gca,'XTick',0:200:800);
% set(gca,'XTickLabel',0:200:800);
ylabel(['\phi_{au} (% ',char(197),'^{-1})']);

%% making a movie
v = VideoWriter('HMMC_XWFH_100fps_2000iterations.avi');
v.Quality = 10;
v.FrameRate = 100;
open(v);
hfig = figure('pos',[500,200,400,300]);
for ii=1:nidx
    I_tmp = (Is_cal_all(:,idx_plot(ii))-Ibk)/Inorm_mean/yplot_scale;    
    if ii==1
        ha1 = subplot(2,1,1);
        hold on; box on;
        plot(alphaf_list,(Is-Ibk)/Inorm_mean/yplot_scale,'o','markersize',4,'color',[0 0.4470 0.7410,1]);
        hline1 = plot(alphaf_list,I_tmp,'linewidth',1.5,'Parent',ha1);        
        set(gca,'XTickLabelMode','auto')
        set(gca,'YTickLabelMode','auto')
        set(gca,'TickLength',[0.025,0.025]);
        set(gca,'xlim',[0,0.8]);
        set(gca,'ylim',[0,7]);
        set(gca,'xminortick','on')
        set(gca,'yminortick','on')
        ylabel('XWFH (a.u.)');
        xlabel('Exit Angle (deg)');
        
        ha2 = subplot(2,1,2);
        box on;
        hline2 = plot(zlist,phi_interp(:,idx_plot(ii)),'linewidth',1.5,'Parent',ha2);
        set(gca,'xminortick','on')
        set(gca,'yminortick','on')
        set(gca,'TickLength',[0.025,0.025]);
        set(gca,'xlim',[0,800]);
        set(gca,'ylim',[0,0.025]);
        xlabel(['Height above Pd mirror (',char(197),')']);
        % set(gca,'XTick',0:200:800);
        % set(gca,'XTickLabel',0:200:800);
        ylabel(['\phi_{au} (% ',char(197),'^{-1})']);
    else
        hline1.YData = I_tmp;
        hline2.YData = phi_interp(:,idx_plot(ii));
    end
    title(['Iteration # ',num2str(ii-1)],'Parent',ha1);            
    frame = getframe(hfig);
    writeVideo(v,frame);
end
close(v);