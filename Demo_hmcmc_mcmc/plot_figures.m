load('msplot_MCMC_singlemode_boxplot.mat');

figure('pos',[570,430,400,420]);
m = 5;
n = 2;
subplot(m,n,[1,3]);
boxplot(ESS_all,epsilon_list(1,:),'Notch','on','OutlierSize',3);
set(gca,'XMinorTick','off','YMinorTick','on');
set(gca,'TickLength',[0.04,0.04]);
set(gca,'ylim',[20,120]);
xlabel('\Delta');
ylabel('ESS');
xtickangle(90)

load('msplot_HMCMC_singlemode_boxplot.mat');
subplot(m,n,[2,4]);
boxplot(ESS_all,epsilon_list(1,:),'Notch','on','OutlierSize',3);
set(gca,'XMinorTick','off','YMinorTick','on');
set(gca,'TickLength',[0.04,0.04]);
set(gca,'ylim',[50,500]);
xlabel('\epsilon');
ylabel('ESS');
xtickangle(90);

load('msplot_MCMC_single.mat');
subplot(m,n,[5,7]);
hold on; box on;
[~,hc] = contour(X,Y,Z,5,'b','LineWidth',1);
plot(xp_init(1),xp_init(2),'o');
plot([xp_init(1);param_mcmc(:,1)],[xp_init(2);param_mcmc(:,2)],'r-')
xlabel('x_1');
ylabel('x_2');
set(gca,'TickLength',[0.04,0.04]);
set(gca,'XMinorTick','on','YMinorTick','on');
uistack(hc,'top');  


subplot(m,n,9);
plot(param_mcmc);
xlabel('iteration');
% legend({'x_1','x_2'},'Location','best');
set(gca,'xlim',[1,nmcmc])
set(gca,'TickLength',[0.04,0.04]);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'ylim',[0,15]);

load('msplot_HMCMC_single.mat');
subplot(m,n,[6,8]);
hold on; box on;
[~,hc] = contour(X,Y,Z,5,'b','LineWidth',1);
plot(xp_init(1),xp_init(2),'o');
plot([xp_init(1);param_mcmc(:,1)],[xp_init(2);param_mcmc(:,2)],'r-')
xlabel('x_1');
ylabel('x_2');
set(gca,'TickLength',[0.04,0.04]);
set(gca,'XMinorTick','on','YMinorTick','on');
uistack(hc,'top');  

subplot(m,n,10);
plot(param_mcmc);
xlabel('iteration');
set(gca,'xlim',[1,nmcmc])
set(gca,'TickLength',[0.04,0.04]);
set(gca,'XMinorTick','on','YMinorTick','on');
set(gca,'ylim',[0,15]);
