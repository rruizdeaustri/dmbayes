clear;
gclim=[0 1e31];
redf = 0.75;
bflw = 2.5;
bfsize = 25;
lab_fontsize = 24/redf; axes_fontsize = 24/redf; names_fontsize = 10/redf;  cbar_fontsize = 16/redf;  textbox_fontsize = 30/redf; twolines_fontsize = 30/redf;exp_fontsize = 14/redf;  
 l1='-k'; lw1=0.8; lwaxes = 2.0;
 lwdd= 3.0; 
 l2='.r'; lw1=1.2; 
figure
subplot(    4,    3,    1);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.17188E+03 0.0];
mean = [    0.14763E+03, 0.0];
refpoint = [    0.33800E+02 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p1.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','m_\chi (GeV)','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_1.ps;
figure
subplot(    4,    3,    2);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.24495E+01 0.0];
mean = [    0.24589E+01, 0.0];
refpoint = [    0.27520E+01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p2.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','log[<\sigma v> (10^{27} cm^3 s^{-1})]','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_2.ps;
figure
subplot(    4,    3,    3);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.43897E+00 0.0];
mean = [    0.20122E+00, 0.0];
refpoint = [    0.92400E+00 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p3.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{b \bar{b}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_3.ps;
figure
subplot(    4,    3,    4);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.12562E-01 0.0];
mean = [    0.92674E-01, 0.0];
refpoint = [   -0.40000E+01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p5.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{W^+ W^-}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_4.ps;
figure
subplot(    4,    3,    5);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.44399E-02 0.0];
mean = [    0.47049E-01, 0.0];
refpoint = [   -0.40000E+01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p6.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{Z Z}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_5.ps;
figure
subplot(    4,    3,    6);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.27197E-02 0.0];
mean = [    0.41745E-01, 0.0];
refpoint = [   -0.40000E+01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p7.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{e^+ e^-}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_6.ps;
figure
subplot(    4,    3,    7);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.25679E-01 0.0];
mean = [    0.43472E-01, 0.0];
refpoint = [   -0.40000E+01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p8.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{\mu^+ \mu^-}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_7.ps;
figure
subplot(    4,    3,    8);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.32960E+00 0.0];
mean = [    0.27232E+00, 0.0];
refpoint = [   -0.40000E+01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p9.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{c \bar{c}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_8.ps;
figure
subplot(    4,    3,    9);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.37298E-02 0.0];
mean = [    0.28948E-01, 0.0];
refpoint = [   -0.40000E+01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p10.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{g g}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_9.ps;
figure
subplot(    4,    3,   10);
set(gcf, 'NextPlot', 'replacechildren');
box on; axis square; hold on;
bestfit = [    0.18230E+00 0.0];
mean = [    0.27258E+00, 0.0];
refpoint = [    0.76000E-01 0.0];
tmp = load (fullfile('../../plot_data/sim_rec_logBR_thr3/','sim_rec_logBR_thr3_p20.dat'));
x = tmp(:,1);
pdf = tmp(:,2);
hbar = bar(x,pdf);
set(hbar, 'EdgeColor', 'b','FaceColor', [0.4 0.8 1.0],'BarWidth', 1.0);
colormap(cool); hold on;
 text('String','BR_{\tau^+\tau^-}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('Probability', 'FontSize',lab_fontsize);
axis tight;
set(gca,'YLim', [0.0 1.1]);
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
set(gcf,'paperpositionmode','auto');
set(gcf, 'Units','centimeters');
set(gcf, 'Position',[ 0.0 0.0 16.5 16.5]);
set(gca, 'Units','centimeters');
set(gca,'position',[4.0 4.0 11.0 11.0]);
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperSize',[16.5 16.5]);
hold off;
print -dpsc2 -loose  ../../output_files/sim_rec_logBR_thr3/sim_rec_logBR_thr3_1D_10.ps;