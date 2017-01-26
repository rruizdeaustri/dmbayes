clear;
gclim=[0 1e31];
figure(   1);  clf('reset')
figure(   2);  clf('reset')
figure(   5);  clf('reset')
redf = (1+0.35*   2);
bflw = 2.5;
bfsize = 15;
lab_fontsize = 24/redf; axes_fontsize = 24/redf; names_fontsize = 10/redf;  cbar_fontsize = 16/redf;  textbox_fontsize = 30/redf; twolines_fontsize = 30/redf;exp_fontsize = 14/redf;  
 l1='-k'; lw1=0.8; lwaxes = 2.0;
 lwdd= 3.0; 
 l2='.r'; lw1=1.2; 
load('../../colormaps/greenmap','likegree');
load('../../colormaps/ModInverseHot','modrhot');
load('../../colormaps/yellowblue','yellowblue');
load('../../colormaps/yellowred','yellowred');
figure(   1)
subplot(    2,    2,    1);
box on; axis square; hold on;
bestfit = [    0.21243E+01    0.60887E+02];
mean = [    0.21835E+01    0.67093E+02];
refpoint = [    0.27520E+01    0.33800E+02];
tmp = load (fullfile('../../plot_data/test/','test_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_1_2_cont'));
pdfpts=load (fullfile('../../plot_data/test/','test_2D_1_2_marg'));
contourf(x1,x2,pdfpts,64);
[C h] = contour(x1,x2,pdfpts,cnt(           1,:),l1);set(h,'LineWidth',lw1);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','log[<\sigma v> (10^{27} cm^3 s^{-1})]','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('m_\chi (GeV)','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = modrhot;
description = 'Posterior pdf';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.25], 'HorizontalAlignment', 'right');
description = 'Log BR  priors';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.15], 'HorizontalAlignment', 'right');
text('String', 'CMSSM, \mu>0' ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.05], 'HorizontalAlignment', 'right');
colormap(chosencolor);
figure(   4)
subplot(    2,    2,    1);
box on; axis square; hold on;
bestfit = [    0.21243E+01    0.60887E+02];
mean = [    0.21835E+01    0.67093E+02];
refpoint = [    0.27520E+01    0.33800E+02];
tmp = load (fullfile('../../plot_data/test/','test_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_1_2_cont'));
profpts = load (fullfile('../../plot_data/test/','test_2D_1_2_profl'));
ah=imagesc(x1,x2,profpts);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','log[<\sigma v> (10^{27} cm^3 s^{-1})]','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('m_\chi (GeV)','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = likegree;
description = 'Profile likelihood';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.25], 'HorizontalAlignment', 'right');
description = 'Log BR  priors';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.15], 'HorizontalAlignment', 'right');
text('String', 'CMSSM, \mu>0' ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.05], 'HorizontalAlignment', 'right');
colormap(chosencolor);
figure(   5)
subplot(    2,    2,    1);
box on; axis square; hold on;
bestfit = [    0.21243E+01    0.60887E+02];
mean = [    0.21835E+01    0.67093E+02];
refpoint = [    0.27520E+01    0.33800E+02];
tmp = load (fullfile('../../plot_data/test/','test_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_1_2_cont'));
scpts = load('test/test_single.txt');
scatter(scpts(:,   4),scpts(:,   3),'r', 'o', 'filled','SizeData', 15);
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','log[<\sigma v> (10^{27} cm^3 s^{-1})]','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('m_\chi (GeV)','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
description = 'Uniformly weighted samples';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.25], 'HorizontalAlignment', 'right');
description = 'Log BR  priors';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.15], 'HorizontalAlignment', 'right');
text('String', 'CMSSM, \mu>0' ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.05], 'HorizontalAlignment', 'right');
figure(   1)
subplot(    2,    2,    2);
box on; axis square; hold on;
bestfit = [    0.21507E-02    0.60887E+02];
mean = [    0.47468E-01    0.67093E+02];
refpoint = [    0.92400E+00    0.33800E+02];
tmp = load (fullfile('../../plot_data/test/','test_p3.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_1_3_cont'));
pdfpts=load (fullfile('../../plot_data/test/','test_2D_1_3_marg'));
contourf(x1,x2,pdfpts,64);
[C h] = contour(x1,x2,pdfpts,cnt(           1,:),l1);set(h,'LineWidth',lw1);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','BR_{b \bar{b}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('m_\chi (GeV)','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = modrhot;
description = 'Posterior pdf';
colormap(chosencolor);
figure(   4)
subplot(    2,    2,    2);
box on; axis square; hold on;
bestfit = [    0.21507E-02    0.60887E+02];
mean = [    0.47468E-01    0.67093E+02];
refpoint = [    0.92400E+00    0.33800E+02];
tmp = load (fullfile('../../plot_data/test/','test_p3.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_1_3_cont'));
profpts = load (fullfile('../../plot_data/test/','test_2D_1_3_profl'));
ah=imagesc(x1,x2,profpts);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','BR_{b \bar{b}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('m_\chi (GeV)','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = likegree;
description = 'Profile likelihood';
colormap(chosencolor);
figure(   5)
subplot(    2,    2,    2);
box on; axis square; hold on;
bestfit = [    0.21507E-02    0.60887E+02];
mean = [    0.47468E-01    0.67093E+02];
refpoint = [    0.92400E+00    0.33800E+02];
tmp = load (fullfile('../../plot_data/test/','test_p3.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_1_3_cont'));
scpts = load('test/test_single.txt');
scatter(scpts(:,   5),scpts(:,   3),'r', 'o', 'filled','SizeData', 15);
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','BR_{b \bar{b}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('m_\chi (GeV)','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
description = 'Uniformly weighted samples';
figure(   1)
subplot(    2,    2,    3);
box on; axis square; hold on;
bestfit = [    0.21243E+01    0.72332E+00];
mean = [    0.21835E+01    0.62796E+00];
refpoint = [    0.27520E+01    0.76000E-01];
tmp = load (fullfile('../../plot_data/test/','test_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p20.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_20_2_cont'));
pdfpts=load (fullfile('../../plot_data/test/','test_2D_20_2_marg'));
contourf(x1,x2,pdfpts,64);
[C h] = contour(x1,x2,pdfpts,cnt(           1,:),l1);set(h,'LineWidth',lw1);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','log[<\sigma v> (10^{27} cm^3 s^{-1})]','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('BR_{\tau^+\tau^-}','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = modrhot;
description = 'Posterior pdf';
colormap(chosencolor);
figure(   4)
subplot(    2,    2,    3);
box on; axis square; hold on;
bestfit = [    0.21243E+01    0.72332E+00];
mean = [    0.21835E+01    0.62796E+00];
refpoint = [    0.27520E+01    0.76000E-01];
tmp = load (fullfile('../../plot_data/test/','test_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p20.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_20_2_cont'));
profpts = load (fullfile('../../plot_data/test/','test_2D_20_2_profl'));
ah=imagesc(x1,x2,profpts);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','log[<\sigma v> (10^{27} cm^3 s^{-1})]','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('BR_{\tau^+\tau^-}','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = likegree;
description = 'Profile likelihood';
colormap(chosencolor);
figure(   5)
subplot(    2,    2,    3);
box on; axis square; hold on;
bestfit = [    0.21243E+01    0.72332E+00];
mean = [    0.21835E+01    0.62796E+00];
refpoint = [    0.27520E+01    0.76000E-01];
tmp = load (fullfile('../../plot_data/test/','test_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p20.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_20_2_cont'));
scpts = load('test/test_single.txt');
scatter(scpts(:,   4),scpts(:,  22),'r', 'o', 'filled','SizeData', 15);
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','log[<\sigma v> (10^{27} cm^3 s^{-1})]','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('BR_{\tau^+\tau^-}','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
description = 'Uniformly weighted samples';
figure(   1)
subplot(    2,    2,    4);
box on; axis square; hold on;
bestfit = [    0.21507E-02    0.72332E+00];
mean = [    0.47468E-01    0.62796E+00];
refpoint = [    0.92400E+00    0.76000E-01];
tmp = load (fullfile('../../plot_data/test/','test_p3.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p20.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_20_3_cont'));
pdfpts=load (fullfile('../../plot_data/test/','test_2D_20_3_marg'));
contourf(x1,x2,pdfpts,64);
[C h] = contour(x1,x2,pdfpts,cnt(           1,:),l1);set(h,'LineWidth',lw1);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','BR_{b \bar{b}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('BR_{\tau^+\tau^-}','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = modrhot;
description = 'Posterior pdf';
colormap(chosencolor);
figure(   4)
subplot(    2,    2,    4);
box on; axis square; hold on;
bestfit = [    0.21507E-02    0.72332E+00];
mean = [    0.47468E-01    0.62796E+00];
refpoint = [    0.92400E+00    0.76000E-01];
tmp = load (fullfile('../../plot_data/test/','test_p3.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p20.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_20_3_cont'));
profpts = load (fullfile('../../plot_data/test/','test_2D_20_3_profl'));
ah=imagesc(x1,x2,profpts);
shading flat;
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','BR_{b \bar{b}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('BR_{\tau^+\tau^-}','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
chosencolor = likegree;
description = 'Profile likelihood';
colormap(chosencolor);
figure(   5)
subplot(    2,    2,    4);
box on; axis square; hold on;
bestfit = [    0.21507E-02    0.72332E+00];
mean = [    0.47468E-01    0.62796E+00];
refpoint = [    0.92400E+00    0.76000E-01];
tmp = load (fullfile('../../plot_data/test/','test_p3.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/test/','test_p20.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/test/','test_2D_20_3_cont'));
scpts = load('test/test_single.txt');
scatter(scpts(:,   5),scpts(:,  22),'r', 'o', 'filled','SizeData', 15);
hold on;
text('String', '\it{SuperBayeS v 1.36}','FontSize',names_fontsize, 'Units', 'normalized', 'Position', [1.0 1.05], 'HorizontalAlignment', 'right');
 text('String','BR_{b \bar{b}}','FontSize',lab_fontsize, 'Units', 'normalized', 'Position', [0.5 -0.18], 'HorizontalAlignment', 'center');
ylabel('BR_{\tau^+\tau^-}','FontSize',lab_fontsize);
axis tight;
set(gca,'Layer','top','FontSize',axes_fontsize,  'LineWidth', lwaxes);
set(gca,'ticklength',2*get(gca,'ticklength'));
plot(bestfit(1), bestfit(2),'x','MarkerSize',bfsize,'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(bestfit(1), bestfit(2),'o','MarkerSize',bfsize,'MarkerFaceColor', 'none', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf);
plot(mean(1), mean(2),'o','MarkerSize',bfsize/2,'MarkerFaceColor','k', 'MarkerEdgeColor','k', 'LineWidth', bflw/redf*0.75);
plot(refpoint(1), refpoint(2),'d','MarkerSize',8,'MarkerFaceColor','r', 'MarkerEdgeColor','y', 'LineWidth', 2.0);
description = 'Uniformly weighted samples';
figure(   1)
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPosition',[ 0 0 20 20]);
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperOrientation', 'portrait'); 
print -dpsc2 ../../output_files/test/test_2D_marg.ps;
figure(   4)
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPosition',[ 0 0 20 20]);
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperOrientation', 'portrait'); 
print -dpsc2 ../../output_files/test/test_2D_profl.ps;
figure(   5)
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPosition',[ 0 0 20 20]);
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperOrientation', 'portrait'); 
print -dpsc2 ../../output_files/test/test_2D_single.ps;
