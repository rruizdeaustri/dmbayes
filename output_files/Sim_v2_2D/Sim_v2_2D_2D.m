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
bestfit = [    0.29998E+01    0.10009E+03];
mean = [    0.29996E+01    0.10008E+03];
refpoint = [    0.23500E+01    0.16750E+03];
tmp = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_2D_1_2_cont'));
pdfpts=load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_2D_1_2_marg'));
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
description = 'Flat BR priors';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.15], 'HorizontalAlignment', 'right');
text('String', 'CMSSM, \mu>0' ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.05], 'HorizontalAlignment', 'right');
colormap(chosencolor);
figure(   4)
subplot(    2,    2,    1);
box on; axis square; hold on;
bestfit = [    0.29998E+01    0.10009E+03];
mean = [    0.29996E+01    0.10008E+03];
refpoint = [    0.23500E+01    0.16750E+03];
tmp = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_2D_1_2_cont'));
profpts = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_2D_1_2_profl'));
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
description = 'Flat BR priors';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.15], 'HorizontalAlignment', 'right');
text('String', 'CMSSM, \mu>0' ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.05], 'HorizontalAlignment', 'right');
colormap(chosencolor);
figure(   5)
subplot(    2,    2,    1);
box on; axis square; hold on;
bestfit = [    0.29998E+01    0.10009E+03];
mean = [    0.29996E+01    0.10008E+03];
refpoint = [    0.23500E+01    0.16750E+03];
tmp = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_p2.dat'));
 x1 = tmp(:,1);
tmp = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_p1.dat'));
 x2 = tmp(:,1);
cnt = load (fullfile('../../plot_data/Sim_v2_2D/','Sim_v2_2D_2D_1_2_cont'));
scpts = load('Sim_v2_2D/Sim_v2_2D_single.txt');
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
description = 'Flat BR priors';
text('String', description ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.15], 'HorizontalAlignment', 'right');
text('String', 'CMSSM, \mu>0' ,'FontSize',cbar_fontsize, 'Units', 'normalized', 'Position', [0.97 0.05], 'HorizontalAlignment', 'right');
figure(   1)
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPosition',[ 0 0 20 20]);
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperOrientation', 'portrait'); 
print -dpsc2 ../../output_files/Sim_v2_2D/Sim_v2_2D_2D_marg.ps;
figure(   4)
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPosition',[ 0 0 20 20]);
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperOrientation', 'portrait'); 
print -dpsc2 ../../output_files/Sim_v2_2D/Sim_v2_2D_2D_profl.ps;
figure(   5)
set(gcf, 'PaperUnits','centimeters');
set(gcf, 'PaperPosition',[ 0 0 20 20]);
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperOrientation', 'portrait'); 
print -dpsc2 ../../output_files/Sim_v2_2D/Sim_v2_2D_2D_single.ps;
