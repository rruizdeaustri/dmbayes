
tmp = load('../../chains/bb_grid_1.txt');
like = tmp(:,2); m = tmp(:,3);sigma = tmp(:,4);
% scatter3(m, sigma, like, 'filled')
% xlabel('mchi (GeV)')
% ylabel('log10 sigma v 10^{27}')
% zlabel('-log like')
% set(gca, 'ZLim', [600 800])
% %m = tmp(9:18, 3);
% %sigma = tmp(1:10, 4);
% like=reshape(tmp(:,2), 99, 99);

a = 100:40:300;
C = tmp(1:36,2);
b = 2.0:0.2:3.0;
D = reshape(C,6,6);
 cf =  min(min(D)):20: min(min(D))+10000;
  contourf(a,b,D')
  colorbar
  xlabel('mchi (GeV)')
ylabel('log10 sigma v 10^{27}')

  print -dpsc2 GridTest.ps 