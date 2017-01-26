
tmp = load('../../chains/aa_grid_1.txt');
like = tmp(:,2); m = tmp(:,3);sigma = tmp(:,4);
% scatter3(m, sigma, like, 'filled')
% xlabel('mchi (GeV)')
% ylabel('log10 sigma v 10^{27}')
% zlabel('-log like')
% set(gca, 'ZLim', [600 800])
% %m = tmp(9:18, 3);
% %sigma = tmp(1:10, 4);
% like=reshape(tmp(:,2), 99, 99);

a = 100:10:190
 C = tmp(9:98,2)
 b = 2.1:0.1:2.9
 D = reshape(C,10,9)
 cf =  min(min(D)):2: min(min(D))+10;
  contourf(a,b,D', cf)
  colorbar
  xlabel('mchi (GeV)')
ylabel('log10 sigma v 10^{27}')

  print -dpsc2 GridTest.ps 