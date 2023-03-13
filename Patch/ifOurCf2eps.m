function ifOurCf2eps(fileName,sz)
% If global OurCf2eps variable exists with non-zero/true
% value, then outputs the current figure to fileName.eps in
% the local Figs folder, otherwise nothing happens. 
% However, eps format does not work with isosurfaces, so
% output those otherwise.  Include graphs with scale=0.8
% Optional: sz is of drawn size in cms (two element row).
% AJR, 11 Dec 2020  -- 13 Mar 2023
global OurCf2eps
if ~isempty(OurCf2eps) && OurCf2eps
  if nargin<2, sz=[17 12]; end
  a=gca;
  a.XLabel.Interpreter='latex';
  a.YLabel.Interpreter='latex';
  a.ZLabel.Interpreter='latex';
  set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 sz] ...
      ,'renderer','Painters')
  disp(['***** Making file Figs/' fileName '.eps'])
  print('-depsc2',['Figs/' fileName])
end
end
