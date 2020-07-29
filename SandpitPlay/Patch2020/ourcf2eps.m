function ourcf2eps(fileName)
% If uNFC3x6L variable exists in the Workspace (any value),
% then outputs the current figure to fileName.eps in the
% local Figs folder. If uNFC3x6L does not exist, then
% nothing happens.  Use 'clear uNFC3x6L' to remove variable
% if necessary.    AJR, 29 July 2020
global uNFC3x6L
if exist('uNFC3x6L')
  set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 14 10] ...
      ,'renderer','Painters')
  disp(['Making file Figs/' fileName '.eps'])
  print('-depsc2',['Figs/' fileName])
end
end
