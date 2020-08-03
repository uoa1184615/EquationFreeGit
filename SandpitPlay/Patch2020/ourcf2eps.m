function ourcf2eps(fileName)
% If global uNFC3x6L variable exists with non-zero value,
% then outputs the current figure to fileName.eps in the
% local Figs folder, otherwise nothing happens.    
% AJR, 31 July 2020      
global uNFC3x6L
if ~isempty(uNFC3x6L) && uNFC3x6L
  set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 14 10] ...
      ,'renderer','Painters')
  disp(['***** Making file Figs/' fileName '.eps'])
  print('-depsc2',['Figs/' fileName])
end
end
