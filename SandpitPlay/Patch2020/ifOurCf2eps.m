function ifOurCf2eps(fileName)
% If global OurCf2eps variable exists with non-zero value,
% then outputs the current figure to fileName.eps in the
% local Figs folder, otherwise nothing happens.    
% AJR, 7 Aug 2020      
global OurCf2eps
if ~isempty(OurCf2eps) && OurCf2eps
  set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 14 10] ...
      ,'renderer','Painters')
  disp(['***** Making file Figs/' fileName '.eps'])
  print('-depsc2',['Figs/' fileName])
end
end
