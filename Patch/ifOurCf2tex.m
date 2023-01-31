function ifOurCf2tex(fileName)
% If global OurCf2eps variable exists with non-zero/true
% value, then outputs the current figure to fileName.tex in
% the local Figs folder, otherwise nothing happens.  To
% be used with \usepackage{pgfplots} \pgfplotsset{compat=newest}
% AJR, 31 Jan 2023
global OurCf2eps
if ~isempty(OurCf2eps) && OurCf2eps
    ffn = ['Figs/' fileName];
    disp(['***** Making file ' ffn '.tex  for pgfplots'])
    matlab2tikz([ffn '.tex'],'showInfo',false ...
    ,'noSize',true,'parseStrings',false,'showWarnings',false ...
    ,'extraCode',['\tikzsetnextfilename{' ffn '}'] ...
    ,'extraAxisOptions','\extraAxisOptions' )
end
end
