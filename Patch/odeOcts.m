% Provides Matlab-like front-end to Octave ODE solver. Uses
% non-stiff integrator as stiff ones are, surprisingly, far
% too slow. But cannot use lsode, and hence this function,
% recursively.  Used by configPatches1.m, configPatches2.m,
% ensembleAverageExample.m, homogenisationExample.m,
% waterWaveExample.m, wave2D.m, and so on. AJR, 17 Aug 2020
%{
\begin{matlab}
%}
function [ts,xs] = odeOcts(dxdt,tSpan,x0)
    if length(tSpan)>2, ts = tSpan;
    else ts = linspace(tSpan(1),tSpan(end),21)';
    end
    lsode_options('integration method','non-stiff');
    xs = lsode(@(x,t) dxdt(t,x),x0,ts);
end
%{
\end{matlab}
%}
