% Provides Matlab-like front-end to Octave ODE solver.
% Default lsode() will be a stiff integrator when necessary.
% But cannot use lsode, and hence this, recursively.  Used
% by configPatches1.m, configPatches2.m,
% ensembleAverageExample.m, homogenisationExample.m,
% waterWaveExample.m, wave2D.m. AJR, 4 Apr 2019
%{
\begin{matlab}
%}
function [ts,xs] = odeOcts(dxdt,tSpan,x0)
    if length(tSpan)>2, ts = tSpan;
    else ts = linspace(tSpan(1),tSpan(end),21)';
    end
    lsode_options('integration method','stiff');
    xs = lsode(@(x,t) dxdt(t,x),x0,ts);
end
%{
\end{matlab}
%}
