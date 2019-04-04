% Provides Matlab-like front-end to Octave ODE solver. But
% cannot use lsode, and hence this, recursively.  Used by
% MMburst.m, PIG.m, PIGExample.m, PIGExplore.m
% AJR, 4 Apr 2019
%{
\begin{matlab}
%}
function [ts,xs] = odeOct(dxdt,tSpan,x0)
    if length(tSpan)>2, ts = tSpan;
    else ts = linspace(tSpan(1),tSpan(end),21);
    end
    % mimic ode45 and ode23, but much slower for non-PI
    lsode_options('integration method','non-stiff');
    xs = lsode(@(x,t) dxdt(t,x),x0,ts);
end
%{
\end{matlab}
%}
