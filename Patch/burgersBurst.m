% Simulates a burst in time of a microscale map that is
% applied on patches in space.  Used by BurgersExample.m
% AJR, 4 Apr 2019
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{burgerBurst()}: code a burst of the patch map}
\label{sec:burgerBurst}
\begin{matlab}
%}
function [ts, us] = burgersBurst(ti, ui, bT) 
%{
\end{matlab}
First find and set the number of microscale time-steps.
\begin{matlab}
%}
  global patches
  dt = diff(patches.x(2:3))^2/2;
  ndt = ceil(bT/dt -0.2);
  ts = ti+(0:ndt)'*dt;
%{
\end{matlab}
Use \verb|patchSmooth1()| (\cref{sec:patchSmooth1}) to apply
the microscale map over all time-steps in the burst. The
\verb|patchSmooth1()| interface provides the interpolated
edge-values of each patch.  Store the results in rows to be
consistent with \ode\ and projective integrators.
\begin{matlab}
%}
  us = nan(ndt+1,numel(ui)); 
  us(1,:) = reshape(ui,1,[]);
  for j = 1:ndt
    ui = patchSmooth1(ts(j),ui);
    us(j+1,:) = reshape(ui,1,[]);
  end
%{
\end{matlab}
Linearly interpolate (extrapolate) to get the field values
at the precise final time of the burst.  Then return.
\begin{matlab}
%}
  ts(ndt+1) = ti+bT;
  us(ndt+1,:) = us(ndt,:) ...
    + diff(ts(ndt:ndt+1))/dt*diff(us(ndt:ndt+1,:));
end
%{
\end{matlab}
%}
