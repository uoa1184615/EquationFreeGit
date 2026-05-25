% time-homogenisation by projective integration of a rapidly
% fluctuating ODE.  AJR, 26 May 2026
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{timeHomogenise}: example
time-homogenisation by projective integration of a rapidly
fluctuating ODE}
\label{sec:ExPIRFode}
\localtableofcontents

Any one of three Projective Integration schemes are invoked
to integrate an ordinary differential equation which has a
rapidly fluctuating component (\cref{FRFode}).  Here
\begin{equation}
\de tu = (1-t/2)u+\frac1\tau u\cos(2\pi t/\tau),\quad 
u(0)=1\,,\quad
\text{over }0\leq t\leq 6\,.
\label{ERFode}
\end{equation}
The fluctuation, of microscale period~\(\tau\), is of large
strength~\(1/\tau\) in order for the fluctuations to have a
strong effect on the solutions (\cref{FRFode}). The trick to
this application is to sub-sample the micro-burst to one
data-point each microscale period (aka the micro-period
return map): then the projective integration functions just
operates on the `smooth' return-map data.
\begin{figure}
\centering
\caption{\label{FRFode}Time homogenisation of rapidly
fluctuating \ode~\cref{ERFode} by projective integration of
short burst: blue, the macroscale homogenisation solution; 
yellow, the micro-scale bursts of rapid fluctuations.}
\input{../ProjInt/Figs/timeHomogenise}
\end{figure}


Set the period~\(\tau\) as global variable~\verb|mPeriod|. 
Also (optional) define and reset global cells in order to
store the details of the micro-bursts \emph{before}
sub-sampling.
\begin{matlab}
%}
global ums tms mPeriod
ums={}; tms={}; 
%{
\end{matlab}
Select the projective integration method with~\verb|pirk|,
set macro-time-steps to be one (low accuracy), and set
micro-period~\(\tau\) here to be relatively big so we can
see the microscale fluctuations in the plots
\begin{matlab}
%}
pirk = 4    % select 1=PIG or 2nd or 4th order method
ts = 0:6    % set macroscale time steps
mPeriod = 0.1   % microscale period of fluctuations
%{
\end{matlab}
Execute selected projective integration scheme using
\verb|timeHomogeniseBurst| function.  Use a burst of length
just \emph{one period} long because here there are no
exp-decaying transients to overcome, and because one period
is sufficient for the projections.
\begin{matlab}
%}
switch pirk
case 2, u = PIRK2(@timeHomogeniseBurst, ts, 1, mPeriod);
case 4, u = PIRK4(@timeHomogeniseBurst, ts, 1, mPeriod);
case 1, mPeriod = 0.01    % shorter microscale
    [ts,u] = PIG(@ode23 ...
        ,@(t0, u0) timeHomogeniseBurst(t0,u0,mPeriod) ...
        ,ts([1 end]),1);
end
%{
\end{matlab}
Plot the PI macroscale solution (blue).  Add plots of the
stored microscale bursts (dark yellow).
\begin{matlab}
%}
figure(1), plot(ts,u,'o:')
title('Projective integration of rapidly fluctuating ODE')
xlabel('time $t$'), ylabel('$u(t)$')
hold on
arrayfun(@(x,y) plot(x{:},y{:},'Color',[0.71 0.54 0]) ...
    ,tms([1 2:pirk:end]),ums([1 2:pirk:end]))
hold off
%{
\end{matlab}

%Optionally output graphic for \verb|pgfplots|
%\begin{matlab}
%%}
%matlab2tikz('Figs/timeHomogenise.tex'...
%,'showInfo',false,'noSize',true...
%,'parseStrings',false,'showWarnings',false ...
%);
%%,'extraCode','\tikzsetnextfilename{xxxx}' ...
%%,'extraAxisOptions',string ...
%%{
%\end{matlab}


\paragraph{Code a burst of a rapidly fluctuating ODE}
Integrate a burst of length~\verb|bT| of the
\ode~\cref{ERFode} adapted from \cite{Bunder2013a}. 
\begin{matlab}
%}
function [ts, us] = timeHomogeniseBurst(ti, ui, bT) 
    global ums tms mPeriod
%{
\end{matlab}
Code \ode\ in internal function~\verb|dudt| with variable
\(u=\verb|u(1)|\).  For fluctuating ODE need the phase of
micro-fluctuation to be synchronised with each burst.
\begin{matlab}
%}
    dudt = @(t,u) (1-t/2)*u+u*cos(2*pi/mPeriod*(t-ti))/mPeriod;
%{
\end{matlab}
Use \verb|rk2Int| for the burst because it empowers us to
force that a burst is an integral number of micro-periods.
Ensure there are at least~10 time-steps per micro-period
(\verb|rk2Int| enforces at least~10).
\begin{matlab}
%}
    nPeriods = ceil(bT/mPeriod);
    ts = ti+linspace(0,nPeriods*mPeriod,nPeriods*10+1);
    [ts, us] = rk2Int(dudt, ts, ui);
%{
\end{matlab}
We then sub-sample the micro-burst at one phase in each
period.   To avoid discarding sub-period details, here store
sub-period details in global cell-arrays.  Since
\verb|rkInt| may increase the number of (equal-spaced)
micro-steps until accurate, we have to recompute the number
of micro-steps per micro-period.
\begin{matlab}
%}
    tms{end+1}=ts;  ums{end+1}=us;  
    mSteps = (length(ts)-1)/nPeriods;
    assert(round(mSteps)==mSteps,"not integer micro-steps")
    ts = ts(1:mSteps:end);
    us = us(1:mSteps:end);
end
%{
\end{matlab}
%}
