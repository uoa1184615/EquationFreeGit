% Simulates a burst of the system linked to by the
% configuration of patches. Used by ??.m
% AJR, 4 Apr 2019 -- 21 Oct 2022
%!TEX root = doc.tex
%{
\subsection{\texttt{heteroBurstF()}: a burst of
heterogeneous diffusion}
\label{sec:heteroBurstF}
This code integrates in time the derivatives computed by
\verb|heteroDiff| from within the patch coupling of
\verb|patchSys1|.  Try~\verb|ode23|, although \verb|ode45|
may give smoother results. Sample every period of the
microscale time fluctuations (or, at least, close to the
period).
\begin{matlab}
%}
function [ts, ucts] = heteroBurstF(ti, ui, bT) 
  global microTimePeriod
  [ts,ucts] = ode45( @patchSys1,ti+(0:microTimePeriod:bT),ui(:));
end
%{
\end{matlab}
%}
