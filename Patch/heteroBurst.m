% Simulates a burst of the system linked to by the
% configuration of patches. Used by homogenisationExample.m,
% homoDiffEdgy1.m, and maybe homoLanLif1D.m
% AJR, 4 Apr 2019 -- Sep 2021
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroBurst()}: a burst of heterogeneous diffusion}
\label{sec:heteroBurst}
This code integrates in time the derivatives computed by
\verb|heteroDiff| from within the patch coupling of
\verb|patchSmooth1|.  Try~\verb|ode23| or \verb|rk2Int|,
although \verb|ode45| may give smoother results.
\begin{matlab}
%}
function [ts, ucts] = heteroBurst(ti, ui, bT) 
	if ~exist('OCTAVE_VERSION','builtin')
	[ts,ucts] = ode23( @patchSmooth1,[ti ti+bT],ui(:));
	else % octave version
	[ts,ucts] = rk2Int(@patchSmooth1,[ti ti+bT],ui(:));
	end
end
%{
\end{matlab}
%}
