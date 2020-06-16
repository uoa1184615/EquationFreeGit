% Computes the time derivatives of heterogeneous diffusion
% in 1D on patches.  Used by homogenisationExample.m,
% ensembleAverageExample.m
% AJR, 4 Apr 2019 -- 7 Feb 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff()}: heterogeneous diffusion}
\label{sec:heteroDiff}
This function codes the lattice heterogeneous diffusion
inside the patches.  For 2D input arrays~\verb|u|
and~\verb|x| (via edge-value interpolation of
\verb|patchSmooth1|, \cref{sec:patchSmooth1}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
column vector (or possibly array) of diffusion
coefficients~\(c_i\) have previously been stored in
struct~\verb|patches|.
\begin{matlab}
%}
function ut = ensHeteroDiff(t,u,x)
  global patches
  dx = diff(x(2:3)); % space step
  i = 2:size(u,1)-1; % interior points in a patch
  m = size(u,3);
  j=1:m; jp=[2:m 1]; jm=[m 1:m-1]; % ensemble indices
  ut = nan(size(u)); % preallocate output array
  ut(i,:,j) = (patches.c(1,1,jp).*(u(i+1,:,jp)-u(i,:,j)) ...
               +patches.c(1,1,j).*(u(i-1,:,jm)-u(i,:,j)))/dx^2; 
end% function
%{
\end{matlab}
%}
