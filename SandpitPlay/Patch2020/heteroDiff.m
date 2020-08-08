% Computes the time derivatives of heterogeneous diffusion
% in 1D on patches.  Used by homogenisationExample.m,
% ensembleAverageExample.m  Optionally becomes Burgers PDE
% with heterogeneous advection.
% AJR, 4 Apr 2019 -- Aug 2020
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
column vector, or possibly array of phase-shifts, of
diffusivities~\(c_i\) have previously been stored in
struct~\verb|patches.cs|.
\begin{matlab}
%}
function ut = heteroDiff(t,u,x)
  global patches
  dx = diff(x(2:3));   % space step
  i = 2:size(u,1)-1;   % interior points in a patch
  ut = nan(size(u)); % preallocate output array
  ut(i,:,:,:) = diff(patches.cs.*diff(u))/dx^2; 
  % if set, include heterogeneous Burgers' advection 
  if isfield(patches,'b') % check for advection coeffs
      buu = patches.b.*u.^2;
      ut(i,:) = ut(i,:)-(buu(i+1,:)-buu(i-1,:))/(dx*2);   
  end
end% function
%{
\end{matlab}
%}
