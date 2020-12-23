% Computes the time derivatives of heterogeneous diffusion
% in 1D on patches.  Used by homogenisationExample.m,
% homoDiffEdgy1.m  Optionally becomes Burgers PDE with
% heterogeneous advection. AJR, Apr 2019 -- Nov 2020
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
column vector of diffusivities~\(c_i\), and possibly
Burgers' advection coefficients~\(b_i\), have previously
been stored in struct~\verb|patches.cs|.
\begin{matlab}
%}
function ut = heteroDiff(t,u,patches)
  dx = diff(patches.x(2:3));   % space step
  i = 2:size(u,1)-1;   % interior points in a patch
  ut = nan+u;          % preallocate output array
  ut(i,:,:,:) = diff(patches.cs(:,1,:).*diff(u))/dx^2; 
  % possibly include heterogeneous Burgers' advection 
  if size(patches.cs,2)>1 % check for advection coeffs
      buu = patches.cs(:,2,:).*u.^2;
      ut(i,:) = ut(i,:)-(buu(i+1,:)-buu(i-1,:))/(dx*2);   
  end
end% function
%{
\end{matlab}
%}
