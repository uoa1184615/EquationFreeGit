% Computes the time derivatives of heterogeneous diffusion
% in 2D on patches.  
% JEB May 2020, adapted from 1D heterogeneous diffusion (AJR, 4 Apr 2019)
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff2()}: heterogeneous diffusion}
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
function ut = heteroDiff2(t,u,x,y)
  global patches
  dx = diff(x(2:3)); % x space step
  dy = diff(y(2:3)); % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan(size(u)); % preallocate output array
  ut(ix,iy,:,:,:) = diff(patches.cx.*diff(u(:,iy,:,:,:),1),1)/dx^2 ...
      + diff(patches.cy.*diff(u(ix,:,:,:,:),1,2),1,2)/dy^2; 
end% function
%{
\end{matlab}
%}


% ix=2:size(u0,1)-1;
% iy = 2:size(u0,2)-1;
%  ut = nan(size(u0));
%   ut(ix,iy,:,:,:) = diff(patches.cx.*diff(u0(:,iy,:,:,:),1),1)/dx^2 ...
%       + diff(patches.cy.*diff(u0(ix,:,:,:,:),1,2),1,2)/dy^2; 