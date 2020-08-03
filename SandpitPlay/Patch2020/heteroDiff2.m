% Computes the time derivatives of heterogeneous diffusion
% in 2D on patches.  Adapted from 1D heterogeneous diffusion.
% JEB & AJR, May 2020 -- July 2020 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff2()}: heterogeneous diffusion}
\label{sec:heteroDiff2}
This function codes the lattice heterogeneous diffusion
inside the patches.  For 6D input arrays~\verb|u|, \verb|x|,
and~\verb|y| (via edge-value interpolation of
\verb|patchSmooth2|, \cref{sec:patchSmooth2}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
two 2D array of diffusivities,~$c^x_{ij}$ and~$c^y_{ij}$,
have previously been stored in struct~\verb|patches|. 
Although u and ut are 6D arrays, it seems we can just refer
to them as 3D since computations are independent of the last
four dimensions. 
\begin{matlab}
%}
function ut = heteroDiff2(t,u,x,y)
  global patches
  dx = diff(x(2:3));  % x space step
  dy = diff(y(2:3));  % y space step
  ix = 2:size(u,1)-1; % x interior points in a patch
  iy = 2:size(u,2)-1; % y interior points in a patch
  ut = nan(size(u)); % preallocate output array
  ut(ix,iy,:) = diff(patches.cx.*diff(u(:,iy,:),1),1)/dx^2 ...
               +diff(patches.cy.*diff(u(ix,:,:),1,2),1,2)/dy^2; 
end% function
%{
\end{matlab}
%}
