% Computes the time derivatives of heterogeneous diffusion
% in 1D on patches.  Used by homogenisationExample.m,
% homoDiffEdgy1.m  Optionally becomes Burgers PDE with
% heterogeneous advection. AJR, Apr 2019 -- Jun 2026
% AutoDiff prefers  ,'like',u  in preallocation.
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff()}: heterogeneous diffusion}
\label{sec:heteroDiff}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 2D input arrays~\verb|u|
and~\verb|x| (via edge-value interpolation of
\verb|patchSys1|, \cref{sec:patchSys1}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
column vector of diffusivities~\(c_i\), and possibly
Burgers' advection coefficients~\(b_i\), have previously
been stored in struct~\verb|patches.cs|.
\begin{matlab}
%}
function ut = heteroDiff(t,u,patches)
    if ~patches.periodic  % apply simple Dirichlet BCs
        u( 1 ,:,:, 1 )=0; % left-edge of leftmost is zero
        u(end,:,:,end)=0; % right-edge of rightmost is zero
    end%if not periodic

  dx = diff(patches.x(2:3));  % space step
  i = 2:size(u,1)-1;          % interior points in a patch
  ut = nan(size(u),'like',u); % preallocate output array, nan+u poss
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
