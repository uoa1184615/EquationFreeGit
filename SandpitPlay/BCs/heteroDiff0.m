% Computes the time derivatives of heterogeneous diffusion
% in 1D on patches.  Encodes Dirichlet zero BCS.  Optionally
% becomes Burgers PDE with heterogeneous advection. AJR, Apr
% 2019 -- Sep 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff0()}: heterogeneous diffusion}
\label{sec:heteroDiff0}

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
function ut = heteroDiff0(t,u,patches)
  dx = diff(patches.x(2:3));   % space step
  i = 2:size(u,1)-1;   % interior points in a patch
  ut = nan+u;          % preallocate output array
%{
\end{matlab}
The macroscale Dirichlet boundary conditions are zero at the
extreme edges of the two extreme patches.
\begin{matlab}
%}
  u( 1 ,:,:, 1 )=0; % left-edge of leftmost is zero
  u(end,:,:,end)=0; % right-edge of rightmost is zero
%{
\end{matlab}
Code the microscale diffusion, with possible advection.
\begin{matlab}
%}
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
