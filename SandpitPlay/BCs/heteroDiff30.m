% Computes the time derivatives of heterogeneous diffusion
% in 3D on patches with macroscale BCs.  AJR, Oct 2023 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroDiff30()}: heterogeneous diffusion}
\label{sec:heteroDiff30}

This function codes the lattice heterogeneous diffusion
inside the patches.  For 8D input array~\verb|u| (via
edge-value interpolation of \verb|patchEdgeInt3|, such as by
\verb|patchSys3|, \cref{sec:patchSys3}), computes the
time derivative~\cref{eq:HomogenisationExample} at each
point in the interior of a patch, output in~\verb|ut|.  The
three 3D array of diffusivities,~$c^x_{ijk}$, $c^y_{ijk}$
and~$c^z_{ijk}$, have previously been stored
in~\verb|patches.cs| (4+D). 

Supply patch information as a third argument (required by
parallel computation), or otherwise by a global variable.
\begin{matlab}
%}
function ut = heteroDiff30(t,u,patches)
  if nargin<3, global patches, end
%{
\end{matlab}
Microscale space-steps.  
\begin{matlab}
%}
  dx = diff(patches.x(2:3));  % x micro-scale step
  dy = diff(patches.y(2:3));  % y micro-scale step
  dz = diff(patches.z(2:3));  % z micro-scale step
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
  k = 2:size(u,3)-1; % y interior points in a patch
%{
\end{matlab}
The macroscale boundary conditions are Dirichlet zero at the
extreme edges of left-right extreme patches, and Neumann zero 
at extreme edges of top-bottom extreme patches.
\begin{matlab}
%}
  u( 1 ,:,:,:,:, 1 ,:,:)=0; % left-edge of leftmost is zero
  u(end,:,:,:,:,end,:,:)=0; % right-edge of rightmost is zero
  u(:, 1 ,:,:,:,:, 1 ,:)=u(:,  2  ,:,:,:,:, 1 ,:); % bottom-edge of bottommost
  u(:,end,:,:,:,:,end,:)=u(:,end-1,:,:,:,:,end,:); % top-edge of topmost
  u(:,:, 1 ,:,:,:,:, 1 )=0; % back-edge of rearmost is one
  u(:,:,end,:,:,:,:,end)=1; % front-edge of frontmost is zero
%{
\end{matlab}
Reserve storage and then assign interior patch values to the
heterogeneous diffusion time derivatives. Using \verb|nan+u|
appears quicker than \verb|nan(size(u),patches.codist)|
\begin{matlab}
%}
  ut = nan+u; % reserve storage
  ut(i,j,k,:,:,:,:,:) ...
  = diff(patches.cs(:,j,k,1,:).*diff(u(:,j,k,:,:,:,:,:),1),1)/dx^2 ...
   +diff(patches.cs(i,:,k,2,:).*diff(u(i,:,k,:,:,:,:,:),1,2),1,2)/dy^2 ...
   +diff(patches.cs(i,j,:,3,:).*diff(u(i,j,:,:,:,:,:,:),1,3),1,3)/dz^2; 
end% function
%{
\end{matlab}
%}
