% Computes the time derivatives of forced heterogeneous
% diffusion in 1D on patches.  AJR, Apr 2019 -- May 2023
%!TEX root = doc.tex
%{
\subsection{\texttt{aHeteroDiffF()}: heterogeneous diffusion}
\label{sec:aHeteroDiffF}

This function codes the lattice heterogeneous diffusion
inside the patches with microscale boundary
conditions on the macroscale boundaries.  
\begin{matlab}
%}
function ut = aHeteroDiffF(t,u,patches)
  global Diri bcShft
%{
\end{matlab}
Two basic parameters, and initialise result array to NaNs.
\begin{matlab}
%}
  dx = diff(patches.x(2:3));   % space step
  i = 2:size(u,1)-1;   % interior points in a patch
  ut = nan+u;          % preallocate output array
%{
\end{matlab}
The macroscale Dirichlet boundary conditions are zero at the
extreme edges of the two extreme patches.
\begin{matlab}
%}
if Diri
  u( 1+bcShft ,:,:, 1 )=0; % left-edge of leftmost is zero
  u(end-bcShft,:,:,end)=0; % right-edge of rightmost is zero
else
  u( 1+bcShft ,:,:, 1 )=u( 2+bcShft ,:,:, 1 ); % left-edge of leftmost 
  u(end-bcShft,:,:,end)=u(end-1-bcShft,:,:,end); % right-edge of rightmost
end
%{
\end{matlab}
Code the microscale forced diffusion.
\begin{matlab}
%}
  ut(i,:,:,:) = diff(patches.cs(:,1,:).*diff(u))/dx^2; 
if Diri
  ut( 1+bcShft ,:,:, 1 )=0; % left-edge of leftmost is zero
  ut(end-bcShft,:,:,end)=0; % right-edge of rightmost is zero
else % ????????????
  ut( 1+bcShft ,:,:, 1 )=ut( 2+bcShft ,:,:, 1 ); % left-edge of leftmost 
  ut(end-bcShft,:,:,end)=ut(end-1-bcShft,:,:,end); % right-edge of end
end% function
%{
\end{matlab}
%}