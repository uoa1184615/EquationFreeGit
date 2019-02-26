% Test spectral interpolation in 2D with many random parameters
% AJR, Nov 2018
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt2test}: tests 2D spectral interpolation}
\label{sec:...}
\localtableofcontents
Try 99 realisations of random tests.
\begin{matlab}
%}
clear all, close all
global patches
for realisation=1:99
%{
\end{matlab}
Choose and configure random sized domains, random sub-patch resolution, random size-ratios, random number of periodic-patches.
\begin{matlab}
%}
Lx=1+3*rand, Ly=1+3*rand
nSubP=1+2*randi(3,1,2)
ratios=rand(1,2)/2
nPatch=2+randi(4,1,2)
configPatches2(@sin,[0 Lx 0 Ly],nan,nPatch,0,ratios,nSubP)
%{
\end{matlab}
Choose a random number of fields, then generate trigonometric shape with random wavenumber and random phase shift.
\begin{matlab}
%}
nV=randi(3)
[nx,Nx]=size(patches.x);
[ny,Ny]=size(patches.y);
u0s=nan(nx,ny,Nx,Ny,nV);
for iV=1:nV
  kx=randi([0 ceil((nPatch(1)-1)/2)])
  ky=randi([0 ceil((nPatch(2)-1)/2)])
  phix=pi*rand*(2*kx~=nPatch(1))
  phiy=pi*rand*(2*ky~=nPatch(2))
  % generate 2D array via auto-replication
  u0=sin(2*pi*kx*patches.x(:) /Lx+phix) ...
   .*sin(2*pi*ky*patches.y(:)'/Ly+phiy);
  % reshape into 4D array
  u0=reshape(u0,[nx Nx ny Ny]);
  u0=permute(u0,[1 3 2 4]);
  % store into 5D array
  u0s(:,:,:,:,iV)=u0;
end
%{
\end{matlab}
Copy and NaN the edges, then interpolate
\begin{matlab}
%}
u=u0s; u([1 end],:,:,:,:)=nan; u(:,[1 end],:,:,:)=nan;
u=patchEdgeInt2(u(:));
%{
\end{matlab}
If there is an error in the interpolation then abort the script for checking: record parameter values and inform.
\begin{matlab}
%}
err=u-u0s;
normerr=norm(err(:))
if normerr>1e-12, error('2D interpolation failed'), end
end
%{
\end{matlab}
%}