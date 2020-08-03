% Test spectral interpolation in 2D with many random parameters
% AJR, Nov 2018 -- 17 Apr 2020 -- July 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{patchEdgeInt2test}: tests 2D patch coupling}
\label{sec:patchEdgeInt2test}

A script to test the spectral and finite-width interpolation
of function \verb|patchEdgeInt2()|. Tests one or several
variables,  normal grids, and also tests centre and edge
interpolation.  But does not yet test staggered grids, core
averaging, etc.

Start by establishing global data struct for the range of
various cases.
\begin{matlab}
%}
clear all, close all
global patches
%{
\end{matlab}




\subsubsection{Test standard spectral interpolation}
Test over various numbers of patches, random domain lengths
and random ratios.  Try 99 realisations of random tests.
\begin{matlab}
%}
for realisation=1:99
%{
\end{matlab}
Choose and configure random sized domains, random sub-patch
resolution, random size-ratios, random number of
periodic-patches, randomly edge or mid-patch interpolation.
\begin{matlab}
%}
Lx = 1+3*rand, Ly = 1+3*rand
nSubP = 1+2*randi(3,1,2)
ratios = rand(1,2)/2
nPatch = 2+randi(4,1,2)
edgyInt = (rand>0.5)
configPatches2(@sin,[0 Lx 0 Ly],nan,nPatch,0,ratios,nSubP ...
    ,'EdgyInt',edgyInt)
%{
\end{matlab}
Choose a random number of fields, then generate
trigonometric shape with random wavenumber and random phase
shift.  But if an even number of patches in either
direction, then do not test the highest wavenumber because
of aliasing problem.
\begin{matlab}
%}
nV=randi(3)
[nx,Nx]=size(squeeze(patches.x));
[ny,Ny]=size(squeeze(patches.y));
u0=nan(nx,ny,nV,1,Nx,Ny);
for iV=1:nV
  kx=randi([0 floor((nPatch(1)-1)/2)])
  ky=randi([0 floor((nPatch(2)-1)/2)])
  phix=pi*rand*(2*kx~=nPatch(1))
  phiy=pi*rand*(2*ky~=nPatch(2))
  % generate 6D array via auto-replication
  u0(:,:,iV,1,:,:)=sin(2*pi*kx*patches.x/Lx+phix) ...
                 .*sin(2*pi*ky*patches.y/Ly+phiy);
end
%{
\end{matlab}
Copy and NaN the edges, then interpolate
\begin{matlab}
%}
u=u0; u([1 end],:,:)=nan; u(:,[1 end],:)=nan;
u=patchEdgeInt2(u(:));
%{
\end{matlab}
Compute difference, ignoring the nans which should only be
in the corners. If there is an error in the interpolation,
then abort the script for checking: please record parameter
values and inform us.
\begin{matlab}
%}
err=u-u0;
normerr=norm(err(~isnan(err)))
assert(normerr<1e-12, '2D interpolation failed')
%{
\end{matlab}

End the for-loop over realisations
\begin{matlab}
%}
disp('***** This test passed')
end
%{
\end{matlab}









\subsubsection{Check standard finite width interpolation}
Check over various types and orders of interpolation,
numbers of patches, random domain lengths and random ratios.
(The \verb|@sin| is a dummy.)
\begin{matlab}
%}
for ordCC=2:2:8
for realisations=1:9
    ordCC=ordCC
    nPatch=ordCC+1+randi(3,1,2)
    edgyInt = (rand>0.5)
    nSubP = 1+2*randi(3,1,2)
    Domain=5*[-rand rand -rand rand]
    ratios=0.5*rand(1,2)
    configPatches2(@sin,Domain,nan,nPatch,ordCC,ratios,nSubP ...
        ,'EdgyInt',edgyInt);
%{
\end{matlab}

\paragraph{Check multiple fields simultaneously}
Set profiles to be various powers of~\(x\), \verb|ps|, and
store as different `variables' at each point.  
\begin{matlab}
%}
    [ps,qs]=meshgrid(0:ordCC);
    ps=reshape(ps,1,1,[]); qs=reshape(qs,1,1,[]);
    u0=ones(size(ps)).*patches.x.^ps.*patches.y.^qs;
%{
\end{matlab}
Then evaluate the interpolation.
\begin{matlab}
%}
    ui=patchEdgeInt2(u0(:));
%{
\end{matlab}
The interior patches should have zero error. Appear to need
error tolerance of~\(10^{-8}\) because of the siae of the
domain and the high order of interpolation.
\begin{matlab}
%}
    I=ordCC/2+1:nPatch(1)-ordCC/2;
    J=ordCC/2+1:nPatch(2)-ordCC/2;
    iError=ui(:,:,:,:,I,J)-u0(:,:,:,:,I,J);
    normError=norm(iError(~isnan(iError)))
    assert(normError<5e-8 ...
    ,'failed finite stencil interpolation')
%{
\end{matlab}

End the for-loops over various parameters.
\begin{matlab}
%}
end,end %order and realisations
disp('Passed all standard polynomial interpolation')
%{
\end{matlab}







\subsubsection{Finished}
If no error messages, then all OK.
\begin{matlab}
%}
disp('**** If you read this, then all the tests were successful')
%{
\end{matlab}
%}