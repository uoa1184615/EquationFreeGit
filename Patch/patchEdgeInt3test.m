% Test inter-patch interpolation in 3D with many random parameters
% AJR, Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{patchEdgeInt3test}: tests 3D patch coupling}
\label{sec:patchEdgeInt3test}

A script to test the spectral and finite-order polynomial
interpolation of function \verb|patchEdgeInt3()|. Tests one
or several variables,  normal grids, and also tests centre
and edge interpolation.  But does not yet test staggered
grids, core averaging, etc.

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
Lx = 1+3*rand, Ly = 1+3*rand, Lz = 1+3*rand
nSubP = 1+2*randi(3,1,3)
ratios = 0.5*rand(1,3)
nPatch = randi([3 6],1,3)
edgyInt = (rand>0.5)
configPatches3(@sin,[0 Lx 0 Ly 0 Lz],nan,nPatch,0,ratios,nSubP ...
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
[nz,Nz]=size(squeeze(patches.z));
u0=nan(nx,ny,nz,nV,1,Nx,Ny,Nz);
for iV=1:nV
  ks=floor( rand(1,3).*floor(([Nx Ny Nz]+1)/2) )
  phis = pi*rand(1,3).*(2*ks~=[Nx Ny Nz])
  % generate 8D array via auto-replication
  u0(:,:,:,iV,1,:,:,:)=sin(2*pi*ks(1)*patches.x/Lx+phis(1)) ...
                     .*sin(2*pi*ks(2)*patches.y/Ly+phis(2)) ...
                     .*sin(2*pi*ks(3)*patches.z/Lz+phis(3));
end
%{
\end{matlab}
Copy and NaN the edges, then interpolate
\begin{matlab}
%}
u=u0; 
u([1 end],:,:,:)=nan; u(:,[1 end],:,:)=nan; u(:,:,[1 end],:)=nan;
u=patchEdgeInt3(u(:));
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
assert(normerr<1e-12, '3D interpolation failed')
%{
\end{matlab}

End the for-loop over realisations
\begin{matlab}
%}
disp('***** This test passed')
end
disp('If you read this, then all spectral tests passed')
%{
\end{matlab}









\subsubsection{Check polynomial finite width interpolation}
Check over various types and orders of interpolation,
numbers of patches, random domain lengths and random ratios.
(The \verb|@sin| is a dummy.)
\begin{matlab}
%}
for ordCC=2:2:6
for realisations=1:9
    ordCC=ordCC
    nPatch=ordCC+randi([1 4],1,3)
    edgyInt = (rand>0.5)
    nSubP = 1+2*randi(3,1,3)
    Domain=5*[-rand rand -rand rand -rand rand]
    ratios=0.5*rand(1,3)
    configPatches3(@sin,Domain,nan,nPatch,ordCC,ratios,nSubP ...
        ,'EdgyInt',edgyInt);
%{
\end{matlab}

\paragraph{Check multiple fields simultaneously}
Set profiles to be various powers of~\(x,y,z\), namely
\verb|ps|, \verb|qs|, \verb|rs|, and store as different
`variables' at each point.  
\begin{matlab}
%}
    [ps,qs,rs]=meshgrid(0:ordCC);
    ps=reshape(ps,1,1,1,[]); 
    qs=reshape(qs,1,1,1,[]);
    rs=reshape(rs,1,1,1,[]);
    u0=ones(size(ps)).*patches.x.^ps ...
        .*patches.y.^qs.*patches.z.^rs;
%{
\end{matlab}
Then evaluate the interpolation.
\begin{matlab}
%}
    ui=patchEdgeInt3(u0(:));
%{
\end{matlab}
The interior patches should have zero error. Appear to need
error tolerance of~\(10^{-8}\) because of the size of the
domain and the high order of interpolation.
\begin{matlab}
%}
    I=ordCC/2+1:nPatch(1)-ordCC/2;
    J=ordCC/2+1:nPatch(2)-ordCC/2;
    K=ordCC/2+1:nPatch(3)-ordCC/2;
    iError=ui(:,:,:,:,:,I,J,K)-u0(:,:,:,:,:,I,J,K);
    normError=norm(iError(~isnan(iError)))
    assert(normError<5e-8 ...
    ,'failed finite stencil polynomial interpolation')
%{
\end{matlab}

End the for-loops over various parameters.
\begin{matlab}
%}
end,end %order and realisations
disp('Passed all polynomial interpolation tests')
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