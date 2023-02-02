% Test spectral, finite-width, and divided difference
% interpolation in 2D with many random parameters. 
% Including edgy interpolation. AJR, Nov 2018 -- 2 Feb 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{patchEdgeInt2test}: tests 2D patch coupling}
\label{sec:patchEdgeInt2test}

A script to test the spectral, finite-order, and divided
difference, polynomial interpolation of function
\verb|patchEdgeInt2()|. Tests one or several variables,
normal grids, and also tests centre and edge interpolation. 
But does not yet test staggered grids, core averaging, etc
as they are not yet implemented.

Start by establishing global data struct for the range of
various cases.  Choose a number of realisations for every
type.
\begin{matlab}
%}
clear all, close all
global patches
nRealise = 19
%{
\end{matlab}







\subsubsection{Check divided difference interpolation}
Check over various types and orders of interpolation,
numbers of patches, random domain lengths, random ratios,
and randomised distribution of patches. (The \verb|@sin| is
a dummy.) 
\begin{matlab}
%}
maxErrors=[];
for realisation = 1:nRealise
    Lx = 1+3*rand; Ly = 1+3*rand;
    Domain = [0 Lx 0 Ly]-[rand*[1 1] rand*[1 1]]
    nSubP = 1+2*randi(3,1,2)
    nPatch = randi([3 8],1,2)
    dx = [Lx Ly]./nPatch./nSubP.*rand(1,2)/2
    ordCC = 2*randi([1 4])
    edgyInt = (rand>0.5)
    configPatches2(@sin,Domain,'equispace' ...
        ,nPatch,ordCC,dx,nSubP,'EdgyInt',edgyInt);
%{
\end{matlab}
Second, displace patches to a random non-uniform spacing.
\begin{matlab}
%}
Hx = diff(patches.x(1,1,:,:,1:2,1));
patches.x = patches.x+0.8*Hx*(rand(1,1,1,1,nPatch(1),1)-0.5);
Hx = squeeze( diff(patches.x(1,1,:,:,:,1)) );% for information only
Hy = diff(patches.y(1,1,:,:,:,1:2));
patches.y = patches.y+0.8*Hy*(rand(1,1,1,1,1,nPatch(2))-0.5);
Hy = squeeze( diff(patches.y(1,1,:,:,1,:)) );% for information only
%{
\end{matlab}



\paragraph{Check multiple fields simultaneously} Set
profiles to be various powers of~\(x\) and~\(y\), \verb|ps|
and~\verb|qs|, and store as different `variables' at each
point. First, limit the order of test polynomials by the
order of interpolation and by the number of patches.
\begin{matlab}
%}
    ox=min(ordCC,nPatch(1)-1); 
    oy=min(ordCC,nPatch(2)-1);
    [ps,qs]=ndgrid(0:ox,0:oy);
    ps=reshape(ps,1,1,[]);
    qs=reshape(qs,1,1,[]);
    cs=2*rand(size(ps))-1;
    u0=cs.*patches.x.^ps.*patches.y.^qs;
%    sizeu0=size(u0)
%{
\end{matlab}
Then evaluate the interpolation
\begin{matlab}
%}
    u=u0; u([1 end],:)=inf; u(:,[1 end],:)=inf;
    ui=patchEdgeInt2(u(:));
%{
\end{matlab}
All patches should have zero error: but need to either in
\verb|patchEdgeInt2| comment out \verb|NaN| assignment of
boundary values, or not test the two extreme patches here,
or add code to omit NaNs here.  High-order interpolation
seems to be more affected by round-off so relax error size.
\begin{matlab}
%}
    error = ui-u0;
%    hist(log10(abs(error(abs(error)>1e-20))),-20:-8)
%    xlabel('log10 error')
%    pause(0.1)%??
    maxError=max(abs(error(:)))
    maxErrors=[maxErrors maxError];
    assert(maxError<1e-9 ...
    ,'failed divided difference interpolation')
    disp('***** This divided difference test passed')
%{
\end{matlab}

End the for-loops over various parameters.
\begin{matlab}
%}
end% for realisation
maxMaxErrorDividedDiffs = max(maxErrors)
disp('***** Passed all divided difference interpolation')
pause(1)
%{
\end{matlab}








\subsubsection{Test standard spectral interpolation}
Test over various numbers of patches, random domain lengths
and random ratios.  Try realisations of random tests.
\begin{matlab}
%}
for realisation=1:nRealise
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
nPatch = randi([3 6],1,2)
edgyInt = (rand>0.5)
configPatches2(@sin,[0 Lx 0 Ly],nan ...
    ,nPatch,0,ratios,nSubP,'EdgyInt',edgyInt);
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
error=u-u0;
assert(all(~isnan(error(:))),'found nans in the error!')
normError=norm(error(:))
assert(normError<1e-12, '2D spectral interpolation failed')
disp('***** This spectral test passed')
%{
\end{matlab}

End the for-loop over realisations
\begin{matlab}
%}
end
disp('***** All the spectral tests passed')
pause(1)
%{
\end{matlab}









\subsubsection{Check standard finite width interpolation}
Check over various types and orders of interpolation,
numbers of patches, random domain lengths and random ratios.
(The \verb|@sin| is a dummy.)
\begin{matlab}
%}
for realisations=1:nRealise
    ordCC=2*randi([1 4])
    nPatch=ordCC+randi([2 4],1,2)
    edgyInt = (rand>0.5)
    nSubP = 1+2*randi(3,1,2)
    Domain=5*[-rand rand -rand rand]
    ratios=0.5*rand(1,2)
    configPatches2(@sin,Domain,nan ...
        ,nPatch,ordCC,ratios,nSubP,'EdgyInt',edgyInt);
%{
\end{matlab}

\paragraph{Check multiple fields simultaneously} Set
profiles to be various powers of~\(x\), \verb|ps|, and store
as different `variables' at each point.  
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
error tolerance of~\(10^{-8}\) because of the size of the
domain and the high order of interpolation.
\begin{matlab}
%}
    I=ordCC/2+1:nPatch(1)-ordCC/2;
    J=ordCC/2+1:nPatch(2)-ordCC/2;
    error=ui(:,:,:,:,I,J)-u0(:,:,:,:,I,J);
    assert(all(~isnan(error(:))),'found nans in the error!')
    normError=norm(error(:))
    assert(normError<5e-8 ...
    ,'failed finite stencil interpolation')
    disp('***** This finite stencil test passed')
%{
\end{matlab}

End the for-loops over various parameters.
\begin{matlab}
%}
end %for realisations
disp('***** Passed all standard polynomial interpolation')
%{
\end{matlab}







\subsubsection{Finished}
If no error messages, then all OK.
\begin{matlab}
%}
disp('***** All the interpolation tests successful')
%{
\end{matlab}
%}