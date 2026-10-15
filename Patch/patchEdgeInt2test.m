% Test spectral, finite-width, and divided difference
% interpolation in 2D with many random parameters. 
% AJR, Nov 2018 -- 16 Oct 2026
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
clear all
global patches
nRealise = 20
%{
\end{matlab}







\subsubsection{Check divided difference interpolation}
\begin{matlab}
%}
fprintf('\n\n**** Check divided difference interpolation\n')
pause(1)
%{
\end{matlab}
Check over various types and orders of interpolation,
numbers of patches, random domain lengths, random ratios,
and randomised distribution of patches. (The \verb|@sin| is
a dummy.) 
\begin{matlab}
%}
maxErrors=[];
for realisation = 1:nRealise
    nEdge = randi(3,1,2)% =1,2, or 3
    edgyInt = (rand>0.5)
    Lx = 1+3*rand; Ly = 1+3*rand;
    xyLim = [0 Lx 0 Ly]-[rand*[1 1] rand*[1 1]]
    nSubP = nEdge.*( (2-edgyInt)*randi(3,1,2)+1+edgyInt )
    ordCC = 2*randi(4)
    nPatch = ordCC+randi(4,1,2)
    dx = [Lx Ly]./nPatch./nSubP.*rand(1,2)/2
    configPatches2(@sin,xyLim,'equispace',nPatch,ordCC ...
        ,dx,nSubP,'EdgyInt',edgyInt,'nEdge',nEdge);
%{
\end{matlab}
Second, displace patches to a random non-uniform spacing.
\begin{matlab}
%}
Hx = diff(patches.x(1,1,:,:,1:2,1));
patches.x = patches.x+0.9*Hx*(rand(1,1,1,1,nPatch(1),1)-0.5);
Hx = squeeze( diff(patches.x(1,1,:,:,:,1)) );% for information only
Hy = diff(patches.y(1,1,:,:,:,1:2));
patches.y = patches.y+0.9*Hy*(rand(1,1,1,1,1,nPatch(2))-0.5);
Hy = squeeze( diff(patches.y(1,1,:,:,1,:)) );% for information only
%{
\end{matlab}



\paragraph{Check multiple fields simultaneously} 
Set profiles to be products of all Chebyshev polynomials
in~\(x,y\), optionally multiplied by random coefficients,
and store as different `variables' at each point. Limit the
order of test polynomials by the order of interpolation and
by the number of patches.
\begin{matlab}
%}
    nox=min(ordCC+1,nPatch(1)); 
    noy=min(ordCC+1,nPatch(2));
    u0=nan([nSubP,nox*noy,1,nPatch]);
    % first fillin Chebyshev in x, constant in y
    xx=quantile(patches.x(:),[0 1]);
    xi=(patches.x-mean(xx))*2/diff(xx);
    assert(abs(1-max(abs(xi(:))))<1e-10,"x-Chebyshev scaling incorrect")
    u0(:,:,1,:,:,:)=1;
    u0(:,:,2,:,:,:)=repmat(xi,[1,nSubP(2),1,1,1,nPatch(2)]);
    for n=2:nox-1
        u0(:,:,n+1,:,:,:)=2*xi.*u0(:,:,n,:,:,:)-u0(:,:,n-1,:,:,:);
    end;
    % 2nd fillin by multiplying by y-recurrence
    yy=quantile(patches.y(:),[0 1]);
    xi=(patches.y-mean(yy))*2/diff(yy);
    assert(abs(1-max(abs(xi(:))))<1e-10,"y-Chebyshev scaling incorrect")
    n=1:nox;  u0(:,:,n+nox,:,:,:)=xi.*u0(:,:,n,:,:,:);
    for n=nox+1:(noy-1)*nox
        u0(:,:,n+nox,:,:,:)=2*xi.*u0(:,:,n,:,:,:)-u0(:,:,n-nox,:,:,:);
    end;
    u0=u0.*randn(1,1,nox*noy);
%{
\end{matlab}
Then evaluate the interpolation, setting edges to \verb|inf|
for error checking.
\begin{matlab}
%}
    u=u0; 
    u([1:nEdge(1)  end-nEdge(1)+1:end],:,:)=inf; 
    u(:,[1:nEdge(2)  end-nEdge(2)+1:end],:)=inf;
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
    hist(log10(abs(error(abs(error)>1e-20))),-20:-7)
    xlabel('log10 error'), pause(0.3)%??
    maxError=max(abs(error(:)))
    maxErrors=[maxErrors maxError];
    assert(maxError<1e-11 ...
    ,'failed divided difference interpolation')
    disp('*** This divided difference test passed')
%{
\end{matlab}

End the for-loops over various parameters.
\begin{matlab}
%}
end% for realisation
maxMaxErrorDividedDiffs = max(maxErrors)
disp('***** Passed all divided difference interpolation')
pause(1)
%assert(0,"stop here---only testing finite domain divided differences")
%{
\end{matlab}








\subsubsection{Test standard spectral interpolation}
\begin{matlab}
%}
fprintf('\n\n**** Test standard spectral interpolation\n')
pause(1)
%{
\end{matlab}
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
nEdge=randi(3,1,2)% =1,2, or 3
edgyInt = (rand>0.5)
Lx = 1+3*rand, Ly = 1+3*rand
xyLim = [0 Lx 0 Ly]-[rand*[1 1] rand*[1 1]]
nSubP = nEdge.*( (2-edgyInt)*randi(3,1,2)+1+edgyInt )
nPatch = randi([3 6],1,2)
dx = [Lx Ly]./nPatch./nSubP.*rand(1,2)/2
configPatches2(@sin,xyLim,'periodic',nPatch,0 ...
    ,dx,nSubP,'EdgyInt',edgyInt,'nEdge',nEdge);
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
Copy and \verb|nan| the edges, then interpolate
\begin{matlab}
%}
u=u0; 
u([1:nEdge(1)  end-nEdge(1)+1:end],:,:)=nan; 
u(:,[1:nEdge(2)  end-nEdge(2)+1:end],:)=nan;
u=patchEdgeInt2(u(:));
%{
\end{matlab}
Compute difference. If there is an error in the interpolation,
then abort the script for checking: please record parameter
values and inform us.
\begin{matlab}
%}
error = u-u0;
assert(all(~isnan(error(:))),'found nans in the error!')
hist(log10(abs(error(abs(error)>1e-20))),-20:-7)
xlabel('log10 error'), pause(0.3)%??
normError=norm(error(:))
assert(normError<1e-12, '2D spectral interpolation failed')
disp('*** This spectral test passed')
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









\subsubsection{Check polynomial finite width interpolation}
Check over various types and orders of interpolation,
numbers of patches, random domain lengths and random ratios.
(The \verb|@sin| is a dummy.)
\begin{matlab}
%}
for realisations=1:nRealise
    nEdge = randi(3,1,2)% =1,2, or 3
    edgyInt = (rand>0.5)
    nSubP = nEdge.*( (2-edgyInt)*randi(3,1,2)+1+edgyInt )
    ordCC = 2*randi(4)
    nPatch = ordCC+randi(4,1,2)
    xyLim=5*[-rand(1,2); rand(1,2)]
    dx = diff(xyLim)./nPatch./nSubP.*rand(1,2)/2
    configPatches2(@sin,xyLim,'periodic',nPatch,ordCC ...
        ,dx,nSubP,'EdgyInt',edgyInt,'nEdge',nEdge);
%{
\end{matlab}

\paragraph{Check multiple fields simultaneously} 
Set profiles to be products of all Chebyshev polynomials in~\(x,y\), 
optionally multiplied by random coefficients, and
store as different `variables' at each point.  
Limit the order of test polynomials by the
order of interpolation and by the number of patches.  
\begin{matlab}
%}
    nox=ordCC+1;	noy=ordCC+1;
    u0=nan([nSubP,nox*noy,1,nPatch]);
    % first fillin Chebyshev in x, constant in y
    xx=quantile(patches.x(:),[0 1]);
    xi=(patches.x-mean(xx))*2/diff(xx);
    assert(abs(1-max(abs(xi(:))))<1e-10,"x-Chebyshev scaling incorrect")
    u0(:,:,1,:,:,:)=1;
    u0(:,:,2,:,:,:)=repmat(xi,[1,nSubP(2),1,1,1,nPatch(2)]);
    for n=2:nox-1
        u0(:,:,n+1,:,:,:)=2*xi.*u0(:,:,n,:,:,:)-u0(:,:,n-1,:,:,:);
    end;
    % 2nd fillin by multiplying by y-recurrence
    yy=quantile(patches.y(:),[0 1]);
    xi=(patches.y-mean(yy))*2/diff(yy);
    assert(abs(1-max(abs(xi(:))))<1e-10,"y-Chebyshev scaling incorrect")
    n=1:nox;  u0(:,:,n+nox,:,:,:)=xi.*u0(:,:,n,:,:,:);
    for n=nox+1:(noy-1)*nox
        u0(:,:,n+nox,:,:,:)=2*xi.*u0(:,:,n,:,:,:)-u0(:,:,n-nox,:,:,:);
    end;
    u0=u0.*randn(1,1,nox*noy);
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
    hist(log10(abs(error(abs(error)>1e-20))),-20:-7)
    xlabel('log10 error'), pause(0.3)%??
    normError=norm(error(:))
    assert(normError<5e-9 ...
    ,'failed finite stencil interpolation')
    disp('*** This finite stencil test passed')
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
disp('******* All the interpolation tests successful')
%{
\end{matlab}
%}