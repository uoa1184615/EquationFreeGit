% Test the spectral, finite-width, and divided difference 1D
% interpolation of patchEdgeInt1(), for one or several
% variables, including normal grid and staggered grid, and
% also edgy interpolation. All tests passed in 2019--2020,
% except not testing the zig-zag modes.
% AJR, 26 Sep 2018 -- Jan 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{patchEdgeInt1test}: test the 1D patch coupling}
\label{sec:patchEdgeInt1test}

A script to test the spectral and finite-order polynomial
interpolation of function \verb|patchEdgeInt1()|. Tests one
or several variables,  normal and staggered grids, and also
tests centre and edge interpolation.  But does not yet test
core averaging, nor divided differences on staggered, etc.

Start by establishing global data struct for the range of
various cases.
\begin{matlab}
%}
clear all
global patches
nSubP=5
%{
\end{matlab}






\subsubsection{Check divided difference interpolation}
But not yet implemented staggered grid version?? Check over
various types and orders of interpolation, numbers of
patches, random domain lengths, random ratios, and
randomised distribution of patches. (The \verb|@sin| is a
dummy.)
\begin{matlab}
%}
for edgyInt=[mod(nSubP,2)==0 true]
for ordCC=2:2:8
for nPatch=ordCC+(2:4)
    edgyInt=edgyInt
    ordCC=ordCC
    nPatch=nPatch
    Domain=5*[-rand rand]
    ratio=0.5*rand
    configPatches1(@sin,Domain,nan,nPatch,ordCC,ratio,nSubP ...
        ,'EdgyInt',edgyInt);
%{
\end{matlab}
Until \verb|configPatches1| updated, change the data
structure here: first, specify general divided differences.
\begin{matlab}
%}
patches.periodic=false;
%{
\end{matlab}
Second, displace patches to a random non-uniform spacing.
\begin{matlab}
%}
H = diff(patches.x(1,:,:,1:2));
patches.x = patches.x+0.8*H*(rand(1,1,1,nPatch)-0.5);
H = squeeze( diff(patches.x(1,:,:,:)) )% for information only
%{
\end{matlab}



\paragraph{Check multiple fields simultaneously}
Set profiles to be various powers of~\(x\), \verb|ps|, and
store as different `variables' at each point.  
\begin{matlab}
%}
    ps=1:ordCC
    cs=randn(size(ps))
    u0=patches.x.^ps.*cs+randn;
%{
\end{matlab}
Then evaluate the interpolation and squeeze the singleton
dimension of an `ensemble'.
\begin{matlab}
%}
    ui=patchEdgeInt1(u0(:));
    ui=squeeze(ui);
%{
\end{matlab}
All patches should have zero error: but need to either in
\verb|patchEdgeInt1| comment out \verb|NaN| assignment of
boundary values, or not test the two extreme patches here,
or add code to omit NaNs here.  High-order interpolation
seems to be more affected by round-off so relax error size.
\begin{matlab}
%}
    j=1:nPatch;
    iError=ui(:,:,j)-u0(:,:,j);
    hist(log10(abs(iError(abs(iError)>0))),-17:-10)
    xlabel('log10 iError'),pause(0.3)%??
    normError=norm(iError(:))
    assert(normError<5e-10 ...
    ,'failed divided difference interpolation')
%{
\end{matlab}

End the for-loops over various parameters.
\begin{matlab}
%}
end,end,end
disp('Passed all divided difference interpolation')
%error('Not testing any further')% temporary halt
%{
\end{matlab}










\subsubsection{Test standard spectral interpolation}
Test over various numbers of patches, random domain lengths
and random ratios.  Say do three realisations for each
number of patches.
\begin{matlab}
%}
nReal=3 
for nPatch=repmat(5:10,1,nReal) 
nPatch=nPatch
Len=10*rand
ratio=0.5*rand
configPatches1(@sin,[0,Len],nan,nPatch,0,ratio,nSubP ...
    ,'EdgyInt',rand<0.5); % random Edgy or not
if mod(nPatch,2)==0, disp('Avoiding highest wavenumber'), end
kMax=floor((nPatch-1)/2);
if patches.EdgyInt, i0=[2 nSubP-1], else i0=(nSubP+1)/2, end
patches.periodic=true; % temporary
%{
\end{matlab}

\paragraph{Test single field}
Set a profile, and evaluate the interpolation.
\begin{matlab}
%}
for k=-kMax:kMax
  u0=exp(1i*k*patches.x*2*pi/Len);
  ui=patchEdgeInt1(u0(:));
  normError=rms(ui(:)-u0(:));
  if abs(normError)>5e-14
    normError=normError, k=k
    error(['failed single var interpolation k=' num2str(k)])
  end
end
%{
\end{matlab}

\paragraph{Test multiple fields}
Use this to measure some of the errors in order to omit
singleton dimensions, and also squish any errors if the
second argument is essential zero (to cater for cosine
aliasing errors).
\begin{matlab}
%}
  normDiff=@(u,v) ...
  norm(squeeze(u)-squeeze(v))*norm(squeeze(v(i0,:,:,:)));
%{
\end{matlab}

Set a profile, and evaluate the interpolation. For the case
of the highest wavenumber, squash the error when the
centre-patch values are all zero by multiplying by result norm.
Not yet working for edgy interpolation.
\begin{matlab}
%}
for k=1:(nPatch-1)/2 % not checking the highest wavenumber
  u0=sin(k*patches.x*2*pi/Len);
  v0=cos(k*patches.x*2*pi/Len);
  uvi=patchEdgeInt1( reshape([u0 v0],[],1) );
  normuError=normDiff(uvi(:,1,:,:),u0);
  normvError=normDiff(uvi(:,2,:,:),v0);
  if abs(normuError)+abs(normvError)>2e-13
    normuError=normuError, normvError=normvError
    error(['failed double field interpolation k=' num2str(k)])
  end
end
%{
\end{matlab}

End the for-loop over various geometries.
\begin{matlab}
%}
end
disp('Passed standard spectral interpolation tests')
%{
\end{matlab}









\subsubsection{Now test spectral interpolation on staggered grid}
Must have even number of patches for a staggered grid.
\begin{matlab}
%}
disp('*** spectral interpolation on staggered grid')
for nPatch=repmat(6:2:20,1,nReal)
nPatch=nPatch
ratio=0.5*rand
nSubP=7; % of form 4*N-1
Len=10*rand
configPatches1(@simpleWavepde,[0,Len],nan,nPatch,-1,ratio,nSubP ...
    ,'EdgyInt',rand<0.5);
if mod(nPatch,4)==0, disp('Avoiding highest wavenumber'), end
kMax=floor((nPatch/2-1)/2) 
if patches.EdgyInt, i0=[2 nSubP-1], else i0=(nSubP+1)/2, end%??
patches.periodic=true; % temporary
%{
\end{matlab}
Identify which microscale grid points are \(h\)~or~\(u\) values.
\begin{matlab}
%}
uPts=mod( (1:nSubP)'+(1:nPatch) ,2);
hPts=find(1-uPts);
uPts=find(uPts);
%{
\end{matlab}
Set a profile for various wavenumbers. The capital
letter~\verb|U| denotes an array of values merged from
both~\(u\) and~\(h\) fields on the staggered grids.
\begin{matlab}
%}
fprintf('Staggered: single field-pair test.\n')
for k=-kMax:kMax
  U0=nan(nSubP,nPatch);
  U0(hPts)=rand*exp(+1i*k*patches.x(hPts)*2*pi/Len);
  U0(uPts)=rand*exp(-1i*k*patches.x(uPts)*2*pi/Len);
  Ui=patchEdgeInt1(U0(:));
  normError=norm(Ui(:)-U0(:));
  if abs(normError)>5e-14
    normError=normError
    patches=patches
    error(['staggered: failed single sys interpolation k=' num2str(k)])
  end
end
%{
\end{matlab}

\paragraph{Test multiple fields}
Use this to measure some of the errors in order to omit
singleton dimensions, and also squish any errors if the
third argument is essential zero (to cater for cosine
aliasing errors).
\begin{matlab}
%}
  normDiff=@(u,v,w) ...
  norm(squeeze(u)-squeeze(v))*norm(squeeze(w(i0,:,:,:)));
%{
\end{matlab}
Set a profile, and evaluate the interpolation. For the case
of the highest wavenumber zig-zag, squash the error when the
alternate centre-patch values are all zero. First shift the
\(x\)-coordinates so that the zig-zag mode is centred on a
patch.
\begin{matlab}
%}
fprintf('Staggered: Two field-pairs test.\n')
x0=patches.x((nSubP+1)/2,1);
patches.x=patches.x-x0;
oddP=1:2:nPatch; evnP=2:2:nPatch;
for k=1:kMax
  U0=nan(nSubP,1,1,nPatch); V0=U0;
  U0(hPts)=rand*sin(k*patches.x(hPts)*2*pi/Len);
  U0(uPts)=rand*sin(k*patches.x(uPts)*2*pi/Len);
  V0(hPts)=rand*cos(k*patches.x(hPts)*2*pi/Len);
  V0(uPts)=rand*cos(k*patches.x(uPts)*2*pi/Len);
  UVi=patchEdgeInt1([U0 V0]);
%  U0i0odd=squeeze(U0(i0,:,:,oddP)), U0i0evn=squeeze(U0(i0,:,:,evnP)) 
%  V0i0odd=squeeze(V0(i0,:,:,oddP)), V0i0evn=squeeze(V0(i0,:,:,evnP)) 
  normuError=[normDiff(UVi(:,1,:,oddP),U0(:,:,:,oddP),U0(:,:,:,evnP)) 
      normDiff(UVi(:,1,:,evnP),U0(:,:,:,evnP),U0(:,:,:,oddP))]';
  normvError=[normDiff(UVi(:,2,:,oddP),V0(:,:,:,oddP),V0(:,:,:,evnP)) 
      normDiff(UVi(:,2,:,evnP),V0(:,:,:,evnP),V0(:,:,:,oddP))]';
  if norm(normuError)+norm(normvError)>2e-13
    normuError=normuError, normvError=normvError
%    [U0 UVi(:,1,:,:) V0 UVi(:,2,:,:)]
    patches=patches
    error(['staggered: failed double field interpolation k=' num2str(k)])
  end
end
%{
\end{matlab}

End for-loop over patches
\begin{matlab}
%}
end
%{
\end{matlab}








\subsubsection{Check standard finite width interpolation}
Check over various types and orders of interpolation,
numbers of patches, random domain lengths and random ratios.
(The \verb|@sin| is a dummy.)
\begin{matlab}
%}
for edgyInt=[mod(nSubP,2)==0 true]
for ordCC=2:2:8
for nPatch=ordCC+(2:4)
    edgyInt=edgyInt
    ordCC=ordCC
    nPatch=nPatch
    Domain=5*[-rand rand]
    ratio=0.5*rand
    configPatches1(@sin,Domain,nan,nPatch,ordCC,ratio,nSubP ...
        ,'EdgyInt',edgyInt);
    patches.periodic=true; % temporary
%{
\end{matlab}

\paragraph{Check multiple fields simultaneously}
Set profiles to be various powers of~\(x\), \verb|ps|, and
store as different `variables' at each point.  
\begin{matlab}
%}
    ps=1:ordCC
    cs=randn(size(ps))
    u0=patches.x.^ps.*cs+randn;
%{
\end{matlab}
Then evaluate the interpolation and squeeze the singleton
dimension of an `ensemble'.
\begin{matlab}
%}
    ui=patchEdgeInt1(u0(:));
    ui=squeeze(ui);
%{
\end{matlab}
The interior patches should have zero error.
\begin{matlab}
%}
    j=ordCC/2+1:nPatch-ordCC/2;
    iError=ui(:,:,j)-u0(:,:,j);
    normError=norm(iError(:))
    assert(normError<5e-12 ...
    ,'failed finite stencil interpolation')
%{
\end{matlab}

End the for-loops over various parameters.
\begin{matlab}
%}
end,end,end
disp('Passed all standard polynomial interpolation')
%{
\end{matlab}
\endinput%%%%%%%%%%%%%%%%%%%%%%%%%%??












\subsubsection{Now test finite width interpolation on staggered grid}
Must have even number of patches for a staggered grid.
\begin{matlab}
%}
for nPatch=6:2:20
nPatch=nPatch
ratio=0.5*rand
nSubP=3; % of form 4*N-1
Len=10*rand
configPatches1(@simpleWavepde,[0,Len],nan,nPatch,-1,ratio,nSubP);
kMax=floor((nPatch/2-1)/2)
patches.periodic=true; % temporary
%{
\end{matlab}
Identify which microscale grid points are \(h\)~or~\(u\) values.
\begin{matlab}
%}
uPts=mod( (1:nSubP)'+(1:nPatch) ,2);
hPts=find(1-uPts);
uPts=find(uPts);
%{
\end{matlab}
Set a profile for various wavenumbers. The capital
letter~\verb|U| denotes an array of values merged from
both~\(u\) and~\(h\) fields on the staggered grids.
\begin{matlab}
%}
fprintf('Single field-pair test.\n')
for k=-kMax:kMax
  U0=nan(nSubP,nPatch);
  U0(hPts)=rand*exp(+1i*k*patches.x(hPts)*2*pi/Len);
  U0(uPts)=rand*exp(-1i*k*patches.x(uPts)*2*pi/Len);
  Ui=squeeze(patchEdgeInt1(U0(:)));
  normError=norm(Ui-U0);
  if abs(normError)>5e-14
    normError=normError
    error(['failed single sys interpolation k=' num2str(k)])
  end
end
%{
\end{matlab}

\paragraph{Test multiple fields}
Set a profile, and evaluate the interpolation. For the case
of the highest wavenumber zig-zag, squash the error when the
alternate centre-patch values are all zero. First shift the
\(x\)-coordinates so that the zig-zag mode is centred on a
patch.
\begin{matlab}
%}
i0=(nSubP+1)/2; % centre-patch index
fprintf('Two field-pairs test.\n')
x0=patches.x((nSubP+1)/2,1);
patches.x=patches.x-x0;
for k=1:nPatch/4
  U0=nan(nSubP,1,1,nPatch); V0=U0;
  U0(hPts)=rand*sin(k*patches.x(hPts)*2*pi/Len);
  U0(uPts)=rand*sin(k*patches.x(uPts)*2*pi/Len);
  V0(hPts)=rand*cos(k*patches.x(hPts)*2*pi/Len);
  V0(uPts)=rand*cos(k*patches.x(uPts)*2*pi/Len);
  UVi=patchEdgeInt1([U0 V0]);
  Ui=squeeze(UVi(:,1,1,:)); 
  Vi=squeeze(UVi(:,2,1,:));
  normuError=norm(Ui(:,1:2:nPatch)-U0(:,1:2:nPatch))*norm(U0(i0,2:2:nPatch)) ...
            +norm(Ui(:,2:2:nPatch)-U0(:,2:2:nPatch))*norm(U0(i0,1:2:nPatch));
  normvError=norm(Vi(:,1:2:nPatch)-V0(:,1:2:nPatch))*norm(V0(i0,2:2:nPatch)) ...
            +norm(Vi(:,2:2:nPatch)-V0(:,2:2:nPatch))*norm(V0(i0,1:2:nPatch));
  if abs(normuError)+abs(normvError)>2e-13
    normuError=normuError, normvError=normvError
    error(['failed double field interpolation k=' num2str(k)])
  end
end
%{
\end{matlab}

End for-loop over the cases
\begin{matlab}
%}
end
%{
\end{matlab}





\subsubsection{Finish}
If no error messages, then all OK.
\begin{matlab}
%}
fprintf('\nIf you read this, then all tests were passed\n')
%{
\end{matlab}
%}
