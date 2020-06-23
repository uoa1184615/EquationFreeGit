% Test the spectral interpolation of patchEdgeInt1()
% All tests passed.  But changing in 2020
% AJR, 26 Sep 2018 -- Jun 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{patchEdgeInt1test}: test the spectral interpolation}
\label{sec:patchEdgeInt1test}

A script to test the spectral interpolation of function
\verb|patchEdgeInt1()| Establish global data struct for the
range of various cases.
\begin{matlab}
%}
clear all
global patches
nSubP=3
i0=(nSubP+1)/2; % centre-patch index
%{
\end{matlab}

\paragraph{Test standard spectral interpolation}
Test over various numbers of patches, random domain lengths
and random ratios.
\begin{matlab}
%}
for nPatch=5:10
nPatch=nPatch
Len=10*rand
ratio=0.5*rand
configPatches1(@sin,[0,Len],nan,nPatch,0,ratio,nSubP ...
    ,'EdgyInt',rand<0.5); % random Edgy or not
kMax=floor((nPatch-1)/2);
if patches.EdgyInt, i0=[2 nSubP-1], else i0=(nSubP+1)/2, end
%{
\end{matlab}

\subparagraph{Test single field}
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

\subparagraph{Test multiple fields}
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
\begin{matlab}
%}
for k=1:nPatch/2
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




\paragraph{Now test spectral interpolation on staggered grid}
Must have even number of patches for a staggered grid.
Currently works for false EdgyInt?? 
\begin{matlab}
%}
disp('*** spectral interpolation on staggered grid')
for nPatch=6:2:20
nPatch=nPatch
ratio=0.5*rand
nSubP=7; % of form 4*N-1
Len=10*rand
configPatches1(@simpleWavepde,[0,Len],nan,nPatch,-1,ratio,nSubP ...
    ,'EdgyInt',true);
kMax=floor((nPatch/2-1)/2)
if patches.EdgyInt, i0=[2 nSubP-1], else i0=(nSubP+1)/2, end%??
%{
\end{matlab}
Identify which microscale grid points are \(h\)~or~\(u\) values.
\begin{matlab}
%}
uPts=mod( bsxfun(@plus,(1:nSubP)',(1:nPatch)) ,2);
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
disp('staggered: passed single sys interp')
%{
\end{matlab}

\subparagraph{Test multiple fields}
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
for k=1:nPatch/4
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




\paragraph{Finish}
If no error messages, then all OK.
\begin{matlab}
%}
fprintf('\nIf you read this, then all tests were passed\n')
%{
\end{matlab}
%}
