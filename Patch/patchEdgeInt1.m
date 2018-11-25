%Provides the interpolation across space for 1D patches of
%simulations of a smooth lattice system such as PDE
%discretisations.
%AJR, Sep 2018
%!TEX root = ../Doc/equationFreeDoc.tex
%{
\subsection{\texttt{patchEdgeInt1()}: sets edge values from interpolation over the macroscale}
\label{sec:patchEdgeInt1}
\localtableofcontents

Couples 1D patches across 1D space by computing their edge values
from macroscale interpolation. Consequently a spatially
discrete system could be integrated in time via the patch or
gap-tooth scheme \cite[]{Roberts06d}. Assumes that the
sub-patch structure is \emph{smooth} so that the patch
centre-values are sensible macroscale variables, and patch
edge values are determined by macroscale interpolation of
the patch-centre values. Communicate patch-design variables
via the global struct~\verb|patches|.
\begin{matlab}
%}
function u=patchEdgeInt1(u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|u| is a vector of length \(\verb|nSubP|\cdot
\verb|nPatch|\cdot \verb|nVars|\) where there are
\verb|nVars| field values at each of the points in the
\(\verb|nSubP|\times \verb|nPatch|\) grid.
\item \verb|patches| a struct set by \verb|configPatches1()|
which includes the following.
\begin{itemize}
\item \verb|.x| is \(\verb|nSubP|\times \verb|nPatch|\)
array of the spatial locations~\(x_{ij}\) of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.
\item \verb|.ordCC| is order of interpolation, currently in
\(\{0,2,4,6,8\}\).
\item \verb|.alt| in \(\{0,1\}\) is one for staggered grid
(alternating) interpolation.
\item \verb|.Cwtsr| and \verb|.Cwtsl|
\end{itemize}
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|u| is \(\verb|nSubP|\times \verb|nPatch|\times
\verb|nVars|\) 2/3D array of the fields with edge values set by
interpolation.
\end{itemize}


\begin{body}
Determine the sizes of things. Any error arising in the
reshape indicates~\verb|u| has the wrong size.
\begin{matlab}
%}
[nSubP,nPatch]=size(patches.x);
nVars=round(numel(u)/numel(patches.x));
if numel(u)~=nSubP*nPatch*nVars
  nSubP=nSubP, nPatch=nPatch, nVars=nVars, sizeu=size(u)
  end
u=reshape(u,nSubP,nPatch,nVars);
%{
\end{matlab}
With Dirichlet patches, the half-length of a patch is
\(h=dx(n_\mu-1)/2\) (or~\(-2\) for specified flux), and the
ratio needed for interpolation is then \(r=h/\Delta X\).
Compute lattice sizes from inside the patches as the edge
values may be \nan{}s, etc.
\begin{matlab}
%}
dx=patches.x(3,1)-patches.x(2,1);
DX=patches.x(2,2)-patches.x(2,1);
r=dx*(nSubP-1)/2/DX;
%{
\end{matlab}

For the moment assume the physical domain is macroscale
periodic so that the coupling formulas are simplest. Should
eventually cater for periodic, odd-mid-gap, even-mid-gap,
even-mid-patch, dirichlet, neumann, ?? These index vectors
point to patches and their two immediate neighbours.
\begin{matlab}
%}
j=1:nPatch; jp=mod(j,nPatch)+1; jm=mod(j-2,nPatch)+1;
%{
\end{matlab}
The centre of each patch (as \verb|nSubP| is odd) is at
\begin{matlab}
%}
i0=round((nSubP+1)/2);
%{
\end{matlab}

\paragraph{Lagrange interpolation gives patch-edge values}
So compute centred differences of the mid-patch values for
the macro-interpolation, of all fields. Assumes the domain
is macro-periodic.
\begin{matlab}
%}
if patches.ordCC>0 % then non-spectral interpolation
  dmu=nan(patches.ordCC,nPatch,nVars);
  if patches.alt % use only odd numbered neighbours
    dmu(1,:,:)=(u(i0,jp,:)+u(i0,jm,:))/2; % \mu
    dmu(2,:,:)= u(i0,jp,:)-u(i0,jm,:); % \delta
    jp=jp(jp); jm=jm(jm); % increase shifts to \pm2
  else % standard
    dmu(1,:,:)=(u(i0,jp,:)-u(i0,jm,:))/2; % \mu\delta
    dmu(2,:,:)=(u(i0,jp,:)-2*u(i0,j,:)+u(i0,jm,:)); % \delta^2
  end% if odd/even
%{
\end{matlab}
Recursively take \(\delta^2\) of these to form higher order
centred differences (could unroll a little to cater for two
in parallel).
\begin{matlab}
%}
  for k=3:patches.ordCC
    dmu(k,:,:)=dmu(k-2,jp,:)-2*dmu(k-2,j,:)+dmu(k-2,jm,:);
  end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet edge values for
each patch \cite[]{Roberts06d}, using weights computed in
\verb|configPatches1()| . Here interpolate to specified order.
\begin{matlab}
%}
  u(nSubP,j,:)=u(i0,j,:)*(1-patches.alt) ...
    +sum(bsxfun(@times,patches.Cwtsr,dmu));
  u( 1,j,:)=u(i0,j,:)*(1-patches.alt) ...
    +sum(bsxfun(@times,patches.Cwtsl,dmu));
%{
\end{matlab}

\paragraph{Case of spectral interpolation}
Assumes the domain is macro-periodic. As the macroscale
fields are \(N\)-periodic, the macroscale Fourier transform
writes the centre-patch values as \(U_j=\sum_{k}C_ke^{ik2\pi
j/N}\). Then the edge-patch values \(U_{j\pm r}
=\sum_{k}C_ke^{ik2\pi/N(j\pm r)} =\sum_{k}C'_ke^{ik2\pi
j/N}\) where \(C'_k=C_ke^{ikr2\pi/N}\). For
\verb|nPatch|~patches we resolve `wavenumbers'
\(|k|<\verb|nPatch|/2\), so set row vector
\(\verb|ks|=k2\pi/N\) for `wavenumbers'
\(k=(0,1, \ldots, k_{\max}, -k_{\max}, \ldots, -1)\) for odd~\(N\), and \(k=(0,1, \ldots, k_{\max}, \pm(k_{\max}+1), -k_{\max}, \ldots, -1)\) for even~\(N\).
\begin{matlab}
%}
else% spectral interpolation
%{
\end{matlab}
Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches1| tests
there are an even number of patches). Then the patch-ratio
is effectively halved. The patch edges are near the middle
of the gaps and swapped.
\begin{matlab}
%}
  if patches.alt % transform by doubling the number of fields
    v=nan(size(u)); % currently to restore the shape of u
    u=cat(3,u(:,1:2:nPatch,:),u(:,2:2:nPatch,:));
    altShift=reshape(0.5*[ones(nVars,1);-ones(nVars,1)],1,1,[]);
    iV=[nVars+1:2*nVars 1:nVars]; % scatter interp to alternate field
    r=r/2;           % ratio effectively halved
    nPatch=nPatch/2; % halve the number of patches
    nVars=nVars*2;   % double the number of fields
  else % the values for standard spectral
    altShift=0;  
    iV=1:nVars;  
  end
%{
\end{matlab}
Now set wavenumbers.
\begin{matlab}
%}
  kMax=floor((nPatch-1)/2); 
  ks=2*pi/nPatch*(mod((0:nPatch-1)+kMax,nPatch)-kMax); 
%{
\end{matlab}
Test for reality of the field values, and define a function
accordingly.
\begin{matlab}
%}
  if imag(u(i0,:,:))==0, uclean=@(u) real(u);
    else                 uclean=@(u) u; 
    end
%{
\end{matlab}
Compute the Fourier transform of the patch centre-values for
all the fields. If there are an even number of points, then
zero the zig-zag mode in the \textsc{ft} and add it in later
as cosine.
\begin{matlab}
%}
  Ck=fft(u(i0,:,:));
  if mod(nPatch,2)==0
    Czz=Ck(1,nPatch/2+1,:)/nPatch; 
    Ck(1,nPatch/2+1,:)=0; 
  end 
%{
\end{matlab}
The inverse Fourier transform gives the edge values via a
shift a fraction~\(r\) to the next macroscale grid point.
Enforce reality when appropriate. 
\begin{matlab}
%}
  u(nSubP,:,iV)=uclean(ifft(bsxfun(@times,Ck ...
      ,exp(1i*bsxfun(@times,ks,altShift+r)))));
  u( 1,:,iV)=uclean(ifft(bsxfun(@times,Ck ...
      ,exp(1i*bsxfun(@times,ks,altShift-r)))));
%{
\end{matlab}
For an even number of patches, add in the cosine mode.
\begin{matlab}
%}
  if mod(nPatch,2)==0
    cosr=cos(pi*(altShift+r+(0:nPatch-1)));
    u(nSubP,:,iV)=u(nSubP,:,iV)+uclean(bsxfun(@times,Czz,cosr));
    cosr=cos(pi*(altShift-r+(0:nPatch-1)));
    u( 1,:,iV)=u( 1,:,iV)+uclean(bsxfun(@times,Czz,cosr));
  end
%{
\end{matlab}
Restore staggered grid when appropriate. Is there a better
way to do this??
\begin{matlab}
%}
if patches.alt
  nVars=nVars/2;  nPatch=2*nPatch;
  v(:,1:2:nPatch,:)=u(:,:,1:nVars);
  v(:,2:2:nPatch,:)=u(:,:,nVars+1:2*nVars);    
  u=v;
end
end% if spectral 
%{
\end{matlab}
Fin, returning the 2/3D array of field values.  
\end{body}
%}
