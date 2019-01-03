% Provides the interpolation across 2D space for 2D patches
% of simulations of a smooth lattice system such as PDE
% discretisations.  AJR, Nov 2018
%!TEX root = ../Doc/equationFreeDoc.tex
%{
\subsection[\texttt{patchEdgeInt2()}: 2D patch edge values from 2D interpolation]{\texttt{patchEdgeInt2()}: sets 2D patch edge values from 2D macroscale interpolation}
\label{sec:patchEdgeInt2}
\localtableofcontents


Couples 2D patches across 2D space by computing their edge
values via macroscale interpolation.  Assumes that the
sub-patch structure is \emph{smooth} so that the patch
centre-values are sensible macroscale variables, and patch
edge values are determined by macroscale interpolation of
the patch-centre values.  Communicate patch-design variables
via the global struct~\verb|patches|.
\begin{matlab}
%}
function u=patchEdgeInt2(u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|u| is a vector of length \(\verb|nx| \cdot
\verb|ny| \cdot \verb|Nx| \cdot \verb|Ny| \cdot
\verb|nVars|\) where there are \verb|nVars| field values at
each of the points in the \(\verb|nx| \times \verb|ny|\times
\verb|Nx| \times \verb|Ny|\) grid on the \(\verb|Nx| \times
\verb|Ny|\) array of patches.
\item \verb|patches| a struct set by \verb|configPatches2()|
which includes the following information.
\begin{itemize}
\item \verb|.x| is \(\verb|nx|\times \verb|Nx|\)
array of the spatial locations~\(x_{ij}\) of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.
\item \verb|.y| is similarly \(\verb|ny|\times \verb|Ny|\)
array of the spatial locations~\(y_{ij}\) of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.
\item \verb|.ordCC| is order of interpolation, currently only
\(\{0\}\).
%\item \verb|.alt| in \(\{0,1\}\) is one for staggered grid
%(alternating) interpolation.
\item \verb|.Cwtsr| and \verb|.Cwtsl|---not yet used
\end{itemize}
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item \verb|u| is \(\verb|nx| \times \verb|ny|\times
\verb|Nx| \times \verb|Ny| \times \verb|nVars|\) array of
the fields with edge values set by interpolation.
\end{itemize}


\begin{body}

Determine the sizes of things. Any error arising in the
reshape indicates~\verb|u| has the wrong size.
\begin{matlab}
%}
[ny,Ny] = size(patches.y);
[nx,Nx] = size(patches.x);
nVars = round(numel(u)/numel(patches.x)/numel(patches.y));
if numel(u) ~= nx*ny*Nx*Ny*nVars
  nSubP=[nx ny], nPatch=[Nx Ny], nVars=nVars, sizeu=size(u)
end
u = reshape(u,[nx ny Nx Ny nVars]);
%{
\end{matlab}
With Dirichlet patches, the half-length of a patch is
\(h=dx(n_\mu-1)/2\) (or~\(-2\) for specified flux), and the
ratio needed for interpolation is then \(r=h/\Delta X\).
Compute lattice sizes from inside the patches as the edge
values may be \nan{}s, etc.
\begin{matlab}
%}
dx = patches.x(3,1)-patches.x(2,1);
DX = patches.x(2,2)-patches.x(2,1);
rx = dx*(nx-1)/2/DX;
dy = patches.y(3,1)-patches.y(2,1);
DY = patches.y(2,2)-patches.y(2,1);
ry = dy*(ny-1)/2/DY;
%{
\end{matlab}

For the moment assume the physical domain is macroscale
periodic so that the coupling formulas are simplest. Should
eventually cater for periodic, odd-mid-gap, even-mid-gap,
even-mid-patch, dirichlet, neumann, ?? These index vectors
point to patches and their two immediate neighbours.
\begin{matlab}
%}
%i=1:Nx; ip=mod(i,Nx)+1; im=mod(j-2,Nx)+1;
%j=1:Ny; jp=mod(j,Ny)+1; jm=mod(j-2,Ny)+1;
%{
\end{matlab}
The centre of each patch (as \verb|nx| and~\verb|ny| are
odd) is at
\begin{matlab}
%}
i0 = round((nx+1)/2);
j0 = round((ny+1)/2);
%{
\end{matlab}

\paragraph{Lagrange interpolation gives patch-edge values}
So compute centred differences of the mid-patch values for
the macro-interpolation, of all fields. Assumes the domain
is macro-periodic.
\begin{matlab}
%}
if patches.ordCC>0 % then non-spectral interpolation
error('non-spectral interpolation not yet implemented')
  dmu=nan(patches.ordCC,nPatch,nVars);
%  if patches.alt % use only odd numbered neighbours
%    dmu(1,:,:)=(u(i0,jp,:)+u(i0,jm,:))/2; % \mu
%    dmu(2,:,:)= u(i0,jp,:)-u(i0,jm,:); % \delta
%    jp=jp(jp); jm=jm(jm); % increase shifts to \pm2
%  else % standard
    dmu(1,:,:)=(u(i0,jp,:)-u(i0,jm,:))/2; % \mu\delta
    dmu(2,:,:)=(u(i0,jp,:)-2*u(i0,j,:)+u(i0,jm,:)); % \delta^2
%  end% if odd/even
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
\verb|configPatches2()| . Here interpolate to specified order.
\begin{matlab}
%}
  u(nSubP,j,:)=u(i0,j,:)*(1-patches.alt) ...
    +sum(bsxfun(@times,patches.Cwtsr,dmu));
  u( 1,j,:)=u(i0,j,:)*(1-patches.alt) ...
    +sum(bsxfun(@times,patches.Cwtsl,dmu));
%{
\end{matlab}

\paragraph{Case of spectral interpolation}
Assumes the domain is macro-periodic. We interpolate in
terms of the patch index~\(j\), say, not directly in space. 
As the macroscale fields are \(N\)-periodic in the patch
index~\(j\), the macroscale Fourier transform writes the
centre-patch values as \(U_j=\sum_{k}C_ke^{ik2\pi j/N}\).
Then the edge-patch values \(U_{j\pm r}
=\sum_{k}C_ke^{ik2\pi/N(j\pm r)} =\sum_{k}C'_ke^{ik2\pi
j/N}\) where \(C'_k=C_ke^{ikr2\pi/N}\). For \(N\)~patches we
resolve `wavenumbers' \(|k|<N/2\), so set row vector
\(\verb|ks|=k2\pi/N\) for `wavenumbers' \(\mathcode`\,="213B
k=(0,1, \ldots, k_{\max}, -k_{\max}, \ldots, -1)\) for
odd~\(N\), and \(\mathcode`\,="213B k=(0,1, \ldots,
k_{\max}, \pm(k_{\max}+1) -k_{\max}, \ldots, -1)\) for
even~\(N\).
\begin{matlab}
%}
else% spectral interpolation
%{
\end{matlab}
Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches2|
tests there are an even number of patches). Then the
patch-ratio is effectively halved. The patch edges are near
the middle of the gaps and swapped.
\begin{matlab}
%}
%  if patches.alt % transform by doubling the number of fields
%  error('staggered grid not yet implemented')
%    v=nan(size(u)); % currently to restore the shape of u
%    u=cat(3,u(:,1:2:nPatch,:),u(:,2:2:nPatch,:));
%    altShift=reshape(0.5*[ones(nVars,1);-ones(nVars,1)],1,1,[]);
%    iV=[nVars+1:2*nVars 1:nVars]; % scatter interp to alternate field
%    r=r/2;           % ratio effectively halved
%    nPatch=nPatch/2; % halve the number of patches
%    nVars=nVars*2;   % double the number of fields
%  else % the values for standard spectral
    altShift = 0;  
    iV = 1:nVars;
%  end
%{
\end{matlab}
Now set wavenumbers in the two directions.  In the case of
even~\(N\) these compute the \(+\)-case for the highest
wavenumber zig-zag mode, \(\mathcode`\,="213B k=(0,1,
\ldots, k_{\max}, +(k_{\max}+1) -k_{\max}, \ldots, -1)\).
\begin{matlab}
%}
  kMax = floor((Nx-1)/2); 
  krx = rx*2*pi/Nx*(mod((0:Nx-1)+kMax,Nx)-kMax); 
  kMay = floor((Ny-1)/2); 
  kry = ry*2*pi/Ny*(mod((0:Ny-1)+kMay,Ny)-kMay); 
%{
\end{matlab}
Test for reality of the field values, and define a function
accordingly.
\begin{matlab}
%}
  if imag(u(i0,j0,:,:,:))==0, uclean = @(u) real(u);
     else uclean = @(u) u; end
%{
\end{matlab}
Compute the Fourier transform of the patch centre-values for
all the fields. If there are an even number of points, then
zero the zig-zag mode in the \textsc{ft} and add it in later
as cosine.
\begin{matlab}
%}
  Ck = fft2(squeeze(u(i0,j0,:,:,:)));
%{
\end{matlab}
The inverse Fourier transform gives the edge values via a
shift a fraction~\(\verb|rx|/\verb|ry|\) to the next
macroscale grid point. Initially preallocate storage for all
the \textsc{ifft}s that we need to cater for the zig-zag
modes when there are an even number of patches in the
directions.
\begin{matlab}
%}
nFTx = 2-mod(Nx,2);
nFTy = 2-mod(Ny,2);
unj = nan(1,ny,Nx,Ny,nVars,nFTx*nFTy);
u1j = nan(1,ny,Nx,Ny,nVars,nFTx*nFTy);
uin = nan(nx,1,Nx,Ny,nVars,nFTx*nFTy);
ui1 = nan(nx,1,Nx,Ny,nVars,nFTx*nFTy);
%{
\end{matlab}
Loop over the required \textsc{ifft}s.
\begin{matlab}
%}
iFT = 0;
for iFTx = 1:nFTx
for iFTy = 1:nFTy
iFT = iFT+1;
%{
\end{matlab}
First interpolate onto \(x\)-limits of the patches. (It may
be more efficient to product exponentials of vectors,
instead of exponential of array---only for \(N>100\).  Can
this be vectorised further??)
\begin{matlab}
%}
for jj = 1:ny 
  ks = (jj-j0)*2/(ny-1)*kry; % fraction of kry along the edge
  unj(1,jj,:,:,iV,iFT) = ifft2( bsxfun(@times,Ck ...
      ,exp(1i*bsxfun(@plus,altShift+krx',ks))));
  u1j(1,jj,:,:,iV,iFT) = ifft2( bsxfun(@times,Ck ...
      ,exp(1i*bsxfun(@plus,altShift-krx',ks))));
end
%{
\end{matlab}
Second interpolate onto \(y\)-limits of the patches.
\begin{matlab}
%}
for i = 1:nx 
  ks = (i-i0)*2/(nx-1)*krx; % fraction of krx along the edge
  uin(i,1,:,:,iV,iFT) = ifft2( bsxfun(@times,Ck ...
      ,exp(1i*bsxfun(@plus,ks',altShift+kry))));
  ui1(i,1,:,:,iV,iFT) = ifft2( bsxfun(@times,Ck ...
      ,exp(1i*bsxfun(@plus,ks',altShift-kry))));
end
%{
\end{matlab}
When either direction have even number of patches then swap
the zig-zag wavenumber to the conjugate.
\begin{matlab}
%}
if nFTy==2, kry(Ny/2+1) = -kry(Ny/2+1); end
end% iFTy-loop
if nFTx==2, krx(Nx/2+1) = -krx(Nx/2+1); end
end% iFTx-loop
%{
\end{matlab}
Put edge-values into the \(u\)-array, using \verb|mean()| to
treat a zig-zag mode as cosine. Enforce reality when
appropriate via \verb|uclean()|. 
\begin{matlab}
%}
u(end,:,:,:,iV) = uclean( mean(unj,6) );
u( 1 ,:,:,:,iV) = uclean( mean(u1j,6) );
u(:,end,:,:,iV) = uclean( mean(uin,6) );
u(:, 1 ,:,:,iV) = uclean( mean(ui1,6) );
%{
\end{matlab}
Restore staggered grid when appropriate. Is there a better
way to do this??
\begin{matlab}
%}
%if patches.alt
%  nVars=nVars/2;  nPatch=2*nPatch;
%  v(:,1:2:nPatch,:)=u(:,:,1:nVars);
%  v(:,2:2:nPatch,:)=u(:,:,nVars+1:2*nVars);    
%  u=v;
%end
end% if spectral 
end% function patchEdgeInt2
%{
\end{matlab}
Fin, returning the 4/5D array of field values with
interpolated edges. 
\end{body} 
%}
