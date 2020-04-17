% Provides the interpolation across 2D space for 2D patches
% of simulations of a smooth lattice system such as PDE
% discretisations.  AJR, Nov 2018 -- 17 Apr 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section[\texttt{patchEdgeInt2()}: 2D patch edge values from 2D interpolation]
{\texttt{patchEdgeInt2()}: sets 2D patch edge values from 2D macroscale interpolation}
\label{sec:patchEdgeInt2}
%\localtableofcontents


Couples 2D patches across 2D space by computing their edge
values via macroscale interpolation.  Assumes that the
sub-patch structure is \emph{smooth} so that the patch
centre-values are sensible macroscale variables, and patch
edge values are determined by macroscale interpolation of
the patch-centre values.  Communicate patch-design variables
via the global struct~\verb|patches|.
\begin{matlab}
%}
function u = patchEdgeInt2(u)
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
equi-spaced lattice on both macro- and microscales.

\item \verb|.y| is similarly \(\verb|ny|\times \verb|Ny|\)
array of the spatial locations~\(y_{ij}\) of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and microscales.

\item \verb|.ordCC| is order of interpolation, currently only
\(\{0\}\).
%\item \verb|.alt| in \(\{0,1\}\) is one for staggered grid
%(alternating) interpolation.

\item \verb|.Cwtsr| and \verb|.Cwtsl|---not yet used

\item \verb|.EdgyInt| in \(\{0,1\}\) is one for
interpolating patch-edge values from opposite next-to-edge
values (often preserves symmetry).
\end{itemize}
\end{itemize}

\paragraph{Output}
\begin{itemize}
\item \verb|u| is \(\verb|nx| \times \verb|ny|\times
\verb|Nx| \times \verb|Ny| \times \verb|nVars|\) array of
the fields with edge values set by interpolation.
\end{itemize}





\begin{devMan}

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
\(h=dx(n_\mu-1)/2\) (or~\(-2\) for specified flux), unless
we are interpolating from next-to-edge values, and the
ratio needed for interpolation is then \(r=h/\Delta X\).
Compute lattice sizes from inside the patches as the edge
values may be \nan{}s, etc.
\begin{matlab}
%}
dx = patches.x(3,1)-patches.x(2,1);
DX = patches.x(2,2)-patches.x(2,1);
if patches.EdgyInt==0, rx = dx*(nx-1)/2/DX;
    else               rx = dx*(nx-2)/DX;
    end
dy = patches.y(3,1)-patches.y(2,1);
DY = patches.y(2,2)-patches.y(2,1);
if patches.EdgyInt==0, ry = dy*(ny-1)/2/DY;
    else               ry = dy*(ny-2)/DY;
    end

%{
\end{matlab}

For the moment assume the physical domain is macroscale
periodic so that the coupling formulas are simplest. Should
eventually cater for periodic, odd-mid-gap, even-mid-gap,
even-mid-patch, Dirichlet, Neumann, Robin?? These index vectors
point to patches and their two immediate neighbours---currently 
not needed.
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

\paragraph{Lagrange interpolation gives patch-edge values
---not yet implemented}
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
\verb|configPatches2()|. Here interpolate to specified order.
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
Now set wavenumbers in the two directions into two row
vectors.  In the case of even~\(N\) these compute the
\(+\)-case for the highest wavenumber zig-zag mode,
\(\mathcode`\,="213B k=(0,1, \ldots, k_{\max}, +(k_{\max}+1)
-k_{\max}, \ldots, -1)\).
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
  if max(abs(imag(u(:))))<1e-9*max(abs(u(:)))
       uclean=@(u) real(u);
  else uclean=@(u) u; 
  end
%{
\end{matlab}
Compute the Fourier transform of the mid-patch values for
all the fields.  Unless doing patch-edgy interpolation when
FT the next-to-edge values.  If there are an even number of
patches, then treat highest wavenumber as positive, but then
take only real-part when real data.
\begin{matlab}
%}
ix=2:nx-1;  iy=2:ny-1; % indices of interior
if patches.EdgyInt==0
     Cle = fft2(shiftdim(u(i0,j0,:,:,:),2)); 
     Cbo=Cle;
else Cle = fft2(shiftdim(u(   2,iy ,:,:,:),2));
     Cri = fft2(shiftdim(u(nx-1,iy ,:,:,:),2));
     Cbo = fft2(shiftdim(u(ix,2    ,:,:,:),2));
     Cto = fft2(shiftdim(u(ix,ny-1 ,:,:,:),2));
end
% fill in the cross of Fourier-shifted mid-values
if patches.EdgyInt==0 
  % y-fraction of kry along left/right edges
  ks = (shiftdim(iy ,-3)-j0)*2/(ny-1).*kry; 
  Cle = bsxfun(@times,Cle,exp(1i*ks)); 
  Cri = Cle;
  % x-fraction of krx along bottom/top edges
  ks = (shiftdim(ix',-3)-i0)*2/(nx-1).*krx'; 
  Cbo = bsxfun(@times,Cbo,exp(1i*ks)); 
  Cto = Cbo;
end
% put edge values into the patch field
u(end,iy,:,:,:) = uclean( shiftdim( ifft2( ...
    bsxfun(@times,Cle,exp(1i*(altShift+krx')))  ),2+(nVars>1)) );
u( 1 ,iy,:,:,:) = uclean( shiftdim( ifft2( ...
    bsxfun(@times,Cri,exp(1i*(altShift-krx')))  ),2+(nVars>1)) );
u(ix,end,:,:,:) = uclean( shiftdim( ifft2( ...
    bsxfun(@times,Cbo,exp(1i*(altShift+kry)))  ),2+(nVars>1)) );
u(ix, 1 ,:,:,:) = uclean( shiftdim( ifft2( ...
    bsxfun(@times,Cto,exp(1i*(altShift-kry)))  ),2+(nVars>1)) );

end% if spectral 
end% function patchEdgeInt2
%{
\end{matlab}
Fin, returning the 4/5D array of field values with
interpolated edges. 
\end{devMan} 
%}
