% patchEdgeInt3() provides the interpolation across 3D space
% for 3D patches of simulations of a smooth lattice system
% such as PDE discretisations.  AJR, Aug--Dec 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt3()}: sets 3D patch
face values from 3D macroscale interpolation}
\label{sec:patchEdgeInt3}

Couples 3D patches across 3D space by computing their face
values via macroscale interpolation.  Assumes that the patch
centre-values are sensible macroscale variables, and patch
face values are determined by macroscale interpolation of
the patch centre-plane values \cite[]{Roberts2011a,
Bunder2019c}, or patch next-to-face values which appears
better \cite[]{Bunder2020a}. This function is primarily used
by \verb|patchSys3()| but is also useful for user
graphics. 
\footnote{Script \texttt{patchEdgeInt3test.m} verifies this code.}

Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|), or otherwise via the global
struct \verb|patches|.
\begin{matlab}
%}
function u = patchEdgeInt3(u,patches)
if nargin<2, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector\slash array of length
$\verb|prod(nSubP)|  \cdot \verb|nVars| \cdot \verb|nEnsem|
\cdot \verb|prod(nPatch)|$ where there are $\verb|nVars|
\cdot \verb|nEnsem|$ field values at each of the points in
the $\verb|nSubP1| \cdot \verb|nSubP2| \cdot \verb|nSubP3|
\cdot \verb|nPatch1| \cdot \verb|nPatch2| \cdot
\verb|nPatch3|$ multiscale spatial grid on the
$\verb|nPatch1| \cdot \verb|nPatch2| \cdot \verb|nPatch3|$
array of patches.

\item \verb|patches| a struct set by \verb|configPatches3()|
which includes the following information.
\begin{itemize}

\item \verb|.x| is $\verb|nSubP1| \times1 \times1 \times1
\times1 \times \verb|nPatch1| \times1 \times1 $ array of
the spatial locations~$x_{iI}$ of the microscale grid
points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.

\item \verb|.y| is similarly $1\times \verb|nSubP2| \times1
\times1 \times1 \times1 \times \verb|nPatch2| \times1$
array of the spatial locations~$y_{jJ}$ of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.

\item \verb|.z| is similarly $1 \times1 \times \verb|nSubP3|
\times1 \times1 \times1 \times1 \times \verb|nPatch3|$
array of the spatial locations~$z_{kK}$ of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.

\item \verb|.ordCC| is order of interpolation, currently
only $\{0,2,4,\ldots\}$

\item \verb|.stag| in $\{0,1\}$ is one for staggered grid
(alternating) interpolation.  Currently must be zero.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the coupling
coefficients for finite width interpolation in each of the
$x,y,z$-directions.

\item \verb|.EdgyInt| true/false is true for interpolating
patch-face values from opposite next-to-face values (often
preserves symmetry).

\item \verb|.nEnsem| the number of realisations in the ensemble.

\item \verb|.parallel| whether serial or parallel.

\end{itemize}
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 8D array, $\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nSubP3| \cdot \verb|nVars| \cdot
\verb|nEnsem| \cdot \verb|nPatch1| \cdot \verb|nPatch2|
\cdot \verb|nPatch3|$, of the fields with face values set
by interpolation (edge and corner vales set to~\verb|NaN|).
\end{itemize}







\begin{devMan}

Test for reality of the field values, and define a function
accordingly.  Could be problematic if some variables are
real and some are complex, or if variables are of quite
different sizes. 
\begin{matlab}
%}
  if max(abs(imag(u(:))))<1e-9*max(abs(u(:)))
       uclean=@(u) real(u);
  else uclean=@(u) u; 
  end
%{
\end{matlab}

Determine the sizes of things. Any error arising in the
reshape indicates~\verb|u| has the wrong size.
\begin{matlab}
%}
[~,~,nz,~,~,~,~,Nz] = size(patches.z);
[~,ny,~,~,~,~,Ny,~] = size(patches.y);
[nx,~,~,~,~,Nx,~,~] = size(patches.x);
nEnsem = patches.nEnsem;
nVars = round( numel(u)/numel(patches.x) ...
    /numel(patches.y)/numel(patches.z)/nEnsem );
assert(numel(u) == nx*ny*nz*Nx*Ny*Nz*nVars*nEnsem ...
  ,'patchEdgeInt3: input u has wrong size for parameters')
u = reshape(u,[nx ny nz nVars nEnsem Nx Ny Nz]);
%{
\end{matlab}
Get the size ratios of the patches in each direction.
\begin{matlab}
%}
rx = patches.ratio(1);
ry = patches.ratio(2);
rz = patches.ratio(3);
%{
\end{matlab}

For the moment assume the physical domain is macroscale
periodic so that the coupling formulas are simplest. Should
eventually cater for periodic, odd-mid-gap, even-mid-gap,
even-mid-patch, Dirichlet, Neumann, Robin?? These index
vectors point to patches and their six immediate neighbours.
\begin{matlab}
%}
I=1:Nx; Ip=mod(I,Nx)+1; Im=mod(I-2,Nx)+1;
J=1:Ny; Jp=mod(J,Ny)+1; Jm=mod(J-2,Ny)+1;
K=1:Nz; Kp=mod(K,Nz)+1; Km=mod(K-2,Nz)+1;
%{
\end{matlab}
The centre of each patch (as \verb|nx|, \verb|ny|
and~\verb|nz| are odd for centre-patch interpolation) is at
indices
\begin{matlab}
%}
i0 = round((nx+1)/2);
j0 = round((ny+1)/2);
k0 = round((nz+1)/2);
%{
\end{matlab}



\paragraph{Lagrange interpolation gives patch-face values}
Compute centred differences of the mid-patch values for the
macro-interpolation, of all fields. Assumes the domain is
macro-periodic.
\begin{matlab}
%}
ordCC=patches.ordCC;
if ordCC>0 % then finite-width polynomial interpolation
%{
\end{matlab}
The patch-edge values are either interpolated from the
next-to-edge-face values, or from the centre-cross-plane
values (not the patch-centre value itself as that seems to
have worse properties in general).  Have not yet implemented
core averages??
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    Ux = u([2 nx-1],2:(ny-1),2:(nz-1),:,:,I,J,K);
    Uy = u(2:(nx-1),[2 ny-1],2:(nz-1),:,:,I,J,K);
    Uz = u(2:(nx-1),2:(ny-1),[2 nz-1],:,:,I,J,K);
  else % interpolate centre-cross values
    Ux = u(i0,2:(ny-1),2:(nz-1),:,:,I,J,K);
    Uy = u(2:(nx-1),j0,2:(nz-1),:,:,I,J,K);
    Uz = u(2:(nx-1),2:(ny-1),k0,:,:,I,J,K);
  end;
%{
\end{matlab}
Just in case the last array dimension(s) are one, we have to
force a padding of the sizes, then adjoin the extra
dimension for the subsequent array of differences.
\begin{matlab}
%}
szUxO=size(Ux); szUxO=[szUxO ones(1,8-length(szUxO)) ordCC];
szUyO=size(Uy); szUyO=[szUyO ones(1,8-length(szUyO)) ordCC];
szUzO=size(Uz); szUzO=[szUzO ones(1,8-length(szUzO)) ordCC];
%{
\end{matlab}
Use finite difference formulas for the interpolation, so
store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in these arrays.  When parallel, in order
to preserve the distributed array structure we use an index
at the end for the differences.
\begin{matlab}
%}
  if patches.parallel
    dmux = zeros(szUxO,patches.codist); % 9D
    dmuy = zeros(szUyO,patches.codist); % 9D
    dmuz = zeros(szUzO,patches.codist); % 9D
  else
    dmux = zeros(szUxO); % 9D
    dmuy = zeros(szUyO); % 9D
    dmuz = zeros(szUzO); % 9D
  end
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$ in
both space directions.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
    dmux(:,:,:,:,:,I,:,:,1) = (Ux(:,:,:,:,:,Ip,:,:) +Ux(:,:,:,:,:,Im,:,:))/2; % \mu
    dmux(:,:,:,:,:,I,:,:,2) = (Ux(:,:,:,:,:,Ip,:,:) -Ux(:,:,:,:,:,Im,:,:)); % \delta
    Ip = Ip(Ip); Im = Im(Im); % increase shifts to \pm2
    dmuy(:,:,:,:,:,:,J,:,1) = (Ux(:,:,:,:,:,:,Jp,:)+Ux(:,:,:,:,:,:,Jm,:))/2; % \mu
    dmuy(:,:,:,:,:,:,J,:,2) = (Ux(:,:,:,:,:,:,Jp,:)-Ux(:,:,:,:,:,:,Jm,:)); % \delta
    Jp = Jp(Jp); Jm = Jm(Jm); % increase shifts to \pm2
    dmuz(:,:,:,:,:,:,:,K,1) = (Ux(:,:,:,:,:,:,:,Kp)+Ux(:,:,:,:,:,:,:,Km))/2; % \mu
    dmuz(:,:,:,:,:,:,:,K,2) = (Ux(:,:,:,:,:,:,:,Kp)-Ux(:,:,:,:,:,:,:,Km)); % \delta
    Kp = Kp(Kp); Km = Km(Km); % increase shifts to \pm2
  else    
    dmux(:,:,:,:,:,I,:,:,1) = (Ux(:,:,:,:,:,Ip,:,:) ...
                              -Ux(:,:,:,:,:,Im,:,:))/2; %\mu\delta 
    dmux(:,:,:,:,:,I,:,:,2) = (Ux(:,:,:,:,:,Ip,:,:) ...
       -2*Ux(:,:,:,:,:,I,:,:) +Ux(:,:,:,:,:,Im,:,:));   %\delta^2    
    dmuy(:,:,:,:,:,:,J,:,1) = (Uy(:,:,:,:,:,:,Jp,:) ...
                              -Uy(:,:,:,:,:,:,Jm,:))/2; %\mu\delta 
    dmuy(:,:,:,:,:,:,J,:,2) = (Uy(:,:,:,:,:,:,Jp,:) ...
       -2*Uy(:,:,:,:,:,:,J,:) +Uy(:,:,:,:,:,:,Jm,:));   %\delta^2
    dmuz(:,:,:,:,:,:,:,K,1) = (Uz(:,:,:,:,:,:,:,Kp) ...
                              -Uz(:,:,:,:,:,:,:,Km))/2; %\mu\delta 
    dmuz(:,:,:,:,:,:,:,K,2) = (Uz(:,:,:,:,:,:,:,Kp) ...
       -2*Uz(:,:,:,:,:,:,:,K) +Uz(:,:,:,:,:,:,:,Km));   %\delta^2
  end% if stag
%{
\end{matlab}
Recursively take $\delta^2$ of these to form successively higher order
centred differences in all three space directions.
\begin{matlab}
%}
   for k = 3:ordCC    
    dmux(:,:,:,:,:,I,:,:,k) =     dmux(:,:,:,:,:,Ip,:,:,k-2) ...
    -2*dmux(:,:,:,:,:,I,:,:,k-2) +dmux(:,:,:,:,:,Im,:,:,k-2);    
    dmuy(:,:,:,:,:,:,J,:,k) =     dmuy(:,:,:,:,:,:,Jp,:,k-2) ...
    -2*dmuy(:,:,:,:,:,:,J,:,k-2) +dmuy(:,:,:,:,:,:,Jm,:,k-2);
    dmuz(:,:,:,:,:,:,:,K,k) =     dmuz(:,:,:,:,:,:,:,Kp,k-2) ...
    -2*dmuz(:,:,:,:,:,:,:,K,k-2) +dmuz(:,:,:,:,:,:,:,Km,k-2);
  end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet face values for
each patch \cite[]{Roberts06d, Bunder2013b}, using the weights
pre-computed by \verb|configPatches3()|. Here interpolate to
specified order.

For the case where next-to-face values interpolate to the
opposite face-values: when we have an ensemble of
configurations, different configurations might be coupled to
each other, as specified by \verb|patches.le|,
\verb|patches.ri|, \verb|patches.to|, \verb|patches.bo|,
\verb|patches.fr| and \verb|patches.ba|.
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two faces
u(nx,2:(ny-1),2:(nz-1),:,patches.ri,I,:,:) ...
  = Ux(1,:,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,1),-8).*dmux(1,:,:,:,:,:,:,:,:) ,9);  
u(1 ,2:(ny-1),2:(nz-1),:,patches.le,I,:,:) ...
  = Ux(k,:,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,1),-8).*dmux(k,:,:,:,:,:,:,:,:) ,9);
u(2:(nx-1),ny,2:(nz-1),:,patches.to,:,J,:) ...
  = Uy(:,1,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,2),-8).*dmuy(:,1,:,:,:,:,:,:,:) ,9);
u(2:(nx-1),1 ,2:(nz-1),:,patches.bo,:,J,:) ...
  = Uy(:,k,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,2),-8).*dmuy(:,k,:,:,:,:,:,:,:) ,9);
u(2:(nx-1),2:(ny-1),nz,:,patches.fr,:,:,K) ...
  = Uz(:,:,1,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,3),-8).*dmuz(:,:,1,:,:,:,:,:,:) ,9);
u(2:(nx-1),2:(ny-1),1 ,:,patches.ba,:,:,K) ...
  = Uz(:,:,k,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,3),-8).*dmuz(:,:,k,:,:,:,:,:,:) ,9);
%{
\end{matlab}



\paragraph{Case of spectral interpolation}
Assumes the domain is macro-periodic.
\begin{matlab}
%}
else% spectral interpolation
%{
\end{matlab}
We interpolate in terms of the patch index, $j$~say, not
directly in space. As the macroscale fields are $N$-periodic
in the patch index~$j$, the macroscale Fourier transform
writes the centre-patch values as $U_j=\sum_{k}C_ke^{ik2\pi
j/N}$. Then the face-patch values $U_{j\pm r}
=\sum_{k}C_ke^{ik2\pi/N(j\pm r)} =\sum_{k}C'_ke^{ik2\pi
j/N}$ where $C'_k =C_ke^{ikr2\pi/N}$. For $N$~patches we
resolve `wavenumbers' $|k|<N/2$, so set row vector
$\verb|ks| =k2\pi/N$ for `wavenumbers' $\mathcode`\,="213B
k=(0,1, \ldots, k_{\max}, -k_{\max}, \ldots, -1)$ for
odd~$N$, and $\mathcode`\,="213B k=(0,1, \ldots, k_{\max},
\pm(k_{\max}+1) -k_{\max}, \ldots, -1)$ for even~$N$.

Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches3|
tests there are an even number of patches). Then the
patch-ratio is effectively halved. The patch faces are near
the middle of the gaps and swapped.
\begin{matlab}
%}
 if patches.stag % transform by doubling the number of fields
 error('staggered grid not yet implemented??')
   v=nan(size(u)); % currently to restore the shape of u
   u=cat(3,u(:,1:2:nPatch,:),u(:,2:2:nPatch,:));
   stagShift=reshape(0.5*[ones(nVars,1);-ones(nVars,1)],1,1,[]);
   iV=[nVars+1:2*nVars 1:nVars]; % scatter interp to alternate field
   r=r/2;           % ratio effectively halved
   nPatch=nPatch/2; % halve the number of patches
   nVars=nVars*2;   % double the number of fields
 else % the values for standard spectral
    stagShift = 0;  
    iV = 1:nVars;
 end
%{
\end{matlab}
Now set wavenumbers in the three directions into three vectors
at the correct dimension.  In the case of even~$N$ these
compute the $+$-case for the highest wavenumber zig-zag
mode, $\mathcode`\,="213B k=(0,1, \ldots, k_{\max},
+(k_{\max}+1) -k_{\max}, \ldots, -1)$.
\begin{matlab}
%}
  kMax = floor((Nx-1)/2); 
  krx = shiftdim( rx*2*pi/Nx*(mod((0:Nx-1)+kMax,Nx)-kMax) ,-4);
  kMay = floor((Ny-1)/2); 
  kry = shiftdim( ry*2*pi/Ny*(mod((0:Ny-1)+kMay,Ny)-kMay) ,-5);
  kMaz = floor((Nz-1)/2); 
  krz = shiftdim( rz*2*pi/Nz*(mod((0:Nz-1)+kMaz,Nz)-kMaz) ,-6);
%{
\end{matlab}

Compute the Fourier transform of the patch values on the
centre-planes for all the fields.  Unless doing patch-edgy
interpolation when FT the next-to-face values.  If there are
an even number of points, then if complex, treat as positive
wavenumber, but if real, treat as cosine. When using an
ensemble of configurations, different configurations might
be coupled to each other, as specified by \verb|patches.le|,
\verb|patches.ri|, \verb|patches.to|, \verb|patches.bo|,
\verb|patches.fr| and \verb|patches.ba|.
\begin{matlab}
%}
% indices of interior
ix=(2:nx-1)';  iy=2:ny-1;  iz=shiftdim(2:nz-1,-1); 
if ~patches.EdgyInt
     Cle = fft(fft(fft( u(i0,iy,iz,:,:,:,:,:) ...
         ,[],6),[],7),[],8); 
     Cbo = fft(fft(fft( u(ix,j0,iz,:,:,:,:,:) ...
         ,[],6),[],7),[],8); 
     Cba = fft(fft(fft( u(ix,iy,k0,:,:,:,:,:) ...
         ,[],6),[],7),[],8); 
     Cri = Cle;    Cto = Cbo;    Cfr = Cba;
else 
     Cle = fft(fft(fft( u(   2,iy,iz ,:,patches.le,:,:,:) ...
         ,[],6),[],7),[],8);
     Cri = fft(fft(fft( u(nx-1,iy,iz ,:,patches.ri,:,:,:) ...
         ,[],6),[],7),[],8);
     Cbo = fft(fft(fft( u(ix,2   ,iz ,:,patches.bo,:,:,:) ...
         ,[],6),[],7),[],8);
     Cto = fft(fft(fft( u(ix,ny-1,iz ,:,patches.to,:,:,:) ...
         ,[],6),[],7),[],8);
     Cba = fft(fft(fft( u(ix,iy,2    ,:,patches.ba,:,:,:) ...
         ,[],6),[],7),[],8);
     Cfr = fft(fft(fft( u(ix,iy,nz-1 ,:,patches.fr,:,:,:) ...
         ,[],6),[],7),[],8);
end     
%{
\end{matlab}
Now invert the triple Fourier transforms to complete
interpolation.   (Should stagShift be multiplied by
rx/ry/rz??) Enforce reality when appropriate. 
\begin{matlab}
%}
u(nx,iy,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cle.*exp(1i*(stagShift+krx))  ,[],6),[],7),[],8) );
u( 1,iy,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cri.*exp(1i*(stagShift-krx))  ,[],6),[],7),[],8) );
u(ix,ny,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cbo.*exp(1i*(stagShift+kry))  ,[],6),[],7),[],8) );
u(ix, 1,iz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cto.*exp(1i*(stagShift-kry))  ,[],6),[],7),[],8) );
u(ix,iy,nz,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cba.*exp(1i*(stagShift+krz))  ,[],6),[],7),[],8) );
u(ix,iy, 1,:,:,:,:,:) = uclean( ifft(ifft(ifft( ...
    Cfr.*exp(1i*(stagShift-krz))  ,[],6),[],7),[],8) );
end% if spectral 
%{
\end{matlab}
Nan the values in corners and edges, of every 3D patch.
\begin{matlab}
%}
u(:,[1 ny],[1 nz],:,:,:,:,:)=nan; 
u([1 nx],:,[1 nz],:,:,:,:,:)=nan; 
u([1 nx],[1 ny],:,:,:,:,:,:)=nan; 
end% function patchEdgeInt3
%{
\end{matlab}
Fin, returning the 8D array of field values with
interpolated faces. 
\end{devMan} 
%}
