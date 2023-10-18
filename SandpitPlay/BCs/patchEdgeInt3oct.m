% patchEdgeInt3old() provides the interpolation across 3D space
% for 3D patches of simulations of a smooth lattice system
% such as PDE discretisations.  AJR, Aug 2020 -- 12 Apr 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt3old()}: sets 3D patch
face values from 3D macroscale interpolation}
\label{sec:patchEdgeInt3old}


Couples 3D patches across 3D space by computing their face
values via macroscale interpolation.  Assumes patch face
values are determined by macroscale interpolation of the
patch centre-plane values \cite[]{Roberts2011a,
Bunder2019d}, or patch next-to-face values which appears
better \cite[]{Bunder2020a}.  This function is primarily
used by \verb|patchSys3()| but is also useful for user
graphics. \footnote{Script \texttt{patchEdgeInt3oldtest.m}
verifies this code.}

Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|), or otherwise via the global
struct~\verb|patches|.
\begin{matlab}
%}
function u = patchEdgeInt3old(u,patches)
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
\times1 \times \verb|nPatch1| \times1 \times1 $ array of the
spatial locations~$x_{iI}$ of the microscale grid points in
every patch. Currently it \emph{must} be an equi-spaced
lattice on the microscale index~$i$, but may be variable
spaced in macroscale index~$I$.

\item \verb|.y| is similarly $1\times \verb|nSubP2| \times1
\times1 \times1 \times1 \times \verb|nPatch2| \times1$ array
of the spatial locations~$y_{jJ}$ of the microscale grid
points in every patch. Currently it \emph{must} be an
equi-spaced lattice on the microscale index~$j$, but may be
variable spaced in macroscale index~$J$.

\item \verb|.z| is similarly $1 \times1 \times \verb|nSubP3|
\times1 \times1 \times1 \times1 \times \verb|nPatch3|$ array
of the spatial locations~$z_{kK}$ of the microscale grid
points in every patch. Currently it \emph{must} be an
equi-spaced lattice on the microscale index~$k$, but may be
variable spaced in macroscale index~$K$.

\item \verb|.ordCC| is order of interpolation, currently
only $\{0,2,4,\ldots\}$

\item \verb|.periodic| indicates whether macroscale is
periodic domain, or alternatively that the macroscale has
left, right, top, bottom, front and back boundaries so
interpolation is via divided differences. 

\item \verb|.stag| in $\{0,1\}$ is one for staggered grid
(alternating) interpolation.  Currently must be zero.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the coupling
coefficients for finite width interpolation in each of the
$x,y,z$-directions---when invoking a periodic domain.

\item \verb|.EdgyInt|, true/false, for determining
patch-edge values by interpolation: true, from opposite-edge
next-to-edge values (often preserves symmetry); false, from
centre cross-patch values (near original scheme).

\item \verb|.nEdge|, three elements, the width of edge
values set by interpolation at the \(x,y,z\)-face regions, 
respectively, of each patch (default is one all 
\(x,y,z\)-faces).

\item \verb|.nEnsem| the number of realisations in the
ensemble.

\item \verb|.parallel| whether serial or parallel.

\end{itemize}
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 8D array, $\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nSubP3| \cdot \verb|nVars| \cdot
\verb|nEnsem| \cdot \verb|nPatch1| \cdot \verb|nPatch2|
\cdot \verb|nPatch3|$, of the fields with face values set by
interpolation.
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
  ,'patchEdgeInt3old: input u has wrong size for parameters')
u = reshape(u,[nx ny nz nVars nEnsem Nx Ny Nz]);
%{
\end{matlab}

For the moment assume the physical domain is either
macroscale periodic or macroscale rectangle so that the
coupling formulas are simplest.  These index vectors point
to patches and, if periodic, their six immediate neighbours.
\begin{matlab}
%}
I=1:Nx; Ip=mod(I,Nx)+1; Im=mod(I-2,Nx)+1;
J=1:Ny; Jp=mod(J,Ny)+1; Jm=mod(J-2,Ny)+1;
K=1:Nz; Kp=mod(K,Nz)+1; Km=mod(K-2,Nz)+1;
%{
\end{matlab}

\paragraph{Implement multiple width edges by folding}
Subsample~\(x,y,z\) coordinates, noting it is only differences
that count \emph{and} the microgrid~\(x,y,z\) spacing must be
uniform.
\begin{matlab}
%}
x = patches.x;
y = patches.y; 
z = patches.z;
if mean(patches.nEdge)>1
  mx = patches.nEdge(1);
  my = patches.nEdge(2);
  mz = patches.nEdge(3);
  x = x(1:mx:nx,:,:,:,:,:,:,:);
  y = y(:,1:my:ny,:,:,:,:,:,:);
  z = z(:,:,1:mz:nz,:,:,:,:,:);
  nx = nx/mx;
  ny = ny/my;
  nz = nz/mz;
  u = reshape(u,mx,nx,my,ny,mz,nz,nVars,nEnsem,Nx,Ny,Nz);
  nVars = nVars*mx*my*mz;
  u = reshape( permute(u,[2:2:6 1:2:5 7:11]) ...
             ,nx,ny,nz,nVars,nEnsem,Nx,Ny,Nz);
end%if patches.nEdge
%{
\end{matlab}

The centre of each patch (as \verb|nx|, \verb|ny| and
\verb|nz| are odd for centre-patch interpolation) is at
indices
\begin{matlab}
%}
i0 = round((nx+1)/2);
j0 = round((ny+1)/2);
k0 = round((nz+1)/2);
%{
\end{matlab}



\subsection{Periodic macroscale interpolation schemes}
\begin{matlab}
%}
if patches.periodic
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


\subsubsection{Lagrange interpolation gives patch-face values}
Compute centred differences of the mid-patch values for the
macro-interpolation, of all fields.  Here the domain is
macro-periodic.
\begin{matlab}
%}
ordCC = patches.ordCC;
if ordCC>0 % then finite-width polynomial interpolation
%{
\end{matlab}
Interpolate the three directions in succession, in this way
we naturally fill-in face-edge and corner values. Start with
\(x\)-direction, and give most documentation for that case
as the others are essentially the same.


\paragraph{\(x\)-normal face values} The patch-edge values
are either interpolated from the next-to-edge-face values,
or from the centre-cross-plane values (not the patch-centre
value itself as that seems to have worse properties in
general).  Have not yet implemented core averages.
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    U = u([2 nx-1],2:(ny-1),2:(nz-1),:,:,I,J,K);
  else % interpolate centre-cross values
    U = u(i0,2:(ny-1),2:(nz-1),:,:,I,J,K);
  end;%if patches.EdgyInt
%{
\end{matlab}
Just in case any last array dimension(s) are one, we force a
padding of the sizes, then adjoin the extra dimension for
the subsequent array of differences.
\begin{matlab}
%}
szUO=size(U); szUO=[szUO ones(1,8-length(szUO)) ordCC];
%{
\end{matlab}
Use finite difference formulas for the interpolation, so
store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in these arrays.  When parallel, in order
to preserve the distributed array structure we use an index
at the end for the differences.
\begin{matlab}
%}
  if ~patches.parallel, dmu = zeros(szUO); % 9D
  else   dmu = zeros(szUO,patches.codist); % 9D
  end%if patches.parallel
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
%    dmux(:,:,:,:,:,I,:,:,1) = (Ux(:,:,:,:,:,Ip,:,:) +Ux(:,:,:,:,:,Im,:,:))/2; % \mu
%    dmux(:,:,:,:,:,I,:,:,2) = (Ux(:,:,:,:,:,Ip,:,:) -Ux(:,:,:,:,:,Im,:,:)); % \delta
%    Ip = Ip(Ip); Im = Im(Im); % increase shifts to \pm2
%    dmuy(:,:,:,:,:,:,J,:,1) = (Ux(:,:,:,:,:,:,Jp,:)+Ux(:,:,:,:,:,:,Jm,:))/2; % \mu
%    dmuy(:,:,:,:,:,:,J,:,2) = (Ux(:,:,:,:,:,:,Jp,:)-Ux(:,:,:,:,:,:,Jm,:)); % \delta
%    Jp = Jp(Jp); Jm = Jm(Jm); % increase shifts to \pm2
%    dmuz(:,:,:,:,:,:,:,K,1) = (Ux(:,:,:,:,:,:,:,Kp)+Ux(:,:,:,:,:,:,:,Km))/2; % \mu
%    dmuz(:,:,:,:,:,:,:,K,2) = (Ux(:,:,:,:,:,:,:,Kp)-Ux(:,:,:,:,:,:,:,Km)); % \delta
%    Kp = Kp(Kp); Km = Km(Km); % increase shifts to \pm2
  else  %disp('starting standard interpolation')  
    dmu(:,:,:,:,:,I,:,:,1) = (U(:,:,:,:,:,Ip,:,:) ...
                             -U(:,:,:,:,:,Im,:,:))/2; %\mu\delta 
    dmu(:,:,:,:,:,I,:,:,2) = (U(:,:,:,:,:,Ip,:,:) ...
       -2*U(:,:,:,:,:,I,:,:) +U(:,:,:,:,:,Im,:,:));   %\delta^2    
  end% if stag
%{
\end{matlab}
Recursively take $\delta^2$ of these to form successively
higher order centred differences in space.
\begin{matlab}
%}
  for k = 3:ordCC    
    dmu(:,:,:,:,:,I,:,:,k) =     dmu(:,:,:,:,:,Ip,:,:,k-2) ...
    -2*dmu(:,:,:,:,:,I,:,:,k-2) +dmu(:,:,:,:,:,Im,:,:,k-2);    
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
  = U(1,:,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,1),-8).*dmu(1,:,:,:,:,:,:,:,:) ,9);  
u(1 ,2:(ny-1),2:(nz-1),:,patches.le,I,:,:) ...
  = U(k,:,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,1),-8).*dmu(k,:,:,:,:,:,:,:,:) ,9);
%{
\end{matlab}




\paragraph{\(y\)-normal face values} Interpolate from either
the next-to-edge-face values, or the centre-cross-plane
values.
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    U = u(:,[2 ny-1],2:(nz-1),:,:,I,J,K);
  else % interpolate centre-cross values
    U = u(:,j0,2:(nz-1),:,:,I,J,K);
  end;%if patches.EdgyInt
%{
\end{matlab}
Adjoin extra dimension for the array of differences.
\begin{matlab}
%}
szUO=size(U); szUO=[szUO ones(1,8-length(szUO)) ordCC];
%{
\end{matlab}
Store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in this array.
\begin{matlab}
%}
  if ~patches.parallel, dmu = zeros(szUO); % 9D
  else   dmu = zeros(szUO,patches.codist); % 9D
  end%if patches.parallel
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
  else  %disp('starting standard interpolation')  
    dmu(:,:,:,:,:,:,J,:,1) = (U(:,:,:,:,:,:,Jp,:) ...
                             -U(:,:,:,:,:,:,Jm,:))/2; %\mu\delta 
    dmu(:,:,:,:,:,:,J,:,2) = (U(:,:,:,:,:,:,Jp,:) ...
       -2*U(:,:,:,:,:,:,J,:) +U(:,:,:,:,:,:,Jm,:));   %\delta^2
  end% if stag
%{
\end{matlab}
Recursively take $\delta^2$.
\begin{matlab}
%}
  for k = 3:ordCC    
    dmu(:,:,:,:,:,:,J,:,k) =     dmu(:,:,:,:,:,:,Jp,:,k-2) ...
    -2*dmu(:,:,:,:,:,:,J,:,k-2) +dmu(:,:,:,:,:,:,Jm,:,k-2);
  end
%{
\end{matlab}
Interpolate macro-values using the weights pre-computed by
\verb|configPatches3()|.  An ensemble of configurations may
have cross-coupling.
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two faces
u(:,ny,2:(nz-1),:,patches.to,:,J,:) ...
  = U(:,1,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,2),-8).*dmu(:,1,:,:,:,:,:,:,:) ,9);
u(:,1 ,2:(nz-1),:,patches.bo,:,J,:) ...
  = U(:,k,:,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,2),-8).*dmu(:,k,:,:,:,:,:,:,:) ,9);
%{
\end{matlab}



\paragraph{\(z\)-normal face values} Interpolate from either
the next-to-edge-face values, or the centre-cross-plane
values.
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    U = u(:,:,[2 nz-1],:,:,I,J,K);
  else % interpolate centre-cross values
    U = u(:,:,k0,:,:,I,J,K);
  end;%if patches.EdgyInt
%{
\end{matlab}
Adjoin extra dimension for the array of differences.
\begin{matlab}
%}
szUO=size(U); szUO=[szUO ones(1,8-length(szUO)) ordCC];
%{
\end{matlab}
Store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in this array.
\begin{matlab}
%}
  if ~patches.parallel, dmu = zeros(szUO); % 9D
  else   dmu = zeros(szUO,patches.codist); % 9D
  end%if patches.parallel
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
  else  %disp('starting standard interpolation')  
    dmu(:,:,:,:,:,:,:,K,1) = (U(:,:,:,:,:,:,:,Kp) ...
                             -U(:,:,:,:,:,:,:,Km))/2; %\mu\delta 
    dmu(:,:,:,:,:,:,:,K,2) = (U(:,:,:,:,:,:,:,Kp) ...
       -2*U(:,:,:,:,:,:,:,K) +U(:,:,:,:,:,:,:,Km));   %\delta^2
  end% if stag
%{
\end{matlab}
Recursively take $\delta^2$.
\begin{matlab}
%}
  for k = 3:ordCC    
    dmu(:,:,:,:,:,:,:,K,k) =     dmu(:,:,:,:,:,:,:,Kp,k-2) ...
    -2*dmu(:,:,:,:,:,:,:,K,k-2) +dmu(:,:,:,:,:,:,:,Km,k-2);
  end
%{
\end{matlab}
Interpolate macro-values using the weights pre-computed by
\verb|configPatches3()|.  An ensemble of configurations may
have cross-coupling.
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two faces
u(:,:,nz,:,patches.fr,:,:,K) ...
  = U(:,:,1,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,3),-8).*dmu(:,:,1,:,:,:,:,:,:) ,9);
u(:,:,1 ,:,patches.ba,:,:,K) ...
  = U(:,:,k,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,3),-8).*dmu(:,:,k,:,:,:,:,:,:) ,9);
%{
\end{matlab}



\subsubsection{Case of spectral interpolation}
Assumes the domain is macro-periodic.
\begin{matlab}
%}
else% patches.ordCC<=0, spectral interpolation
%{
\end{matlab}
We interpolate in terms of the patch index, $I$~say, not
directly in space. As the macroscale fields are $N$-periodic
in the patch index~$I$, the macroscale Fourier transform
writes the centre-patch values as $U_I=\sum_{k}C_ke^{ik2\pi
I/N}$. Then the face-patch values $U_{I\pm r}
=\sum_{k}C_ke^{ik2\pi/N(I\pm r)} =\sum_{k}C'_ke^{ik2\pi
I/N}$ where $C'_k =C_ke^{ikr2\pi/N}$. For $N$~patches we
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
 end%if patches.stag
%{
\end{matlab}
Interpolate the three directions in succession, in this way
we naturally fill-in face-edge and corner values. Start with
\(x\)-direction, and give most documentation for that case
as the others are essentially the same. Need these indices
of patch interior.
\begin{matlab}
%}
ix = 2:nx-1;   iy = 2:ny-1;   iz = 2:nz-1; 
%{
\end{matlab}



\paragraph{\(x\)-normal face values} Now set wavenumbers
into a vector at the correct dimension.  In the case of
even~$N$ these compute the $+$-case for the highest
wavenumber zig-zag mode, $\mathcode`\,="213B k=(0,1, \ldots,
k_{\max}, +(k_{\max}+1) -k_{\max}, \ldots, -1)$.
\begin{matlab}
%}
  kMax = floor((Nx-1)/2); 
  kr = shiftdim( rx*2*pi/Nx*(mod((0:Nx-1)+kMax,Nx)-kMax) ,-4);
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
if ~patches.EdgyInt
     Cm = fft( u(i0,iy,iz,:,:,:,:,:) ,[],6); 
     Cp = Cm;
else 
     Cm = fft( u(   2,iy,iz ,:,patches.le,:,:,:) ,[],6);
     Cp = fft( u(nx-1,iy,iz ,:,patches.ri,:,:,:) ,[],6);
end%if ~patches.EdgyInt  
%{
\end{matlab}
Now invert the Fourier transforms to complete interpolation.
Enforce reality when appropriate. 
\begin{matlab}
%}
u(nx,iy,iz,:,:,:,:,:) = uclean( ifft( ...
    Cm.*exp(1i*(stagShift+kr))  ,[],6) );
u( 1,iy,iz,:,:,:,:,:) = uclean( ifft( ...
    Cp.*exp(1i*(stagShift-kr))  ,[],6) );
%{
\end{matlab}



\paragraph{\(y\)-normal face values} Set wavenumbers into a
vector.
\begin{matlab}
%}
  kMax = floor((Ny-1)/2); 
  kr = shiftdim( ry*2*pi/Ny*(mod((0:Ny-1)+kMax,Ny)-kMax) ,-5);
%{
\end{matlab}
Compute the Fourier transform of the patch values on the
centre-planes for all the fields.  
\begin{matlab}
%}
if ~patches.EdgyInt
     Cm = fft( u(:,j0,iz,:,:,:,:,:) ,[],7); 
     Cp = Cm;
else 
     Cm = fft( u(:,2   ,iz ,:,patches.bo,:,:,:) ,[],7);
     Cp = fft( u(:,ny-1,iz ,:,patches.to,:,:,:) ,[],7);
end%if ~patches.EdgyInt  
%{
\end{matlab}
Invert the Fourier transforms to complete interpolation. 
\begin{matlab}
%}
u(:,ny,iz,:,:,:,:,:) = uclean( ifft( ...
    Cm.*exp(1i*(stagShift+kr))  ,[],7) );
u(:, 1,iz,:,:,:,:,:) = uclean( ifft( ...
    Cp.*exp(1i*(stagShift-kr))  ,[],7) );
%{
\end{matlab}



\paragraph{\(z\)-normal face values} Set wavenumbers into a
vector.
\begin{matlab}
%}
  kMax = floor((Nz-1)/2); 
  kr = shiftdim( rz*2*pi/Nz*(mod((0:Nz-1)+kMax,Nz)-kMax) ,-6);
%{
\end{matlab}
Compute the Fourier transform of the patch values on the
centre-planes for all the fields.  
\begin{matlab}
%}
if ~patches.EdgyInt
     Cm = fft( u(:,:,k0,:,:,:,:,:) ,[],8); 
     Cp = Cm;
else 
     Cm = fft( u(:,:,2    ,:,patches.ba,:,:,:) ,[],8);
     Cp = fft( u(:,:,nz-1 ,:,patches.fr,:,:,:) ,[],8);
end%if ~patches.EdgyInt  
%{
\end{matlab}
Invert the Fourier transforms to complete interpolation. 
\begin{matlab}
%}
u(:,:,nz,:,:,:,:,:) = uclean( ifft( ...
    Cm.*exp(1i*(stagShift+kr))  ,[],8) );
u(:,:, 1,:,:,:,:,:) = uclean( ifft( ...
    Cp.*exp(1i*(stagShift-kr))  ,[],8) );
%{
\end{matlab}

\begin{matlab}
%}
end% if ordCC>0  
%{
\end{matlab}





\subsection{Non-periodic macroscale interpolation}
\begin{matlab}
%}
else% patches.periodic false
assert(~patches.stag, ...
'not yet implemented staggered grids for non-periodic')
%{
\end{matlab}
Determine the order of interpolation~\verb|px|, \verb|py|
and~\verb|pz| (potentially different in the different
directions!), and hence size of the (forward) divided
difference tables in~\verb|F|~(9D) for interpolating to
left/right, top/bottom, and front/back faces. Because of the
product-form of the patch grid, and because we are doing
\emph{only} either edgy interpolation or cross-patch
interpolation (\emph{not} just the centre patch value), the
interpolations are all 1D interpolations.
\begin{matlab}
%}
if patches.ordCC<1
     px = Nx-1;  py = Ny-1;  pz = Nz-1;
else px = min(patches.ordCC,Nx-1); 
     py = min(patches.ordCC,Ny-1); 
     pz = min(patches.ordCC,Nz-1); 
end
% interior indices of faces  (ix n/a)
ix=2:nx-1;  iy=2:ny-1;  iz=2:nz-1; 
%{
\end{matlab}


\subsubsection{\(x\)-direction values}
Set function values in first `column' of the tables for
every variable and across ensemble.  For~\verb|EdgyInt|, the
`reversal' of the next-to-face values are because their
values are to interpolate to the opposite face of each
patch. \todo{Have no plans to implement core averaging as yet.}
\begin{matlab}
%}
  F = nan(patches.EdgyInt+1,ny-2,nz-2,nVars,nEnsem,Nx,Ny,Nz,px+1);
  if patches.EdgyInt % interpolate next-to-face values
    F(:,:,:,:,:,:,:,:,1) = u([nx-1 2],iy,iz,:,:,:,:,:);
    X = x([nx-1 2],:,:,:,:,:,:,:);
  else % interpolate mid-patch cross-patch values 
    F(:,:,:,:,:,:,:,:,1) = u(i0,iy,iz,:,:,:,:,:);
    X = x(i0,:,:,:,:,:,:,:);
  end%if patches.EdgyInt
%{
\end{matlab}

\paragraph{Form tables of divided differences} Compute
tables of (forward) divided differences
\cite[e.g.,][]{DividedDifferences} for every variable, and
across ensemble, and in both directions, and for all three
types of faces (left/right, top/bottom, and front/back).
Recursively find all divided differences in the respective
direction.
\begin{matlab}
%}
for q = 1:px
  i = 1:Nx-q;
  F(:,:,:,:,:,i,:,:,q+1) ...
  = ( F(:,:,:,:,:,i+1,:,:,q)-F(:,:,:,:,:,i,:,:,q)) ...
   ./(X(:,:,:,:,:,i+q,:,:)  -X(:,:,:,:,:,i,:,:));
end
%{
\end{matlab}

\paragraph{Interpolate with divided differences} Now
interpolate to find the face-values on left/right faces
at~\verb|Xface| for every interior~\verb|Y,Z|.
\begin{matlab}
%}
Xface = x([1 nx],:,:,:,:,:,:,:);
%{
\end{matlab}
Code Horner's recursive evaluation of the interpolation
polynomials.  Indices~\verb|i| are those of the left face of
each interpolation stencil, because the table is of forward
differences.  This alternative: the case of order~\(p_x\),
\(p_y\) and~\(p_z\) interpolation across the domain,
asymmetric near the boundaries of the rectangular domain.
\begin{matlab}
%}
  i = max(1,min(1:Nx,Nx-ceil(px/2))-floor(px/2));
  Uface = F(:,:,:,:,:,i,:,:,px+1);
  for q = px:-1:1
    Uface = F(:,:,:,:,:,i,:,:,q) ...
    +(Xface-X(:,:,:,:,:,i+q-1,:,:)).*Uface;
  end
%{
\end{matlab}

Finally, insert face values into the array of field values,
using the required ensemble shifts.  
\begin{matlab}
%}
u(1 ,iy,iz,:,patches.le,:,:,:) = Uface(1,:,:,:,:,:,:,:);
u(nx,iy,iz,:,patches.ri,:,:,:) = Uface(2,:,:,:,:,:,:,:);
%{
\end{matlab}


\subsubsection{\(y\)-direction values}
Set function values in first `column' of the tables for
every variable and across ensemble.
\begin{matlab}
%}
  F = nan(nx,patches.EdgyInt+1,nz-2,nVars,nEnsem,Nx,Ny,Nz,py+1);
  if patches.EdgyInt % interpolate next-to-face values
    F(:,:,:,:,:,:,:,:,1) = u(:,[ny-1 2],iz,:,:,:,:,:);
    Y = y(:,[ny-1 2],:,:,:,:,:,:);
  else % interpolate mid-patch cross-patch values 
    F(:,:,:,:,:,:,:,:,1) = u(:,j0,iz,:,:,:,:,:);
    Y = y(:,j0,:,:,:,:,:,:);
  end%if patches.EdgyInt
%{
\end{matlab}
Form tables of divided differences.
\begin{matlab}
%}
for q = 1:py
  j = 1:Ny-q;
  F(:,:,:,:,:,:,j,:,q+1) ...
  = ( F(:,:,:,:,:,:,j+1,:,q)-F(:,:,:,:,:,:,j,:,q)) ...
   ./(Y(:,:,:,:,:,:,j+q,:)  -Y(:,:,:,:,:,:,j,:));
end
%{
\end{matlab}
Interpolate to find the top/bottom faces~\verb|Yface| for
every~\(x\) and interior~\(z\).
\begin{matlab}
%}
Yface = y(:,[1 ny],:,:,:,:,:,:);
%{
\end{matlab}
Code Horner's recursive evaluation of the interpolation
polynomials.  Indices~\verb|j| are those of the bottom face
of each interpolation stencil, because the table is of
forward differences.  
\begin{matlab}
%}
  j = max(1,min(1:Ny,Ny-ceil(py/2))-floor(py/2));
  Uface = F(:,:,:,:,:,:,j,:,py+1);
  for q = py:-1:1
    Uface = F(:,:,:,:,:,:,j,:,q) ...
    +(Yface-Y(:,:,:,:,:,:,j+q-1,:)).*Uface;
  end
%{
\end{matlab}

Finally, insert face values into the array of field values,
using the required ensemble shifts.  
\begin{matlab}
%}
u(:,1 ,iz,:,patches.bo,:,:,:) = Uface(:,1,:,:,:,:,:,:);
u(:,ny,iz,:,patches.to,:,:,:) = Uface(:,2,:,:,:,:,:,:);
%{
\end{matlab}


\subsubsection{\(z\)-direction values}
Set function values in first `column' of the tables for
every variable and across ensemble.  
\begin{matlab}
%}
  F = nan(nx,ny,patches.EdgyInt+1,nVars,nEnsem,Nx,Ny,Nz,pz+1);
  if patches.EdgyInt % interpolate next-to-face values
    F(:,:,:,:,:,:,:,:,1) = u(:,:,[nz-1 2],:,:,:,:,:);
    Z = z(:,:,[nz-1 2],:,:,:,:,:);
  else % interpolate mid-patch cross-patch values 
    F(:,:,:,:,:,:,:,:,1) = u(:,:,k0,:,:,:,:,:);
    Z = z(:,:,k0,:,:,:,:,:);
  end%if patches.EdgyInt
%{
\end{matlab}
Form tables of divided differences.
\begin{matlab}
%}
for q = 1:pz
  k = 1:Nz-q;
  F(:,:,:,:,:,:,:,k,q+1) ...
  = ( F(:,:,:,:,:,:,:,k+1,q)-F(:,:,:,:,:,:,:,k,q)) ...
   ./(Z(:,:,:,:,:,:,:,k+q)  -Z(:,:,:,:,:,:,:,k));
end
%{
\end{matlab}
Interpolate to find the face-values on front/back
faces~\verb|Zface| for every~\(x,y\).
\begin{matlab}
%}
Zface = z(:,:,[1 nz],:,:,:,:,:);
%{
\end{matlab}
Code Horner's recursive evaluation of the interpolation
polynomials.  Indices~\verb|k| are those of the bottom face
of each interpolation stencil, because the table is of
forward differences.  
\begin{matlab}
%}
  k = max(1,min(1:Nz,Nz-ceil(pz/2))-floor(pz/2));
  Uface = F(:,:,:,:,:,:,:,k,pz+1);
  for q = pz:-1:1
    Uface = F(:,:,:,:,:,:,:,k,q) ...
    +(Zface-Z(:,:,:,:,:,:,:,k+q-1)).*Uface;
  end
%{
\end{matlab}

Finally, insert face values into the array of field values,
using the required ensemble shifts.  
\begin{matlab}
%}
u(:,:,1 ,:,patches.fr,:,:,:) = Uface(:,:,1,:,:,:,:,:);
u(:,:,nz,:,patches.ba,:,:,:) = Uface(:,:,2,:,:,:,:,:);
%{
\end{matlab}


\subsubsection{Optional NaNs for safety}
We want a user to set outer face values on the extreme
patches according to the microscale boundary conditions that
hold at the extremes of the domain.  Consequently, unless testing, override
their computed interpolation values with~\verb|NaN|.
\begin{matlab}
%}
if isfield(patches,'intTest')&&patches.intTest
else % usual case
    u( 1,:,:,:,:, 1,:,:) = nan;
    u(nx,:,:,:,:,Nx,:,:) = nan;
    u(:, 1,:,:,:,:, 1,:) = nan;
    u(:,ny,:,:,:,:,Ny,:) = nan;
    u(:,:, 1,:,:,:,:, 1) = nan;
    u(:,:,nz,:,:,:,:,Nz) = nan;
end%if
%{
\end{matlab}

End of the non-periodic interpolation code.
\begin{matlab}
%}
end%if patches.periodic else
%{
\end{matlab}

\paragraph{Unfold multiple edges}  No need to restore~\(x,y,z\).
\begin{matlab}
%}
if mean(patches.nEdge)>1
  nVars = nVars/(mx*my*mz);
  u = reshape( u ,nx,ny,nz,mx,my,mz,nVars,nEnsem,Nx,Ny,Nz);
  nx = nx*mx;
  ny = ny*my;
  nz = nz*mz;
  u = reshape( permute(u,[4 1 5 2 6 3 7:11]) ...
             ,nx,ny,nz,nVars,nEnsem,Nx,Ny,Nz);
end%if patches.nEdge
%{
\end{matlab}

Fin, returning the 8D array of field values with
interpolated faces. 
\begin{matlab}
%}
end% function patchEdgeInt3old
%{
\end{matlab}
\end{devMan} 
%}
