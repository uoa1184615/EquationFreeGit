% patchEdgeInt2() provides the interpolation across 2D space
% for 2D patches of simulations of a lattice system such as
% a PDE discretisation.  AJR, Nov 2018 -- 12 Apr 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt2()}: sets 2D patch
edge values from 2D macroscale interpolation}
\label{sec:patchEdgeInt2}


Couples 2D patches across 2D space by computing their edge
values via macroscale interpolation.  Research
\cite[]{Roberts2011a, Bunder2019d} indicates the patch
centre-values are sensible macroscale variables, and
macroscale interpolation of these determine patch-edge
values.  However, for computational homogenisation in
multi-D, interpolating patch next-to-edge values appears
better \cite[]{Bunder2020a}.  This function is primarily
used by \verb|patchSys2()| but is also useful for user
graphics. \footnote{Script \texttt{patchEdgeInt2test.m}
verifies this code.}
 
Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|), or otherwise via the global 
struct~\verb|patches|.
\begin{matlab}
%}
function u = patchEdgeInt2(u,patches)
if nargin<2, global patches, end
%{
\end{matlab}



\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector\slash array of length
$\verb|prod(nSubP)|  \cdot \verb|nVars| \cdot \verb|nEnsem|
\cdot \verb|prod(nPatch)|$ where there are $\verb|nVars|
\cdot \verb|nEnsem|$ field values at each of the points in
the $\verb|nSubP1| \cdot \verb|nSubP2| \cdot \verb|nPatch1|
\cdot \verb|nPatch2|$ multiscale spatial grid on the
$\verb|nPatch1| \cdot \verb|nPatch2|$ array of patches.

\item \verb|patches| a struct set by \verb|configPatches2()|
which includes the following information.
\begin{itemize}

\item \verb|.x| is $\verb|nSubP1| \times1 \times1 \times1
\times \verb|nPatch1| \times1 $ array of the spatial
locations~$x_{iI}$ of the microscale grid points in every
patch. Currently it \emph{must} be an equi-spaced lattice on
the microscale index~$i$, but may be variable spaced in
macroscale index~$I$.

\item \verb|.y| is similarly $1 \times \verb|nSubP2| \times1
\times1 \times1 \times \verb|nPatch2|$ array of the spatial
locations~$y_{jJ}$ of the microscale grid points in every
patch. Currently it \emph{must} be an equi-spaced lattice on
the microscale index~$j$, but may be variable spaced in
macroscale index~$J$.

\item \verb|.ordCC| is order of interpolation, currently 
only $\{0,2,4,\ldots\}$

\item \verb|.periodic| indicates whether macroscale is
periodic domain, or alternatively that the macroscale has
left, right, top and bottom boundaries so interpolation is
via divided differences. 

\item \verb|.stag| in $\{0,1\}$ is one for staggered grid
(alternating) interpolation.  Currently must be zero.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the coupling
coefficients for finite width interpolation in both the
$x,y$-directions---when invoking a periodic domain.

\item \verb|.EdgyInt|, true/false, for determining
patch-edge values by interpolation: true, from opposite-edge
next-to-edge values (often preserves symmetry); false, from
centre cross-patch values (near original scheme).

\item \verb|.nEdge|, two elements, the width of edge
values set by interpolation at the \(x,y\)-edge regions, 
respectively, of each patch (default is one for both 
\(x,y\)-edges).

\item \verb|.nEnsem| the number of realisations in the
ensemble.

\item \verb|.parallel| whether serial or parallel.

\end{itemize}
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 6D array, $\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nVars| \cdot \verb|nEnsem| \cdot
\verb|nPatch1| \cdot \verb|nPatch2|$, of the fields with
edge values set by interpolation.
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
[~,ny,~,~,~,Ny] = size(patches.y);
[nx,~,~,~,Nx,~] = size(patches.x);
nEnsem = patches.nEnsem;
nVars = round(numel(u)/numel(patches.x)/numel(patches.y)/nEnsem);
assert(numel(u) == nx*ny*Nx*Ny*nVars*nEnsem ...
  ,'patchEdgeInt2: input u has wrong size for parameters')
u = reshape(u,[nx ny nVars nEnsem Nx Ny ]);
%{
\end{matlab}

For the moment assume the physical domain is either
macroscale periodic or macroscale rectangle so that the
coupling formulas are simplest.  These index vectors point
to patches and, if periodic, their four immediate neighbours. 
\begin{matlab}
%}
I=1:Nx; Ip=mod(I,Nx)+1; Im=mod(I-2,Nx)+1;
J=1:Ny; Jp=mod(J,Ny)+1; Jm=mod(J-2,Ny)+1;
%{
\end{matlab}

\paragraph{Implement multiple width edges by folding}
Subsample~\(x,y\) coordinates, noting it is only differences
that count \emph{and} the microgrid~\(x,y\) spacing must be
uniform.
\begin{matlab}
%}
%x = patches.x;
%if patches.nEdge(1)>1
%  m = patches.nEdge(1);
%  x = x(1:m:nx,:,:,:,:,:);
%  nx = nx/m;
%  u = reshape(u,m,nx,ny,nVars,nEnsem,Nx,Ny);
%  nVars = nVars*m;
%  u = reshape( permute(u,[2:3 1 4:7]) ...
%             ,nx,ny,nVars,nEnsem,Nx,Ny);
%end%if patches.nEdge(1)
%y = patches.y;
%if patches.nEdge(2)>1
%  m = patches.nEdge(2);
%  y = y(:,1:m:ny,:,:,:,:);
%  ny = ny/m;
%  u = reshape(u,nx,m,ny,nVars,nEnsem,Nx,Ny);
%  nVars = nVars*m;
%  u = reshape( permute(u,[1 3 2 4:7]) ...
%             ,nx,ny,nVars,nEnsem,Nx,Ny);
%end%if patches.nEdge(2)
x = patches.x;
y = patches.y;
if mean(patches.nEdge)>1
  mx = patches.nEdge(1);
  my = patches.nEdge(2);
  x = x(1:mx:nx,:,:,:,:,:);
  y = y(:,1:my:ny,:,:,:,:);
  nx = nx/mx;
  ny = ny/my;
  u = reshape(u,mx,nx,my,ny,nVars,nEnsem,Nx,Ny);
  nVars = nVars*mx*my;
  u = reshape( permute(u,[2 4 1 3 5:8]) ...
             ,nx,ny,nVars,nEnsem,Nx,Ny);
end%if patches.nEdge
%{
\end{matlab}

The centre of each patch (as \verb|nx| and~\verb|ny| are
odd for centre-patch interpolation) is at indices
\begin{matlab}
%}
i0 = round((nx+1)/2);
j0 = round((ny+1)/2);
%{
\end{matlab}



\subsection{Periodic macroscale interpolation schemes}
\begin{matlab}
%}
if patches.periodic
%{
\end{matlab}
Get the size ratios of the patches.
\begin{matlab}
%}
rx = patches.ratio(1);
ry = patches.ratio(2);
%{
\end{matlab}



\subsubsection{Lagrange interpolation gives patch-edge values}
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
we naturally fill-in corner values. Start with
\(x\)-direction, and give most documentation for that case
as the \(y\)-direction is essentially the same.


\paragraph{\(x\)-normal edge values} The patch-edge values
are either interpolated from the next-to-edge values, or
from the centre-cross values (not the patch-centre value
itself as that seems to have worse properties in general). 
Have not yet implemented core averages.
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    U = u([2 nx-1],2:(ny-1),:,:,I,J);
  else % interpolate centre-cross values
    U = u(i0,2:(ny-1),:,:,I,J);
  end;%if patches.EdgyInt
%{
\end{matlab}
Just in case any last array dimension(s) are one, we force a
padding of the sizes, then adjoin the extra dimension for
the subsequent array of differences.
\begin{matlab}
%}
szUO=size(U); szUO=[szUO ones(1,6-length(szUO)) ordCC];
%{
\end{matlab}
Use finite difference formulas for the interpolation, so
store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in these arrays.  When parallel, in order
to preserve the distributed array structure we use an index
at the end for the differences.
\begin{matlab}
%}
  if ~patches.parallel, dmu = zeros(szUO); % 7D
  else   dmu = zeros(szUO,patches.codist); % 7D
  end%if patches.parallel
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
%    dmux(:,:,:,:,I,:,1) = (Ux(:,:,:,:,Ip,:)+Ux(:,:,:,:,Im,:))/2; % \mu
%    dmux(:,:,:,:,I,:,2) = (Ux(:,:,:,:,Ip,:)-Ux(:,:,:,:,Im,:)); % \delta
%    Ip = Ip(Ip); Im = Im(Im); % increase shifts to \pm2
%    dmuy(:,:,:,:,:,J,1) = (Ux(:,:,:,:,:,Jp)+Ux(:,:,:,:,:,Jm))/2; % \mu
%    dmuy(:,:,:,:,:,J,2) = (Ux(:,:,:,:,:,Jp)-Ux(:,:,:,:,:,Jm)); % \delta
%    Jp = Jp(Jp); Jm = Jm(Jm); % increase shifts to \pm2
  else %disp('starting standard interpolation')   
    dmu(:,:,:,:,I,:,1) = (U(:,:,:,:,Ip,:) ...
                         -U(:,:,:,:,Im,:))/2; %\mu\delta 
    dmu(:,:,:,:,I,:,2) = (U(:,:,:,:,Ip,:) ...
       -2*U(:,:,:,:,I,:) +U(:,:,:,:,Im,:));   %\delta^2    
  end% if patches.stag
%{
\end{matlab}
Recursively take $\delta^2$ of these to form successively
higher order centred differences in space.
\begin{matlab}
%}
  for k = 3:ordCC    
    dmu(:,:,:,:,I,:,k) =     dmu(:,:,:,:,Ip,:,k-2) ...
    -2*dmu(:,:,:,:,I,:,k-2) +dmu(:,:,:,:,Im,:,k-2);    
  end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet edge values for
each patch \cite[]{Roberts06d, Bunder2013b}, using weights
computed in \verb|configPatches2()|. Here interpolate to
specified order.

For the case where next-to-edge values interpolate to the
opposite edge-values: when we have an ensemble of
configurations, different configurations might be coupled to
each other, as specified by \verb|patches.le|,
\verb|patches.ri|, \verb|patches.to| and \verb|patches.bo|.
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two edges
u(nx,2:(ny-1),:,patches.ri,I,:) ...
  = U(1,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,1),-6).*dmu(1,:,:,:,:,:,:) ,7);  
u(1 ,2:(ny-1),:,patches.le,I,:,:) ...
  = U(k,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,1),-6).*dmu(k,:,:,:,:,:,:) ,7);
%{
\end{matlab}



\paragraph{\(y\)-normal edge values} Interpolate from either
the next-to-edge values, or the centre-cross-line values.
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    U = u(:,[2 ny-1],:,:,I,J);
  else % interpolate centre-cross values
    U = u(:,j0,:,:,I,J);
  end;%if patches.EdgyInt
%{
\end{matlab}
Adjoin extra dimension for the array of differences.
\begin{matlab}
%}
szUO=size(U); szUO=[szUO ones(1,6-length(szUO)) ordCC];
%{
\end{matlab}
Store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in this array.
\begin{matlab}
%}
  if ~patches.parallel, dmu = zeros(szUO); % 7D
  else   dmu = zeros(szUO,patches.codist); % 7D
  end%if patches.parallel
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    error('polynomial interpolation not yet for staggered patch coupling')
  else  %disp('starting standard interpolation')  
    dmu(:,:,:,:,:,J,1) = (U(:,:,:,:,:,Jp) ...
                         -U(:,:,:,:,:,Jm))/2; %\mu\delta 
    dmu(:,:,:,:,:,J,2) = (U(:,:,:,:,:,Jp) ...
       -2*U(:,:,:,:,:,J) +U(:,:,:,:,:,Jm));   %\delta^2
  end% if stag
%{
\end{matlab}
Recursively take $\delta^2$.
\begin{matlab}
%}
  for k = 3:ordCC    
    dmu(:,:,:,:,:,J,k) =     dmu(:,:,:,:,:,Jp,k-2) ...
    -2*dmu(:,:,:,:,:,J,k-2) +dmu(:,:,:,:,:,Jm,k-2);
  end
%{
\end{matlab}
Interpolate macro-values using the weights pre-computed by
\verb|configPatches2()|.  An ensemble of configurations may
have cross-coupling.
\begin{matlab}
%}
k = 1+patches.EdgyInt; % use centre or two edges
u(:,ny,:,patches.to,:,J) ...
  = U(:,1,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,2),-6).*dmu(:,1,:,:,:,:,:) ,7);
u(:,1 ,:,patches.bo,:,J) ...
  = U(:,k,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,2),-6).*dmu(:,k,:,:,:,:,:) ,7);
%{
\end{matlab}




\subsubsection{Case of spectral interpolation}
Assumes the domain is macro-periodic.
\begin{matlab}
%}
else% patches.ordCC<=0, spectral interpolation
%{
\end{matlab}
We interpolate in terms of the patch index, $j$~say, not
directly in space. As the macroscale fields are $N$-periodic
in the patch index~$I$, the macroscale Fourier transform
writes the centre-patch values as $U_I=\sum_{k}C_ke^{ik2\pi
I/N}$. Then the edge-patch values $U_{I\pm r}
=\sum_{k}C_ke^{ik2\pi/N(I\pm r)} =\sum_{k}C'_ke^{ik2\pi
I/N}$ where $C'_k=C_ke^{ikr2\pi/N}$. For $N$~patches we
resolve `wavenumbers' $|k|<N/2$, so set row vector
$\verb|ks|=k2\pi/N$ for `wavenumbers' $\mathcode`\,="213B
k=(0,1, \ldots, k_{\max}, -k_{\max}, \ldots, -1)$ for
odd~$N$, and $\mathcode`\,="213B k=(0,1, \ldots, k_{\max},
\pm(k_{\max}+1) -k_{\max}, \ldots, -1)$ for even~$N$.

Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches2|
tests there are an even number of patches). Then the
patch-ratio is effectively halved. The patch edges are near
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
Interpolate the two directions in succession, in this way we
naturally fill-in edge-corner values. Start with
\(x\)-direction, and give most documentation for that case
as the other is essentially the same. Need these indices of
patch interior.
\begin{matlab}
%}
ix = 2:nx-1;   iy = 2:ny-1;  
%{
\end{matlab}



\paragraph{\(x\)-normal edge values} Now set wavenumbers
into a vector at the correct dimension.  In the case of
even~$N$ these compute the $+$-case for the highest
wavenumber zig-zag mode, $\mathcode`\,="213B k=(0,1, \ldots,
k_{\max}, +(k_{\max}+1) -k_{\max}, \ldots, -1)$.
\begin{matlab}
%}
  kMax = floor((Nx-1)/2); 
  kr = shiftdim( rx*2*pi/Nx*(mod((0:Nx-1)+kMax,Nx)-kMax) ,-3);
%{
\end{matlab}
Compute the Fourier transform of the centre-cross values.
Unless doing patch-edgy interpolation when FT the
next-to-edge values.  If there are an even number of points,
then if complex, treat as positive wavenumber, but if real,
treat as cosine. When using an ensemble of configurations,
different configurations might be coupled to each other, as
specified by \verb|patches.le|, \verb|patches.ri|,
\verb|patches.to| and \verb|patches.bo|.
\begin{matlab}
%}
if ~patches.EdgyInt
     Cm = fft( u(i0,iy,:,:,:,:) ,[],5); 
     Cp = Cm;
else 
     Cm = fft( u(   2,iy ,:,patches.le,:,:) ,[],5);
     Cp = fft( u(nx-1,iy ,:,patches.ri,:,:) ,[],5);
end%if ~patches.EdgyInt  
%{
\end{matlab}
Now invert the Fourier transforms to complete interpolation.
Enforce reality when appropriate. 
\begin{matlab}
%}
u(nx,iy,:,:,:,:) = uclean( ifft( ...
    Cm.*exp(1i*(stagShift+kr))  ,[],5) );
u( 1,iy,:,:,:,:) = uclean( ifft( ...
    Cp.*exp(1i*(stagShift-kr))  ,[],5) );
%{
\end{matlab}



\paragraph{\(y\)-normal edge values} Set wavenumbers into a
vector.
\begin{matlab}
%}
  kMax = floor((Ny-1)/2); 
  kr = shiftdim( ry*2*pi/Ny*(mod((0:Ny-1)+kMax,Ny)-kMax) ,-4);
%{
\end{matlab}
Compute the Fourier transform of the patch values on the
centre-lines for all the fields.  
\begin{matlab}
%}
if ~patches.EdgyInt
     Cm = fft( u(:,j0,:,:,:,:) ,[],6); 
     Cp = Cm;
else 
     Cm = fft( u(:,2    ,:,patches.bo,:,:) ,[],6);
     Cp = fft( u(:,ny-1 ,:,patches.to,:,:) ,[],6);
end%if ~patches.EdgyInt  
%{
\end{matlab}
Invert the Fourier transforms to complete interpolation. 
\begin{matlab}
%}
u(:,ny,:,:,:,:) = uclean( ifft( ...
    Cm.*exp(1i*(stagShift+kr))  ,[],6) );
u(:, 1,:,:,:,:) = uclean( ifft( ...
    Cp.*exp(1i*(stagShift-kr))  ,[],6) );
%{
\end{matlab}

\begin{matlab}
%}
end% if ordCC>0 else, so spectral 
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
Determine the order of interpolation~\verb|px| and~\verb|py|
(potentially different in the different directions!), and
hence size of the (forward) divided difference tables
in~\verb|F|~(7D) for interpolating to left/right, and
top/bottom edges. Because of the product-form of the patch
grid, and because we are doing \emph{only} either edgy
interpolation or cross-patch interpolation (\emph{not} just
the centre patch value), the interpolations are all 1D
interpolations.
\begin{matlab}
%}
if patches.ordCC<1
     px = Nx-1; py = Ny-1;
else px = min(patches.ordCC,Nx-1); 
     py = min(patches.ordCC,Ny-1); 
end
ix=2:nx-1;  iy=2:ny-1; % indices of edge 'interior' (ix n/a)
%{
\end{matlab}

\subsubsection{\(x\)-direction values}
Set function values in first `column' of the tables for
every variable and across ensemble.  For~\verb|EdgyInt|, the
`reversal' of the next-to-edge values are because their
values are to interpolate to the opposite edge of each
patch. \todo{Have no plans to implement core averaging as
yet.}
\begin{matlab}
%}
  F = nan(patches.EdgyInt+1,ny-2,nVars,nEnsem,Nx,Ny,px+1);
  if patches.EdgyInt % interpolate next-to-edge values
    F(:,:,:,:,:,:,1) = u([nx-1 2],iy,:,:,:,:);
    X = x([nx-1 2],:,:,:,:,:);
  else % interpolate mid-patch cross-patch values 
    F(:,:,:,:,:,:,1) = u(i0,iy,:,:,:,:);
    X = x(i0,:,:,:,:,:);
  end%if patches.EdgyInt
%{
\end{matlab}

\paragraph{Form tables of divided differences} Compute
tables of (forward) divided differences
\cite[e.g.,][]{DividedDifferences} for every variable, and
across ensemble, and for left/right edges. Recursively find
all divided differences.
\begin{matlab}
%}
for q = 1:px
  i = 1:Nx-q;
  F(:,:,:,:,i,:,q+1) ...
  = (F(:,:,:,:,i+1 ,:,q)-F(:,:,:,:,i,:,q)) ...
  ./(X(:,:,:,:,i+q,:)   -X(:,:,:,:,i,:));
end
%{
\end{matlab}

\paragraph{Interpolate with divided differences} Now
interpolate to find the edge-values on left/right edges
at~\verb|Xedge| for every interior~\verb|Y|.
\begin{matlab}
%}
Xedge = x([1 nx],:,:,:,:,:);
%{
\end{matlab}
Code Horner's recursive evaluation of the interpolation
polynomials.  Indices~\verb|i| are those of the left edge of
each interpolation stencil, because the table is of forward
differences.  This alternative: the case of order~\(p_x\)
and~\(p_y\) interpolation across the domain, asymmetric near
the boundaries of the rectangular domain.
\begin{matlab}
%}
  i = max(1,min(1:Nx,Nx-ceil(px/2))-floor(px/2));
  Uedge = F(:,:,:,:,i,:,px+1);
  for q = px:-1:1
    Uedge = F(:,:,:,:,i,:,q)+(Xedge-X(:,:,:,:,i+q-1,:)).*Uedge;
  end
%{
\end{matlab}

Finally, insert edge values into the array of field values,
using the required ensemble shifts.  
\begin{matlab}
%}
u(1 ,iy,:,patches.le,:,:) = Uedge(1,:,:,:,:,:);
u(nx,iy,:,patches.ri,:,:) = Uedge(2,:,:,:,:,:);
%{
\end{matlab}

\subsubsection{\(y\)-direction values}
Set function values in first `column' of the tables for
every variable and across ensemble.
\begin{matlab}
%}
  F = nan(nx,patches.EdgyInt+1,nVars,nEnsem,Nx,Ny,py+1);
  if patches.EdgyInt % interpolate next-to-edge values
    F(:,:,:,:,:,:,1) = u(:,[ny-1 2],:,:,:,:);
    Y = y(:,[ny-1 2],:,:,:,:);
  else % interpolate mid-patch cross-patch values 
    F(:,:,:,:,:,:,1) = u(:,j0,:,:,:,:);
    Y = y(:,j0,:,:,:,:);
  end;
%{
\end{matlab}
Form tables of divided differences.
\begin{matlab}
%}
for q = 1:py
  j = 1:Ny-q;
  F(:,:,:,:,:,j,q+1) ...
  = (F(:,:,:,:,:,j+1 ,q)-F(:,:,:,:,:,j,q)) ...
  ./(Y(:,:,:,:,:,j+q)   -Y(:,:,:,:,:,j));
end
%{
\end{matlab}
Interpolate to find the edge-values on top/bottom
edges~\verb|Yedge| for every~\(x\).
\begin{matlab}
%}
Yedge = y(:,[1 ny],:,:,:,:);
%{
\end{matlab}
Code Horner's recursive evaluation of the interpolation
polynomials.  Indices~\verb|j| are those of the bottom edge
of each interpolation stencil, because the table is of
forward differences.
\begin{matlab}
%}
  j = max(1,min(1:Ny,Ny-ceil(py/2))-floor(py/2));
  Uedge = F(:,:,:,:,:,j,py+1);
  for q = py:-1:1
    Uedge = F(:,:,:,:,:,j,q)+(Yedge-Y(:,:,:,:,:,j+q-1)).*Uedge;
  end
%{
\end{matlab}
Finally, insert edge values into the array of field values,
using the required ensemble shifts.  
\begin{matlab}
%}
u(:,1 ,:,patches.bo,:,:) = Uedge(:,1,:,:,:,:);
u(:,ny,:,patches.to,:,:) = Uedge(:,2,:,:,:,:);
%{
\end{matlab}


\subsubsection{Optional NaNs for safety}
We want a user to set outer edge values on the extreme
patches according to the microscale boundary conditions that
hold at the extremes of the domain.  Consequently, unless testing, override
their computed interpolation values with~\verb|NaN|.
\begin{matlab}
%}
if isfield(patches,'intTest')&&patches.intTest
else % usual case
    u( 1,:,:,:, 1,:) = nan;
    u(nx,:,:,:,Nx,:) = nan;
    u(:, 1,:,:,:, 1) = nan;
    u(:,ny,:,:,:,Ny) = nan;
end%if
%{
\end{matlab}

End of the non-periodic interpolation code.
\begin{matlab}
%}
end%if patches.periodic else
%{
\end{matlab}

\paragraph{Unfold multiple edges}  No need to restore~\(x,y\).
\begin{matlab}
%}
if mean(patches.nEdge)>1
  nVars = nVars/(mx*my);
  u = reshape( u ,nx,ny,mx,my,nVars,nEnsem,Nx,Ny);
  nx = nx*mx;
  ny = ny*my;
  u = reshape( permute(u,[3 1 4 2 5:8]) ...
             ,nx,ny,nVars,nEnsem,Nx,Ny);
end%if patches.nEdge
%{
\end{matlab}

Fin, returning the 6D array of field values with
interpolated edges. 
\begin{matlab}
%}
end% function patchEdgeInt2
%{
\end{matlab}
\end{devMan} 
%}
