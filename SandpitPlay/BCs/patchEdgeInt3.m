% patchEdgeInt3() provides the interpolation across 3D space
% for 3D patches of simulations of a smooth lattice system
% such as PDE discretisations.  AJR, Aug 2020 -- 17 Oct 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt3()}: sets 3D patch
face values from 3D macroscale interpolation}
\label{sec:patchEdgeInt3}


Couples 3D patches across 3D space by computing their face
values via macroscale interpolation.  Assumes patch face
values are determined by macroscale interpolation of the
patch centre-plane values \cite[]{Roberts2011a,
Bunder2019d}, or patch next-to-face values which appears
better \cite[]{Bunder2020a}.  This function is primarily
used by \verb|patchSys3()| but is also useful for user
graphics.\footnote{Script \texttt{patchEdgeInt3test.m}
verifies most of this code.}

Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|), or otherwise via the global
struct~\verb|patches|.
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
  ,'patchEdgeInt3: input u has wrong size for parameters')
u = reshape(u,[nx ny nz nVars nEnsem Nx Ny Nz]);
%{
\end{matlab}


\paragraph{Implement multiple width edges by folding}
Subsample~\(x,y,z\) coordinates, noting it is only
differences that count \emph{and} the microgrid~\(x,y,z\)
spacing must be uniform.
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


\paragraph{Staggered grid}
Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches3|
tests there are an even number of patches). Then the
patch-ratio is effectively halved. The patch faces are near
the middle of the gaps and swapped.
\begin{matlab}
%}
 if patches.stag % transform by doubling the number of fields
 error('staggered grid not yet implemented????')
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



Only use the interior values of the fields for interpolating
to the edges.
\begin{matlab}
%}
u = u(2:nx-1,2:ny-1,2:nz-1,:,:,:,:,:); 
%{
\end{matlab}



\subsection{Loop over the three successive directions}
Interpolate in turn, the edge or mid-patch faces normal to
the \(x,y,z\)-directions, in this way we naturally fill-in
face-edge and corner values.
\begin{matlab}
%}
for l=1:3
%{
\end{matlab}
\paragraph{Organise the field into a generic shape} 
For the interpolation use with seven arguments in order of
stuff, l-microscale, stuff, ensemble, stuff, l-macroscale,
stuff \todo{Check whether this wrecks parallel code??}: set
dimensions; set the centre index~\verb|l0| of each patch (as
\verb|nx|, \verb|ny| and \verb|nz| are odd for centre-patch
interpolation).
\begin{matlab}
%}
switch l
case 1 % x-direction
    n1=1; nl=nx; n2=(ny-2)*(nz-2)*nVars;
    N1=1; Nl=Nx; N2=Ny*Nz;
    l0 = round((nx+1)/2);
    lLo = patches.le;  lHi = patches.ri; 
    q = reshape(x,1,nl,1,1,1,Nl,1);
case 2 % y-direction
    n1=nx; nl=ny; n2=(nz-2)*nVars;
    N1=Nx; Nl=Ny; N2=Nz;
    l0 = round((ny+1)/2);
    lLo = patches.bo;  lHi = patches.to; 
    q = reshape(y,1,nl,1,1,1,Nl,1);
case 3 % z-direction
    n1=nx*ny; nl=nz; n2=nVars;
    N1=Nx*Ny; Nl=Nz; N2=1;
    l0 = round((nz+1)/2);
    lLo = patches.ba;  lHi = patches.fr; 
    q = reshape(z,1,nl,1,1,1,Nl,1);
end%switch l
%{
\end{matlab}
Reshape the field values accordingly, fattened for the new
edge values, and set macroscale indices and periodic
neighbours.
\begin{matlab}
%}
u = cat(2,nan(n1,  1 ,n2,nEnsem,N1,Nl,N2) ...
   ,reshape(u,n1,nl-2,n2,nEnsem,N1,Nl,N2) ...
         ,nan(n1,  1 ,n2,nEnsem,N1,Nl,N2) );
Il=1:Nl; Ip=mod(Il,Nl)+1; Im=mod(Il-2,Nl)+1;
%{
\end{matlab}




\subsection{Periodic macroscale interpolation schemes}
First get the size ratios of the patches in the direction.
\begin{matlab}
%}
if patches.periodic(l)
rl = patches.ratio(l);
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


\paragraph{\(l\)-normal face values} The patch-edge values
are either interpolated from the next-to-edge-face values,
or from the centre-cross-plane values (not the patch-centre
value itself as that seems to have worse properties in
general).  Have not yet implemented core averages.
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-face values    
    U = u(:,[2 nl-1],:,:,:,:,:);
  else % interpolate centre-cross values
    U = u(:,   l0   ,:,:,:,:,:);
  end;%if patches.EdgyInt
%{
\end{matlab}
Just in case any last array dimension(s) are one, we force a
padding of the sizes, then adjoin the extra dimension for
the subseluent array of differences.
\begin{matlab}
%}
szUO=size(U); szUO=[szUO ones(1,7-length(szUO)) ordCC];
%{
\end{matlab}
Use finite difference formulas for the interpolation, so
store finite differences ($\mu\delta, \delta^2, \mu\delta^3,
\delta^4, \ldots$) in these arrays.  When parallel, in order
to preserve the distributed array structure we use an index
at the end for the differences.
\begin{matlab}
%}
  if ~patches.parallel, dmu = zeros(szUO); % 8D
  else   dmu = zeros(szUO,patches.codist); % 8D
  end%if patches.parallel
%{
\end{matlab}
First compute differences $\mu\delta$ and $\delta^2$.
\begin{matlab}
%}
  %disp('starting standard interpolation')  
  dmu(:,:,:,:,:,Il,:,1) = (U(:,:,:,:,:,Ip,:) ...
                          -U(:,:,:,:,:,Im,:))/2; %\mu\delta 
  dmu(:,:,:,:,:,Il,:,2) = (U(:,:,:,:,:,Ip,:) ...
     -2*U(:,:,:,:,:,Il,:) +U(:,:,:,:,:,Im,:));   %\delta^2    
%{
\end{matlab}
Recursively take $\delta^2$ of these to form successively
higher order centred differences in space.
\begin{matlab}
%}
  for k = 3:ordCC    
    dmu(:,:,:,:,:,Il,:,k) =     dmu(:,:,:,:,:,Ip,:,k-2) ...
    -2*dmu(:,:,:,:,:,Il,:,k-2) +dmu(:,:,:,:,:,Im,:,k-2);    
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
each other, as specified in \verb|lLo| and~\verb|lHi| (from
\verb|patches.le|, \verb|patches.ri|, \verb|patches.bo|,
\verb|patches.to|, \verb|patches.ba|, and
\verb|patches.fr|).
\begin{matlab}
%}
k=1+patches.EdgyInt; % use centre or two faces
u(:,nl,:,lHi,:,:,:) ...
  = U(:,1,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsr(:,l),-7).*dmu(:,1,:,:,:,:,:,:) ,8);  
u(:,1 ,:,lLo,:,:,:) ...
  = U(:,k,:,:,:,:,:)*(1-patches.stag) ...
  +sum( shiftdim(patches.Cwtsl(:,l),-7).*dmu(:,k,:,:,:,:,:,:) ,8);
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


\paragraph{\(l\)-normal face values} Now set wavenumbers
into a vector at the correct dimension.  In the case of
even~$N$ these compute the $+$-case for the highest
wavenumber zig-zag mode, $\mathcode`\,="213B k=(0,1, \ldots,
k_{\max}, +(k_{\max}+1) -k_{\max}, \ldots, -1)$.
\begin{matlab}
%}
  kMax = floor((Nl-1)/2); 
  kr = shiftdim( rl*2*pi/Nl*(mod((0:Nl-1)+kMax,Nl)-kMax) ,-4);
%{
\end{matlab}
Compute the Fourier transform of the patch values on the
centre-planes for all the fields.  Unless doing patch-edgy
interpolation when FT the next-to-face values.  If there are
an even number of points, then if complex, treat as positive
wavenumber, but if real, treat as cosine. When using an
ensemble of configurations, different configurations might
be coupled to each other, as specified by \verb|lLo|
and~\verb|lHi| (from \verb|patches.le|, \verb|patches.ri|,
\verb|patches.to|, \verb|patches.bo|, \verb|patches.fr| and
\verb|patches.ba|).
\begin{matlab}
%}
if ~patches.EdgyInt
     Cm = fft( u(:,l0,:,:,:,:,:) ,[],6); 
     Cp = Cm;
else 
     Cm = fft( u(:,   2,:,lLo,:,:,:) ,[],6);
     Cp = fft( u(:,nl-1,:,lHi,:,:,:) ,[],6);
end%if ~patches.EdgyInt  
%{
\end{matlab}
Now invert the Fourier transforms to complete interpolation.
Enforce reality when appropriate with \verb|uclean()|. 
\begin{matlab}
%}
u(:,nl,:,:,:,:,:) = uclean( ifft( ...
    Cm.*exp(1i*(stagShift+kr))  ,[],6) );
u(:, 1,:,:,:,:,:) = uclean( ifft( ...
    Cp.*exp(1i*(stagShift-kr))  ,[],6) );
%{
\end{matlab}

End the two periodic cases of spectral or finite-difference
interpolation.
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
'not implementing staggered grids for non-periodic')
%{
\end{matlab}
Determine the order of interpolation~\verb|pl| (potentially
different in the different directions!), and hence size of
the (forward) divided difference tables in~\verb|F|~(8D) for
interpolating to left/right, top/bottom, and front/back
faces. Because of the product-form of the patch grid, and
because we are doing \emph{only} either edgy interpolation
or cross-patch interpolation (\emph{not} just the centre
patch value), the interpolations are all 1D interpolations.
\begin{matlab}
%}
if patches.ordCC<1,  pl = Nl-1;
  else pl = min(patches.ordCC,Nl-1); 
  end
%{
\end{matlab}

\subsubsection{\(l\)-direction values}
Set function values in second column of the tables for every
variable and across ensemble.  For~\verb|EdgyInt|, the
`reversal' of the next-to-face values are because their
values are to interpolate to the opposite face of each
patch. \todo{Have no plans to implement core averaging as
yet.}
\begin{matlab}
%}
  F = nan(n1,patches.EdgyInt+1,n2,nEnsem,N1,Nl,N2,pl+1);
  if patches.EdgyInt % interpolate next-to-face values
    F(:,:,:,:,:,:,:,1) = u(:,[nl-1 2],:,:,:,:,:);
    Q = q(:,[nl-1 2],:,:,:,:,:);
  else % interpolate mid-patch cross-patch values 
    F(:,:,:,:,:,:,:,1) = u(:,   l0   ,:,:,:,:,:);
    Q = q(:,   l0   ,:,:,:,:,:);
  end%if patches.EdgyInt
%{
\end{matlab}

\paragraph{Form tables of divided differences} 
Compute tables of (forward) divided differences
\cite[e.g.,][]{DividedDifferences} for every variable, and
across ensemble, and in both directions, and for all three
types of faces (left/right, top/bottom, and front/back).
Recursively find all divided differences in the respective
direction.
\begin{matlab}
%}
for m = 1:pl
  i = 1:Nl-m;
  F(:,:,:,:,:,i,:,m+1) ...
  = ( F(:,:,:,:,:,i+1,:,m)-F(:,:,:,:,:,i,:,m)) ...
   ./(Q(:,:,:,:,:,i+m,:)  -Q(:,:,:,:,:,i,:));
end
%{
\end{matlab}

\paragraph{Interpolate with divided differences} Now
interpolate to find the face-values on left/right faces
at~\verb|Qface| for every interior in other coordinates.
\begin{matlab}
%}
Qface = q(:,[1 nl],:,:,:,:,:);
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
  i = max(1,min(1:Nl,Nl-ceil(pl/2))-floor(pl/2));
  Uface = F(:,:,:,:,:,i,:,pl+1);
  for m = pl:-1:1
    Uface = F(:,:,:,:,:,i,:,m) ...
    +(Qface-Q(:,:,:,:,:,i+m-1,:)).*Uface;
  end
%{
\end{matlab}

Finally, insert face values into the array of field values,
using the required ensemble shifts.  
\begin{matlab}
%}
u(:,1 ,:,lLo,:,:,:) = Uface(:,1,:,:,:,:,:);
u(:,nl,:,lHi,:,:,:) = Uface(:,2,:,:,:,:,:);
%{
\end{matlab}



\subsubsection{Optional NaNs for safety}
We want a user to set outer face values on the extreme
patches according to the microscale boundary conditions that
hold at the extremes of the domain.  Consequently, unless
testing, override their computed interpolation values
with~\verb|NaN|.
\begin{matlab}
%}
if isfield(patches,'intTest')&&patches.intTest
else %disp('usual case requires user to set bdry values')
    u(:, 1,:,:,:, 1,:) = nan;
    u(:,nl,:,:,:,Nl,:) = nan;
end%if
%{
\end{matlab}

End of the non-periodic interpolation code.
\begin{matlab}
%}
end%if patches.periodic else
%{
\end{matlab}

End the loop over the three spatial directions, and so
restore array~\verb|u| to its original shape.
\begin{matlab}
%}
end%for q 
u = reshape(u,nx,ny,nz,nVars,nEnsem,Nx,Ny,Nz);
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
end% function patchEdgeInt3
%{
\end{matlab}
\end{devMan} 
%}
