% patchEdgeInt1() provides the interpolation across 1D space
% for 1D patches of simulations of a lattice system such as
% PDE discretisations.  AJR & JB, Sep 2018 -- Dec 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt1()}: sets patch-edge values
from interpolation over the 1D macroscale}
\label{sec:patchEdgeInt1}


Couples 1D patches across 1D space by computing their edge
values from macroscale interpolation of either the mid-patch
value \cite[]{Roberts00a, Roberts06d}, or the patch-core
average \cite[]{Bunder2013b}, or the opposite next-to-edge
values \cite[]{Bunder2020a}---this last alternative often
maintains symmetry. This function is primarily used by
\verb|patchSmooth1()| but is also useful for user graphics.
When using core averages, assumes the averages are sensible
macroscale variables: then patch edge values are determined
by macroscale interpolation of the core averages
\citep{Bunder2013b}. 
\footnote{Script \texttt{patchEdgeInt1test.m} verifies this code.}

Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|), or otherwise via the global
struct~\verb|patches|.
\begin{matlab}
%}
function u=patchEdgeInt1(u,patches)
if nargin<2, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector\slash array of length
$\verb|nSubP| \cdot \verb|nVars|\cdot \verb|nEnsem|\cdot
\verb|nPatch|$ where there are $\verb|nVars|\cdot
\verb|nEnsem|$ field values at each of the points in the
$\verb|nSubP| \times \verb|nPatch|$ multiscale spatial grid.

\item \verb|patches| a struct largely set by
\verb|configPatches1()|, and which includes the following.
\begin{itemize}

\item \verb|.x| is $\verb|nSubP| \times1 \times1 \times
\verb|nPatch|$ array of the spatial locations~$x_{iI}$ of
the microscale grid points in every patch. Currently it
\emph{must} be an equi-spaced lattice on both macro- and
microscales.

\item \verb|.ordCC| is order of interpolation, integer~$\geq
-1$.

\item \verb|.stag| in $\{0,1\}$ is one for staggered grid
(alternating) interpolation, and zero for ordinary grid.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the coupling
coefficients for finite width interpolation.

\item \verb|.EdgyInt|, true/false, is true for interpolating
patch-edge values from opposite next-to-edge values (often
preserves symmetry).

\item \verb|.nEnsem| the number of realisations in the ensemble.

\item \verb|.parallel| whether serial or parallel.

\end{itemize}
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 4D array, $\verb|nSubP| \times
\verb|nVars| \times \verb|nEnsem| \times \verb|nPatch|$, of
the fields with edge values set by interpolation.
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
[nx,~,~,Nx] = size(patches.x);
nEnsem = patches.nEnsem;
nVars = round(numel(u)/numel(patches.x)/nEnsem);
assert(numel(u) == nx*nVars*nEnsem*Nx ...
  ,'patchEdgeInt1: input u has wrong size for parameters')
u = reshape(u,nx,nVars,nEnsem,Nx);
%{
\end{matlab}
If the user has not defined the patch core, then we assume
it to be a single point in the middle of the patch, unless
we are interpolating from next-to-edge values. 
Get the size ratios of the patches.
\begin{matlab}
%}
r = patches.ratio(1);
%{
\end{matlab}

For the moment assume the physical domain is macroscale
periodic so that the coupling formulas are simplest. Should
eventually cater for periodic, odd-mid-gap, even-mid-gap,
even-mid-patch, Dirichlet, Neumann, ??. These index vectors
point to patches and their two immediate neighbours.
\begin{matlab}
%}
I = 1:Nx; Ip = mod(I,Nx)+1; Im = mod(I-2,Nx)+1;
%{
\end{matlab}
Calculate centre of each patch and the surrounding core
(\verb|nx| and \verb|nCore| are both odd).
\begin{matlab}
%}
i0 = round((nx+1)/2);
c = round((patches.nCore-1)/2);
%{
\end{matlab}



\paragraph{Lagrange interpolation gives patch-edge values}
Consequently, compute centred differences of the patch
core/edge averages/values for the macro-interpolation of all
fields. Assumes the domain is macro-periodic. 
\begin{matlab}
%}
if patches.ordCC>0 % then finite-width polynomial interpolation
%{
\end{matlab}
\begin{matlab}
%}
  if patches.EdgyInt % interpolate next-to-edge values
    Ux = u([2 nx-1],:,:,I);
  else % interpolate mid-patch values/sums
    Ux = sum( u((i0-c):(i0+c),:,:,I) ,1);
  end;
%{
\end{matlab}
Just in case any last array dimension(s) are one, we have to
force a padding of the sizes, then adjoin the extra
dimension for the subsequent array of differences.
\begin{matlab}
%}
szUxO=size(Ux); 
szUxO=[szUxO ones(1,4-length(szUxO)) patches.ordCC];
%{
\end{matlab}
Use finite difference formulas for the interpolation, so
store finite differences in these arrays.  When parallel, in
order to preserve the distributed array structure we use an
index at the end for the differences.
\begin{matlab}
%}
  if patches.parallel
    dmu = zeros(szUxO,patches.codist); % 5D
  else
    dmu = zeros(szUxO); % 5D
  end
%{
\end{matlab}
First compute differences, either $\mu$ and $\delta$, or
$\mu\delta$ and $\delta^2$ in space.
\begin{matlab}
%}
  if patches.stag % use only odd numbered neighbours
    dmu(:,:,:,I,1) = (Ux(:,:,:,Ip)+Ux(:,:,:,Im))/2; % \mu
    dmu(:,:,:,I,2) = (Ux(:,:,:,Ip)-Ux(:,:,:,Im)); % \delta
    Ip = Ip(Ip); Im = Im(Im); % increase shifts to \pm2
  else % standard
    dmu(:,:,:,I,1) = (Ux(:,:,:,Ip)-Ux(:,:,:,Im))/2; % \mu\delta
    dmu(:,:,:,I,2) = (Ux(:,:,:,Ip)-2*Ux(:,:,:,I) ...
                     +Ux(:,:,:,Im)); % \delta^2
  end% if stag
%{
\end{matlab}
Recursively take $\delta^2$ of these to form successively
higher order centred differences in space.
\begin{matlab}
%}
  for k = 3:patches.ordCC
    dmu(:,:,:,:,k) =       dmu(:,:,:,Ip,k-2) ...
      -2*dmu(:,:,:,I,k-2) +dmu(:,:,:,Im,k-2);
  end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet edge values for
each patch \cite[]{Roberts06d, Bunder2013b}, using weights
computed in \verb|configPatches1()|. Here interpolate to
specified order.

For the case where single-point values interpolate to
patch-edge values: when we have an ensemble of
configurations, different realisations are coupled to each
other as specified by \verb|patches.le| and
\verb|patches.ri|.
\begin{matlab}
%}
  if patches.nCore==1
    k=1+patches.EdgyInt; % use centre/core or two edges
    u(nx,:,patches.ri,I) = Ux(1,:,:,:)*(1-patches.stag) ...
      +sum( shiftdim(patches.Cwtsr,-4).*dmu(1,:,:,:,:) ,5);
    u(1 ,:,patches.le,I) = Ux(k,:,:,:)*(1-patches.stag) ...      
      +sum( shiftdim(patches.Cwtsl,-4).*dmu(k,:,:,:,:) ,5);
%{
\end{matlab}
For a non-trivial core then more needs doing: the core (one
or more) of each patch interpolates to the edge action
regions. When more than one in the core, the edge is set
depending upon near edge values so the average near the edge
is correct.
\begin{matlab}
%}
  else error('not yet considered, july--dec 2020 ??')
    u(nx,:,:,I) = Ux(:,:,I)*(1-patches.stag) ...
      + reshape(-sum(u((nx-patches.nCore+1):(nx-1),:,:,I),1) ...
      + sum( patches.Cwtsr.*dmu ),Nx,nVars);
    u(1,:,:,I) = Ux(:,:,I)*(1-patches.stag) ...      
      + reshape(-sum(u(2:patches.nCore,:,:,I),1)  ...
      + sum( patches.Cwtsl.*dmu ),Nx,nVars);
  end;
%{
\end{matlab}



\paragraph{Case of spectral interpolation}
Assumes the domain is macro-periodic. 
\begin{matlab}
%}
else% spectral interpolation
%{
\end{matlab}
As the macroscale fields are $N$-periodic, the macroscale
Fourier transform writes the centre-patch values as $U_j =
\sum_{k}C_ke^{ik2\pi j/N}$. Then the edge-patch values
$U_{j\pm r} =\sum_{k}C_ke^{ik2\pi/N(j\pm r)}
=\sum_{k}C'_ke^{ik2\pi j/N}$ where $C'_k =
C_ke^{ikr2\pi/N}$. For \verb|Nx|~patches we resolve
`wavenumbers' $|k|<\verb|Nx|/2$, so set row vector
$\verb|ks| = k2\pi/N$ for `wavenumbers'
$\mathcode`\,="213B k = (0,1, \ldots, k_{\max}, -k_{\max},
\ldots, -1)$ for odd~$N$, and $\mathcode`\,="213B k =
(0,1, \ldots, k_{\max}, (k_{\max}+1), -k_{\max}, \ldots,
-1)$ for even~$N$.

Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches1()|
tests that there are an even number of patches). Then the
patch-ratio is effectively halved. The patch edges are near
the middle of the gaps and swapped. Have not yet tested
whether works for Edgy Interpolation??
\begin{matlab}
%}
  if patches.stag % transform by doubling the number of fields
    v = nan(size(u)); % currently to restore the shape of u
    u = [u(:,:,:,1:2:Nx) u(:,:,:,2:2:Nx)];
    stagShift = 0.5*[ones(1,nVars) -ones(1,nVars)];
    iV = [nVars+1:2*nVars 1:nVars]; % scatter interp to alternate field
    r = r/2;           % ratio effectively halved
    Nx = Nx/2; % halve the number of patches
    nVars = nVars*2;   % double the number of fields
  else % the values for standard spectral
    stagShift = 0;  
    iV = 1:nVars;
  end
%{
\end{matlab}
Now set wavenumbers (when \verb|Nx| is even then highest
wavenumber is~$\pi$).
\begin{matlab}
%}
  kMax = floor((Nx-1)/2); 
  ks = shiftdim( ...
      2*pi/Nx*(mod((0:Nx-1)+kMax,Nx)-kMax) ...
      ,-2); 
%{
\end{matlab}

Compute the Fourier transform across patches of the patch
centre or next-to-edge values for all the fields. If there
are an even number of points, then if complex, treat as
positive wavenumber, but if real, treat as cosine. When
using an ensemble of configurations, different
configurations might be coupled to each other, as specified
by \verb|patches.le| and \verb|patches.ri|.
\begin{matlab}
%}
if ~patches.EdgyInt
    Cleft = fft(u(i0  ,:,:,:),[],4); 
    Cright = Cleft;
else
    Cleft = fft(u(2   ,:,:,:),[],4);
    Cright= fft(u(nx-1,:,:,:),[],4);
end
%{
\end{matlab}
The inverse Fourier transform gives the edge values via a
shift a fraction~$r$ to the next macroscale grid point.
\begin{matlab}
%}
  u(nx,iV,patches.ri,:) = uclean( ifft( ...
      Cleft.*exp(1i*ks.*(stagShift+r)) ,[],4));
  u(1 ,iV,patches.le,:) = uclean( ifft( ...
      Cright.*exp(1i*ks.*(stagShift-r)) ,[],4));
%{
\end{matlab}
Restore staggered grid when appropriate. This dimensional
shifting appears to work.
Is there a better way to do this??
\begin{matlab}
%}
if patches.stag
  nVars = nVars/2;  
  u=reshape(u,nx,nVars,2,nEnsem,Nx);
  Nx = 2*Nx;
  v(:,:,:,1:2:Nx) = u(:,:,1,:,:);
  v(:,:,:,2:2:Nx) = u(:,:,2,:,:);    
  u = v;
end
end% if spectral 
%{
\end{matlab}
Fin, returning the 4D array of field values.  
\end{devMan}
%}