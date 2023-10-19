% patchEdgeInt3() provides the interpolation across 3D space
% for 3D patches of simulations of a smooth lattice system
% such as PDE discretisations.  AJR, Aug 2020 -- 19 Oct 2023
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



\subsection{Interpolate over the three successive directions}
Only use the interior values of the fields for interpolating
to the edges.
\begin{matlab}
%}
u = u(2:nx-1,2:ny-1,2:nz-1,:,:,:,:,:); 
%{
\end{matlab}
Interpolate in turn, the edge or mid-patch faces normal to
the \(x,y,z\)-directions, in this way we naturally fill-in
face-edge and corner values.
\begin{matlab}
%}
u = patchEdgeIntCore(1,u,x,patches,stagShift ...
    ,1,nx,(ny-2)*(nz-2)*nVars,nEnsem,1,Nx,Ny*Nz ...
    ,patches.le,patches.ri);
u = patchEdgeIntCore(2,u,y,patches,stagShift ...
    ,nx,ny,(nz-2)*nVars,nEnsem,Nx,Ny,Nz ...
    ,patches.bo,patches.to);
u = patchEdgeIntCore(3,u,z,patches,stagShift ...
    ,nx*ny,nz,    nVars,nEnsem,Nx*Ny,Nz,1 ...
    ,patches.ba,patches.fr);
%{
\end{matlab}
Restore array~\verb|u| to its original shape.
\begin{matlab}
%}
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
