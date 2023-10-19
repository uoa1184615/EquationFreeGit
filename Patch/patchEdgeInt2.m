% patchEdgeInt2() provides the interpolation across 2D space
% for 2D patches of simulations of a smooth lattice system
% such as PDE discretisations.  AJR, Aug 2020 -- 19 Oct 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchEdgeInt2()}: sets 2D patch
face values from 2D macroscale interpolation}
\label{sec:patchEdgeInt2}


Couples 2D patches across 2D space by computing their edge
values via macroscale interpolation.  Assumes patch edge
values are determined by macroscale interpolation of the
patch centre-plane values \cite[]{Roberts2011a,
Bunder2019d}, or patch next-to-edge values which appears
better \cite[]{Bunder2020a}.  This function is primarily
used by \verb|patchSys3()| but is also useful for user
graphics.\footnote{Script \texttt{patchEdgeInt2test.m}
verifies most of this code.}

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
the $\verb|nSubP1| \cdot \verb|nSubP2|
\cdot \verb|nPatch1| \cdot \verb|nPatch2|$ multiscale spatial grid on the
$\verb|nPatch1| \cdot \verb|nPatch2|$
array of patches.

\item \verb|patches| a struct set by \verb|configPatches2()|
which includes the following information.
\begin{itemize}

\item \verb|.x| is $\verb|nSubP1| \times1 \times1
\times1 \times \verb|nPatch1| \times1 $ array of the
spatial locations~$x_{iI}$ of the microscale grid points in
every patch. Currently it \emph{must} be an equi-spaced
lattice on the microscale index~$i$, but may be variable
spaced in macroscale index~$I$.

\item \verb|.y| is similarly $1\times \verb|nSubP2| 
\times1 \times1 \times1 \times \verb|nPatch2|$ array
of the spatial locations~$y_{jJ}$ of the microscale grid
points in every patch. Currently it \emph{must} be an
equi-spaced lattice on the microscale index~$j$, but may be
variable spaced in macroscale index~$J$.

\item \verb|.ordCC| is order of interpolation, currently
only $\{0,2,4,\ldots\}$

\item \verb|.periodic|, logical 2-array, indicates whether macroscale is
periodic domain, or alternatively that the macroscale has
left, right, top, bottom, front and back boundaries so
interpolation is via divided differences. 

\item \verb|.stag| in $\{0,1\}$ is one for staggered grid
(alternating) interpolation.  Currently must be zero.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the coupling
coefficients for finite width interpolation in each of the
$x,y$-directions---when invoking a periodic direction.

\item \verb|.EdgyInt|, true/false, for determining
patch-edge values by interpolation: true, from opposite-edge
next-to-edge values (often preserves symmetry); false, from
centre cross-patch values (near original scheme).

\item \verb|.nEdge|, three elements, the width of edge
values set by interpolation at the \(x,y\)-edge regions, 
respectively, of each patch (default is one all 
\(x,y\)-edges).

\item \verb|.nEnsem| the number of realisations in the
ensemble.

\item \verb|.parallel| whether serial or parallel.

\end{itemize}
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 6D array, $\verb|nSubP1| \cdot
\verb|nSubP2| \cdot \verb|nVars| \cdot
\verb|nEnsem| \cdot \verb|nPatch1| \cdot \verb|nPatch2|$, of the fields with edge values set by
interpolation.
\end{itemize}







\begin{devMan}

Determine the sizes of things. Any error arising in the
reshape indicates~\verb|u| has the wrong size.
\begin{matlab}
%}
[~,ny,~,~,~,Ny] = size(patches.y);
[nx,~,~,~,Nx,~] = size(patches.x);
nEnsem = patches.nEnsem;
nVars = round( numel(u)/numel(patches.x) ...
    /numel(patches.y)/nEnsem );
assert(numel(u) == nx*ny*Nx*Ny*nVars*nEnsem ...
  ,'patchEdgeInt2: input u has wrong size for parameters')
u = reshape(u,[nx ny nVars nEnsem Nx Ny]);
%{
\end{matlab}


\paragraph{Implement multiple width edges by folding}
Subsample~\(x,y\) coordinates, noting it is only
differences that count \emph{and} the microgrid~\(x,y\)
spacing must be uniform.
\begin{matlab}
%}
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
  u = reshape( permute(u,[2:2:4 1:2:3 5:8]) ...
             ,nx,ny,nVars,nEnsem,Nx,Ny);
end%if patches.nEdge
%{
\end{matlab}


\paragraph{Staggered grid}
Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches2|
tests there are an even number of patches). Then the
patch-ratio is effectively halved. The patch edges are near
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




\subsection{Interpolate over the two successive directions}
Only use the interior values of the fields for interpolating
to the edges.
\begin{matlab}
%}
u = u(2:nx-1,2:ny-1,:,:,:,:); 
%{
\end{matlab}
Interpolate in turn, the edge or mid-patch edges normal to
the \(x,y\)-directions, in this way we naturally fill-in
corner values.
\begin{matlab}
%}
u = patchEdgeIntCore(1,u,x,patches,stagShift ...
    ,1,nx,(ny-2)*nVars,nEnsem,1,Nx,Ny,patches.le,patches.ri);
u = patchEdgeIntCore(2,u,y,patches,stagShift ...
    ,nx,ny,      nVars,nEnsem,Nx,Ny,1,patches.bo,patches.to);
%{
\end{matlab}
Restore array~\verb|u| to its original shape.
\begin{matlab}
%}
u = reshape(u,nx,ny,nVars,nEnsem,Nx,Ny);
%{
\end{matlab}





\paragraph{Unfold multiple edges}  No need to restore~\(x,y,z\).
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
