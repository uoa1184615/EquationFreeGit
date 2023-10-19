% patchEdgeInt1() provides the interpolation across 1D space
% for 1D patches of simulations of a lattice system such as
% PDE discretisations.  AJR & JB, Sep 2018 -- 19 Oct 2023
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
maintains symmetry.  This function is primarily used by
\verb|patchSys1()| but is also useful for user graphics.
When using core averages (not fully implemented), assumes
the averages are sensible macroscale variables: then patch
edge values are determined by macroscale interpolation of
the core averages \citep{Bunder2013b}. \footnote{Script
\texttt{patchEdgeInt1test.m} verifies this code.}


Communicate patch-design variables via a second argument
(optional, except required for parallel computing of
\verb|spmd|), or otherwise via the global struct
\verb|patches|.
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
\emph{must} be an equi-spaced lattice on the microscale
index~$i$, but may be variable spaced in macroscale
index~$I$. 

\item \verb|.ordCC| is order of interpolation, integer~$\geq
-1$.

\item \verb|.periodic| indicates whether macroscale is
periodic domain, or alternatively that the macroscale has
left and right boundaries so interpolation is via divided
differences. 

\item \verb|.stag| in $\{0,1\}$ is one for staggered grid
(alternating) interpolation, and zero for ordinary grid.

\item \verb|.Cwtsr| and \verb|.Cwtsl| are the coupling
coefficients for finite width interpolation---when invoking
a periodic domain.

\item \verb|.EdgyInt|, true/false, for determining
patch-edge values by interpolation: 
true, from opposite-edge next-to-edge values (often
preserves symmetry); 
false, from centre-patch values (original scheme).

\item \verb|.nEdge|, for each patch, the number of edge
values set by interpolation at the edge regions of each
patch (default is one).

\item \verb|.nEnsem| the number of realisations in the
ensemble.

\item \verb|.parallel| whether serial or parallel.

\item \verb|.nCore| \todo{introduced sometime but not fully
implemented yet, because prefer ensemble, now disabled}

\end{itemize}
\end{itemize}



\paragraph{Output}
\begin{itemize}
\item \verb|u| is 4D array, $\verb|nSubP| \times
\verb|nVars| \times \verb|nEnsem| \times \verb|nPatch|$, of
the fields with edge values set by interpolation.
\end{itemize}







\begin{devMan}

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


\paragraph{Implement multiple width edges by folding}
Subsample~\(x\) coordinates, noting it is only differences
that count \emph{and} the microgrid~\(x\) spacing must be
uniform.
\begin{matlab}
%}
x = patches.x;
if patches.nEdge>1
  nEdge = patches.nEdge;
  x = x(1:nEdge:nx,:,:,:);
  nx = nx/nEdge;
  u = reshape(u,nEdge,nx,nVars,nEnsem,Nx);
  nVars = nVars*nEdge;
  u = reshape( permute(u,[2 1 3:5]) ,nx,nVars,nEnsem,Nx);
end%if patches.nEdge
%{
\end{matlab}


\paragraph{Staggered grid}
Deal with staggered grid by doubling the number of fields
and halving the number of patches (\verb|configPatches1()|
tests that there are an even number of patches). Then the
patch-ratio is effectively halved. The patch edges are near
the middle of the gaps and swapped. \todo{Have not yet
tested whether works for Edgy Interpolation.}  \todo{Have
not yet implemented multiple edge values for a staggered
grid as I am uncertain whether it makes any sense. }
\begin{matlab}
%}
  if patches.stag % transform by doubling the number of fields
  error('Doubt staggered works following changes made in Oct 2023??')
    v = nan(size(u)); % currently to restore the shape of u
    u = [u(:,:,:,1:2:Nx) u(:,:,:,2:2:Nx)];
    stagShift = 0.5*[ones(1,nVars) -ones(1,nVars)];
    iV = [nVars+1:2*nVars 1:nVars]; % scatter interp to alternate field
    patches.ratio = patches.ratio/2;   % ratio effectively halved
    Nx = Nx/2; % halve the number of patches
    nVars = nVars*2;   % double the number of fields
  else % the values for standard spectral
    stagShift = 0;  
    iV = 1:nVars;
  end%if patches.stag
%{
\end{matlab}


\subsection{Interpolate in the x-direction}
Only use the interior values of the fields for interpolating
to the edges.
\begin{matlab}
%}
u = u(2:nx-1,:,:,:); 
%{
\end{matlab}
Interpolate in turn, the edge or mid-patch edges normal to
the \(x,y\)-directions, in this way we naturally fill-in
corner values.
\begin{matlab}
%}
u = patchEdgeIntCore(1,u,x,patches,stagShift ...
    ,1,nx,nVars,nEnsem,1,Nx,1,patches.le,patches.ri);
%{
\end{matlab}
Restore array~\verb|u| to its original shape.
\begin{matlab}
%}
u = reshape(u,nx,nVars,nEnsem,Nx);
%{
\end{matlab}


\paragraph{Unfold staggered Grid}
\begin{matlab}
%}
if patches.stag % transform by doubling the number of fields
    patches.ratio = patches.ratio*2;   % restore ratio
end
%{
\end{matlab}



\paragraph{Unfold multiple edges}  No need to restore~\(x\).
\begin{matlab}
%}
if patches.nEdge>1
  nVars = nVars/nEdge;
  u = reshape( u ,nx,nEdge,nVars,nEnsem,Nx);
  nx = nx*nEdge;
  u = reshape( permute(u,[2 1 3:5]) ,nx,nVars,nEnsem,Nx);
end%if patches.nEdge
%{
\end{matlab}

Fin, returning the 4D array of field values.  
\begin{matlab}
%}
end% function patchEdgeInt1
%{
\end{matlab}
\end{devMan}
%}