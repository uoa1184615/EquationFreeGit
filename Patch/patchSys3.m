% patchSys3() provides an interface to time integrators
% for the dynamics on patches in 3D coupled across space.
% The system must be a lattice system such as PDE
% discretisations. AJR, Aug 2020 -- 12 Apr 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchSys3()}: interface 3D space to time integrators}
\label{sec:patchSys3}


To simulate in time with 3D spatial patches we often need to
interface a users time derivative function with time
integration routines such as \verb|ode23| or~\verb|PIRK2|.
This function provides an interface.  Communicate
patch-design variables (\cref{sec:configPatches3}) either
via the global struct~\verb|patches| or via an optional
third argument.  \verb|patches| is required for the parallel
computing of \verb|spmd|, or if parameters are to be passed
though to the user microscale function.

\begin{matlab}
%}
function dudt = patchSys3(t,u,patches,varargin)
if nargin<3, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector\slash array of length
$\verb|prod(nSubP)| \cdot \verb|nVars| \cdot \verb|nEnsem|
\cdot \verb|prod(nPatch)|$ where there are $\verb|nVars|
\cdot \verb|nEnsem|$ field values at each of the points in
the $\verb|nSubP(1)| \times \verb|nSubP(2)| \times
\verb|nSubP(3)| \times \verb|nPatch(1)| \times
\verb|nPatch(2)| \times \verb|nPatch(3)|$ spatial grid.

\item \verb|t| is the current time to be passed to the
user's time derivative function.

\item \verb|patches| a struct set by \verb|configPatches3()|
with the following information used here.

\begin{itemize}
\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,patches,...)| that computes the time
derivatives on the patchy lattice.  The array~\verb|u| has
size $\verb|nSubP(1)| \times \verb|nSubP(2)| \times
\verb|nSubP(3)| \times \verb|nVars| \times \verb|nEsem|
\times \verb|nPatch(1)| \times \verb|nPatch(2)| \times
\verb|nPatch(3)|$.  Time derivatives must be computed into
the same sized array, although herein the patch edge-values
are overwritten by zeros.

\item \verb|.x| is $\verb|nSubP(1)| \times1 \times1 \times1
\times1 \verb|nPatch(1)| \times1 \times1$ array of the
spatial locations~$x_{i}$ of the microscale
$(i,j,k)$-grid points in every patch.  Currently it
\emph{must} be an equi-spaced lattice on both macro- and
microscales.

\item \verb|.y| is similarly $1 \times \verb|nSubP(2)|
\times1 \times1 \times1 \times1 \times \verb|nPatch(2)|
\times1$ array of the spatial locations~$y_{j}$ of the
microscale $(i,j,k)$-grid points in every patch. 
Currently it \emph{must} be an equi-spaced lattice on both
macro- and microscales.

\item \verb|.z| is similarly $1 \times1 \times
\verb|nSubP(3)| \times1 \times1 \times1 \times1 \times
\verb|nPatch(3)|$ array of the spatial locations~$z_{k}$
of the microscale $(i,j,k)$-grid points in every patch.
Currently it \emph{must} be an equi-spaced lattice on both
macro- and microscales.

\end{itemize}

\item \verb|varargin|, optional, is arbitrary list of
parameters to be passed onto the users time-derivative
function as specified in configPatches3.
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is a vector\slash array of of time
derivatives, but with patch edge-values set to zero.  It is
of total length $\verb|prod(nSubP)| \cdot \verb|nVars|
\cdot \verb|nEnsem| \cdot \verb|prod(nPatch)|$ and the same
dimensions as~\verb|u|.
\end{itemize}



\begin{devMan}
Sets the edge-face values from macroscale interpolation of
centre-patch values, and if necessary, reshapes the
fields~\verb|u| as a 8D-array.  \cref{sec:patchEdgeInt3}
describes \verb|patchEdgeInt3()|.
\begin{matlab}
%}
sizeu = size(u);
u = patchEdgeInt3(u,patches);
%{
\end{matlab}

Ask the user function for the time derivatives computed in
the array, overwrite its edge\slash face values with the
dummy value of zero (as \verb|ode15s| chokes on NaNs), then
return to the user\slash integrator as same sized array as
input.
\begin{matlab}
%}
dudt = patches.fun(t,u,patches,varargin{:});
m = patches.nEdge(1);
dudt([1:m end-m+1:end],:,:,:) = 0;
m = patches.nEdge(2);
dudt(:,[1:m end-m+1:end],:,:) = 0;
m = patches.nEdge(3);
dudt(:,:,[1:m end-m+1:end],:) = 0;
dudt = reshape(dudt,sizeu);
%{
\end{matlab}
Fin.
\end{devMan}
%}
 