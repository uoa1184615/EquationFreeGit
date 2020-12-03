% patchSmooth2() Provides an interface to time integrators
% for the dynamics on patches in 2D coupled across space.
% The system must be a smooth lattice system such as PDE
% discretisations. AJR, Nov 2018 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchSmooth2()}: interface 2D space to time integrators}
\label{sec:patchSmooth2}


To simulate in time with 2D spatial patches we often need to
interface a users time derivative function with time
integration routines such as \verb|ode15s| or~\verb|PIRK2|.
This function provides an interface. It assumes that the
sub-patch structure is \emph{smooth enough} so that the
patch centre-values are sensible macroscale variables, and
patch edge-values are determined by macroscale interpolation
of the patch-centre or edge values. Nonetheless, microscale
heterogeneous systems may be accurately simulated with this
function via appropriate interpolation. Communicate
patch-design variables (\cref{sec:configPatches2}) either
via the global struct~\verb|patches| or via an optional
third argument (except that this last is required for
parallel computing of \verb|spmd|).

\begin{matlab}
%}
function dudt = patchSmooth2(t,u,patches)
if nargin<3, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector\slash array of length
\(\verb|prod(nSubP)| \cdot \verb|nVars| \cdot \verb|nEnsem|
\cdot \verb|prod(nPatch)|\) where there are \(\verb|nVars|
\cdot \verb|nEnsem|\) field values at each of the points in
the \(\verb|nSubP(1)| \times \verb|nSubP(2)| \times
\verb|nPatch(1)| \times \verb|nPatch(2)|\) grid.

\item \verb|t| is the current time to be passed to the
user's time derivative function.

\item \verb|patches| a struct set by \verb|configPatches2()|
with the following information used here.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,patches)| that computes the time derivatives
on the patchy lattice.  The array~\verb|u| has size
\(\verb|nSubP(1)| \times \verb|nSubP(2)| \times \verb|nVars|
\times \verb|nEsem| \times \verb|nPatch(1)| \times
\verb|nPatch(2)|\).  Time derivatives must be computed into
the same sized array, although herein the patch edge-values
are overwritten by zeros.

\item \verb|.x| is \(\verb|nSubP(1)| \times1 \times1 \times1
\verb|nPatch(1)| \times1\) array of the spatial
locations~\(x_{i}\) of the microscale \((i,j)\)-grid points
in every patch.  Currently it \emph{must} be an equi-spaced
lattice on both macro- and micro-scales.

\item \verb|.y| is similarly \(1 \times \verb|nSubP(2)|
\times1 \times1 \times1 \times \verb|nPatch(2)|\) array of
the spatial locations~\(y_{j}\) of the microscale
\((i,j)\)-grid points in every patch.  Currently it
\emph{must} be an equi-spaced lattice on both macro- and
micro-scales.

\end{itemize}
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is a vector\slash array of of time
derivatives, but with patch edge-values set to zero. It is
of total length \(\verb|prod(nSubP)| \cdot \verb|nVars|
\cdot \verb|nEnsem| \cdot \verb|prod(nPatch)|\) and the same
dimensions as~\verb|u|.
\end{itemize}



\begin{devMan}
Reshape the fields~\verb|u| as a 6D-array, and sets the edge
values from macroscale interpolation of centre-patch values.
 \cref{sec:patchEdgeInt2} describes \verb|patchEdgeInt2()|.
\begin{matlab}
%}
sizeu = size(u);
u = patchEdgeInt2(u,patches);
%{
\end{matlab}

Ask the user function for the time derivatives computed in
the array, overwrite its edge values with the dummy value of
zero, then return to the user\slash integrator as same sized
array as input.
\begin{matlab}
%}
dudt = patches.fun(t,u,patches);
dudt([1 end],:,:,:,:,:) = 0;
dudt(:,[1 end],:,:,:,:) = 0;
dudt = reshape(dudt,sizeu);
%{
\end{matlab}
Fin.
\end{devMan}
%}
 