% Provides an interface to time integrators for the dynamics
% on patches in 3D coupled across space. The system must be
% a smooth lattice system such as PDE discretisations.
% AJR, Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchSmooth3()}: interface to time integrators}
\label{sec:patchSmooth3}
%\localtableofcontents


To simulate in time with spatial patches we often need to
interface a users time derivative function with time
integration routines such as \verb|ode15s| or~\verb|PIRK2|.
This function provides an interface. It assumes that the
sub-patch structure is \emph{smooth enough} so that the patch
centre-values are sensible macroscale variables, and patch
edge-values are determined by macroscale interpolation of
the patch-centre values. Communicate patch-design variables
to this function via the previously established global
struct~\verb|patches|.
\begin{matlab}
%}
function dudt = patchSmooth3(t,u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector of length \(\verb|prod(nSubP)|
\cdot \verb|nVars| \cdot \verb|nEnsem| \cdot
\verb|prod(nPatch)|\) where there are \(\verb|nVars|\cdot
\verb|nEnsem|\) field values at each of the points in the
\(\verb|nSubP(1)| \times \verb|nSubP(2)|\times
\verb|nSubP(3)|\times \verb|nPatch(1)| \times
\verb|nPatch(2) \times \verb|nPatch(3)|\) grid.

\item \verb|t| is the current time to be passed to the
user's time derivative function.

\item \verb|patches| a struct set by \verb|configPatches3()|
with the following information used here.
\begin{itemize}
\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,x,y,z)| that computes the time derivatives on
the patchy lattice.  The array~\verb|u| has size
\(\verb|nSubP(1)| \times\verb|nSubP(2)|
\times\verb|nSubP(3)| \times \verb|nVars| \times
\verb|nEsem| \times \verb|nPatch(1)| \times \verb|nPatch(2)|
\times \verb|nPatch(3)|\).  Time derivatives must be
computed into the same sized array, but herein the patch
edge-values are overwritten by zeros.

\item \verb|.x| is \(\verb|nSubP(1)| \times1 \times1 \times1
\times1 \times1 \verb|nPatch(1)| \times1 \times1\) array of
the spatial locations~\(x_{ij}\) of the microscale grid
points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and microscales.

\item \verb|.y| is similarly \(1 \times \verb|nSubP(2)|
\times1 \times1 \times1 \times1 \times \verb|nPatch(2)|
\times1\) array of the spatial locations~\(y_{ij}\) of the
microscale grid points in every patch.  Currently it
\emph{must} be an equi-spaced lattice on both macro- and
microscales.

\item \verb|.z| is similarly \(1 \times \times1
\verb|nSubP(3)| \times1 \times1 \times1 \times1 \times
\verb|nPatch(3)|\) array of the spatial locations~\(z_{ij}\)
of the microscale grid points in every patch.  Currently it
\emph{must} be an equi-spaced lattice on both macro- and
microscales.

\end{itemize}
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is \(\verb|prod(nSubP)| \cdot \verb|nVars|
\cdot \verb|nEnsem| \cdot \verb|prod(nPatch)|\) vector of
time derivatives, but with patch edge-values set to zero.
\end{itemize}

\begin{devMan}
Reshape the fields~\verb|u| as a 8D-array, and sets the
edge values from macroscale interpolation of centre-patch
values.  \cref{sec:patchEdgeInt3} describes
\verb|patchEdgeInt3()|.
\begin{matlab}
%}
u = patchEdgeInt3(u);
%{
\end{matlab}

Ask the user function for the time derivatives computed in
the array, overwrite its edge values with the dummy value of
zero, then return to a to the user\slash integrator as
column vector.
\begin{matlab}
%}
dudt = patches.fun(t,u,patches.x,patches.y,patches.z);
dudt([1 end],:,:,:) = 0;
dudt(:,[1 end],:,:) = 0;
dudt(:,:,[1 end],:) = 0;
dudt = reshape(dudt,[],1);
%{
\end{matlab}
Fin.
\end{devMan}
%}
 