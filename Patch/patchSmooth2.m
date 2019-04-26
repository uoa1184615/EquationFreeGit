% Provides an interface to time integrators for the dynamics
% on patches in 2D coupled across space. The system must be
% a smooth lattice system such as PDE discretisations.
% AJR, Nov 2018
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchSmooth2()}: interface to time integrators}
\label{sec:patchSmooth2}
%\localtableofcontents


To simulate in time with spatial patches we often need to
interface a users time derivative function with time
integration routines such as \verb|ode15s| or~\verb|PIRK2|.
This function provides an interface. It assumes that the
sub-patch structure is \emph{smooth} so that the patch
centre-values are sensible macroscale variables, and patch
edge-values are determined by macroscale interpolation of
the patch-centre values. Communicate patch-design variables
to this function via the previously established global
struct~\verb|patches|.
\begin{matlab}
%}
function dudt = patchSmooth2(t,u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector of length \(\verb|prod(nSubP)|
\cdot \verb|prod(nPatch)| \cdot \verb|nVars|\) where there
are \verb|nVars| field values at each of the points in the
\(\verb|nSubP(1)| \times \verb|nSubP(2)|\times
\verb|nPatch(1)| \times \verb|nPatch(2)|\) grid.

\item \verb|t| is the current time to be passed to the
user's time derivative function.

\item \verb|patches| a struct set by \verb|configPatches2()|
with the following information used here.
\begin{itemize}
\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,x,y)| that computes the time derivatives on
the patchy lattice.  The array~\verb|u| has size
\(\verb|nSubP(1)| \times\verb|nSubP(2)| \times
\verb|nPatch(1)| \times \verb|nPatch(2)| \times
\verb|nVars|\).  Time derivatives must be computed into the
same sized array, but herein the patch edge-values are
overwritten by zeros.

\item \verb|.x| is \(\verb|nSubP(1)| \times
\verb|nPatch(1)|\) array of the spatial locations~\(x_{ij}\)
of the microscale grid points in every patch. Currently it
\emph{must} be an equi-spaced lattice on both macro- and
microscales.

\item \verb|.y| is similarly \(\verb|nSubP(2)| \times
\verb|nPatch(2)|\) array of the spatial locations~\(y_{ij}\)
of the microscale grid points in every patch.  Currently it
\emph{must} be an equi-spaced lattice on both macro- and
microscales.

\end{itemize}
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is \(\verb|prod(nSubP)| \cdot
\verb|prod(nPatch)| \cdot \verb|nVars|\) vector of time
derivatives, but with patch edge-values set to zero.
\end{itemize}

\begin{devMan}
Reshape the fields~\verb|u| as a 4/5D-array, and sets the
edge values from macroscale interpolation of centre-patch
values.  \cref{sec:patchEdgeInt2} describes
\verb|patchEdgeInt2()|.
\begin{matlab}
%}
u = patchEdgeInt2(u);
%{
\end{matlab}

Ask the user function for the time derivatives computed in
the array, overwrite its edge values with the dummy value of
zero, then return to a to the user\slash integrator as
column vector.
\begin{matlab}
%}
dudt = patches.fun(t,u,patches.x,patches.y);
dudt([1 end],:,:,:,:) = 0;
dudt(:,[1 end],:,:,:) = 0;
dudt = reshape(dudt,[],1);
%{
\end{matlab}
Fin.
\end{devMan}
%}
 