% patchSys1() provides an interface to time integrators
% for the dynamics on patches  coupled across space. The
% system must be a smooth lattice system such as PDE
% discretisations.  AJR, Nov 2017 -- Feb 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchSys1()}: interface 1D space to time integrators}
\label{sec:patchSys1}


To simulate in time with 1D spatial patches we often need to
interface a user's time derivative function with time
integration routines such as \verb|ode23| or~\verb|PIRK2|.
This function provides an interface.  It mostly assumes that
the sub-patch structure is \emph{smooth enough} so that the
patch centre-values are sensible macroscale variables, and
patch edge values are determined by macroscale interpolation
of the patch-centre or edge values. Nonetheless, microscale
heterogeneous systems may be accurately simulated with this
function via appropriate interpolation.  Communicate
patch-design variables (\cref{sec:configPatches1}) either
via the global struct~\verb|patches| or via an optional
third argument (except that this last is required for
parallel computing of \verb|spmd|).

\begin{matlab}
%}
function dudt=patchSys1(t,u,patches,varargin)
if nargin<3, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector\slash array of length
$\verb|nSubP| \cdot \verb|nVars| \cdot \verb|nEnsem| \cdot
\verb|nPatch|$ where there are $\verb|nVars| \cdot
\verb|nEnsem|$ field values at each of the points in the
$\verb|nSubP|\times \verb|nPatch|$ grid.

\item \verb|t| is the current time to be passed to the
user's time derivative function.

\item \verb|patches| a struct set by \verb|configPatches1()|
with the following information  used here.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,patches)| that computes the time derivatives
on the patchy lattice. The array~\verb|u| has size
$\verb|nSubP| \times \verb|nVars| \times \verb|nEnsem|
\times \verb|nPatch|$.  Time derivatives should be computed
into the same sized array, then herein the patch edge
values are overwritten by zeros.

\item \verb|.x| is $\verb|nSubP| \times1 \times1 \times
\verb|nPatch|$ array of the spatial locations~$x_{i}$ of
the microscale grid points in every patch.  Currently it
\emph{must} be an equi-spaced lattice on the microscale.
\end{itemize}

\item \verb|varargin| is arbitrary number of parameters to
be passed onto the users time-derivative function as
specified in configPatches1.
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is a vector\slash array of of time
derivatives, but with patch edge-values set to zero. It is
of total length $\verb|nSubP| \cdot \verb|nVars| \cdot
\verb|nEnsem| \cdot \verb|nPatch|$ and the same dimensions
as~\verb|u|.
\end{itemize}



\begin{devMan}
Reshape the fields~\verb|u| as a 4D-array, and sets the edge
values from macroscale interpolation of centre-patch values.
\cref{sec:patchEdgeInt1} describes \verb|patchEdgeInt1()|.
\begin{matlab}
%}
sizeu = size(u);
u = patchEdgeInt1(u,patches);
%{
\end{matlab}

Ask the user function for the time derivatives computed in
the array, overwrite its edge values with the dummy value of
zero (as \verb|ode15s| chokes on NaNs), then return to the
user\slash integrator as same sized array as input.
\begin{matlab}
%}
dudt=patches.fun(t,u,patches,varargin{:});
dudt([1 end],:,:,:) = 0;
dudt=reshape(dudt,sizeu);
%{
\end{matlab}
Fin.
\end{devMan}
%}
 