%Provides an interface to time integrators for the dynamics
%on patches  coupled across space. The system must be a
%smooth lattice system such as PDE discretisations.
%AJR, Nov 2017 -- Sep 2018
%!TEX root = ../Doc/equationFreeDoc.tex
%{
\subsection{\texttt{patchCoreSmooth1()}: interface to time integrators}
\label{sec:patchSmooth1}
\localtableofcontents

To simulate in time with spatial patches we often need to
interface a users time derivative function with time
integration routines such as \verb|ode15s| or~\verb|PIRK2|.
This function provides an interface. Either the sub-patch
structure is \emph{smooth} so that the patch
centre-values are sensible macroscale variables, or the 
user chooses to average over a \emph{core} of 
values in the centre of each patch with these averages providing
sensible macroscale variables (no core averaging corresponds 
to a core of size one). Patch
edge values are determined by macroscale interpolation of
the patch-centre values or core-averaged values. 
Communicate patch-design variables
to this function using the previously established global
struct~\verb|patches|.
\begin{matlab}
%}
function dudt=patchCoreSmooth1(t,u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|u| is a vector of length \(\verb|nSubP|\cdot
\verb|nPatch|\cdot \verb|nVars|\) where there are
\verb|nVars| field values at each of the points in the
\(\verb|nSubP|\times \verb|nPatch|\) grid.
\item \verb|t| is the current time to be passed to the
user's time derivative function.
\item \verb|patches| a struct set by \verb|configPatches1()|
with the following information  used here.
\begin{itemize}
\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,x)| that computes the time derivatives on the
patchy lattice. The array~\verb|u| has size
\(\verb|nSubP|\times \verb|nPatch|\times \verb|nVars|\).
Time derivatives must be computed into the same sized array,
but herein the patch edge values are overwritten by zeros.
\item \verb|.x| is \(\verb|nSubP|\times \verb|nPatch|\)
array of the spatial locations~\(x_{ij}\) of the microscale
grid points in every patch. Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.
\end{itemize}
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is \(\verb|nSubP|\cdot \verb|nPatch|\cdot
\verb|nVars|\) vector of time derivatives, but with patch
edge values set to zero.
\end{itemize}


\begin{body}
Reshape the fields~\verb|u| as a 2/3D-array, and sets the
edge values from macroscale interpolation of centre-patch
values. \cref{sec:patchEdgeInt1} describes
\verb|patchEdgeInt1()|.
\begin{matlab}
%}
u=patchCoreEdgeInt1(u);
%{
\end{matlab}

Ask the user function for the time derivatives computed in
the array, overwrite its edge values with the dummy value of
zero, then return to a to the user\slash integrator as
column vector.
\begin{matlab}
%}
dudt=patches.fun(t,u,patches.x);
dudt([1 end],:,:)=0;
dudt=reshape(dudt,[],1);
%{
\end{matlab}
Fin.
\end{body}
%}
 