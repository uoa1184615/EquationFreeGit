%Provides the coupling across space for patches of
%simulations of a smooth lattice system such as PDE
%discretisations.
%AJR, Nov 2017 -- Sep 2018
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{patchSmooth1()}}
\label{sec:patchSmooth1}
\localtableofcontents

Couples patches across space so a spatially discrete system can be integrated in time via the patch or gap-tooth scheme \cite[]{Roberts06d}.
Assumes that the sub-patch structure is \emph{smooth} so that the patch centre-values are sensible macroscale variables, and patch edge values are determined by macroscale interpolation of the patch-centre values. 
Need to pass patch-design variables to this function, so use the global struct~\verb|patches|.
\begin{matlab}
%}
function dudt=patchSmooth1(t,u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|u| is a vector of length \(\verb|nSubP|\cdot \verb|nPatch|\cdot \verb|nVars|\) where there are \verb|nVars| field values at each of the points in the \(\verb|nSubP|\times \verb|nPatch|\) grid.
\item \verb|t| is the current time to be passed to the user's time derivative function.
\item \verb|patches| a struct set by \verb|makePatches()| with the following information that is used here.
\begin{itemize}
\item \verb|.fun| is the name of the user's function \verb|fun(t,u,x)| that computes the time derivatives on the patchy lattice. 
The array~\verb|u| has size \(\verb|nSubP|\times \verb|nPatch|\times \verb|nVars|\).
Time derivatives must be computed into the same sized array, but the patch edge values will be overwritten by zeros.
\item \verb|.x| is \(\verb|nSubP|\times \verb|nPatch|\) array of the spatial locations~\(x_{ij}\) of the microscale grid points in every patch.  
Currently it \emph{must} be a regular lattice on both macro- and micro-scales.
\end{itemize}

\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is \(\verb|nSubP|\cdot \verb|nPatch|\cdot \verb|nVars|\) vector of time derivatives, but with zero on patch edges??
\end{itemize}

Reshape the fields~\verb|u| as a 2/3D-array, and sets the edge values from macroscale interpolation of centre-patch values.
\S\ref{sec:patchEdgeInt1} describes function \verb|patchEdgeInt1()|.
\begin{matlab}
%}
u=patchEdgeInt1(u);
%{
\end{matlab}

Ask the user for the time derivatives computed in the array, overwrite its edge values with the dummy value of zero, then return to an integrator as column vector.
\begin{matlab}
%}
dudt=patches.fun(t,u,patches.x);
dudt([1 end],:,:)=0;
dudt=reshape(dudt,[],1);
%{
\end{matlab}
Fin.
%}
 