%Creates a data struct of the design of patches for later
%use by the patch functions such as smoothPatch1() 
%AJR, Nov 2017
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{makePatches()}: makes the spatial patches for the suite}
\label{sec:makePatches}

Constructs the struct~\verb|patches| for use by the patch scheme~\verb|patchDt|.

\begin{matlab}
%}
function makePatches(fun,Xa,Xb,nPatch,ratio,nSubP)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|fun| is the name of the user function, \verb|fun(t,u,x)|, that will compute time derivatives of quantities on the patches.
\item \verb|Xa,Xb| give the macro-space domain of the computation: patches are spread evenly over the interior of the interval~\([\verb|Xa|,\verb|Xb|]\).
Currently the system is assumed macro-periodic in this domain.
\item \verb|nPatch| is the number of evenly spaced patches.
\item \verb|ratio| (real) is the ratio of the half-width of a patch to the spacing of the patch mid-points: so \(\verb|ratio|=\tfrac12\) means the patches abut; and \(\verb|ratio|=1\) is overlapping patches as in holistic discretisation.
\item \verb|nSubP| is the number of microscale lattice points in each patch.  Must be odd so that there is a central lattice point.
\end{itemize}

\paragraph{Output} The \emph{global} struct \verb|patches| is created and set.
\begin{itemize}
\item \verb|patches.fun| is the name of the user's function \verb|fun(u,t,x)| that computes the time derivatives on the patchy lattice. \item \verb|patches.x| is \(\verb|nSubP|\times \verb|nPatch|\) array of the regular spatial locations~\(x_{ij}\) of the microscale grid points in every patch.  
\end{itemize}


First, just store the pointer to the time derivative function in the struct.
\begin{matlab}
%}
patches.fun=fun;
%{
\end{matlab}

Second, set the centre of the patches in a the macroscale grid of patches assuming periodic macroscale domain.
\begin{matlab}
%}
X=linspace(Xa,Xb,nPatch+1);
X=X(1:nPatch)+diff(X)/2;
DX=X(2)-X(1);
%{
\end{matlab}
Construct the microscale in each patch, assuming Dirichlet patch edges, and a half-patch length of~\(\verb|ratio|\cdot\verb|DX|\).
\begin{matlab}
%}
if mod(nSubP,2)==0, error('makePatches: nSubP must be odd'), end
i0=(nSubP+1)/2;
dx=ratio*DX/(i0-1);
patches.x=bsxfun(@plus,dx*(-i0+1:i0-1)',X); % micro-grid
%{
\end{matlab}

Fin.
%}
