%Creates a data struct of the design of patches for later
%use by the patch functions such as smoothPatch1() 
%AJR, Nov 2017
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{makePatches()}: makes the spatial patches for the suite}
\label{sec:makePatches}
\localtableofcontents

Makes the struct~\verb|patches| for use by the patch\slash gap-tooth time derivative function~\verb|patchSmooth1()|.

\begin{matlab}
%}
function makePatches(fun,Xa,Xb,nPatch,ordCC,ratio,nSubP)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|fun| is the name of the user function, \verb|fun(t,u,x)|, that will compute time derivatives of quantities on the patches.
\item \verb|Xa,Xb| give the macro-space domain of the computation: patches are spread evenly over the interior of the interval~\([\verb|Xa|,\verb|Xb|]\).
Currently the system is assumed macro-periodic in this domain.
\item \verb|nPatch| is the number of evenly spaced patches.
\item \verb|ordCC| is the order of interpolation across empty space of the macroscale mid-patch values to the edge of the patches for inter-patch coupling: currently must be in~\(\{2,4,6\}\).
\item \verb|ratio| (real) is the ratio of the half-width of a patch to the spacing of the patch mid-points: so \(\verb|ratio|=\tfrac12\) means the patches abut; and \(\verb|ratio|=1\) is overlapping patches as in holistic discretisation.
\item \verb|nSubP| is the number of microscale lattice points in each patch.  Must be odd so that there is a central lattice point.
\end{itemize}

\paragraph{Output} The \emph{global} struct \verb|patches| is created and set.
\begin{itemize}
\item \verb|patches.fun| is the name of the user's function \verb|fun(u,t,x)| that computes the time derivatives on the patchy lattice. 
\item \verb|patches.ordCC| is the specified order of inter-patch coupling. 
\item \verb|patches.alt| is true for interpolation using only odd neighbouring patches as for staggered grids, and false for the usual case of all neighbour coupling.
\item \verb|patches.Cwtsr| and \verb|.Cwtsl| are the \(\verb|ordCC|\)-vector of weights for the inter-patch interpolation onto the right and left edges (respectively) with patch:macroscale ratio as specified.
\item \verb|patches.x| is \(\verb|nSubP|\times \verb|nPatch|\) array of the regular spatial locations~\(x_{ij}\) of the microscale grid points in every patch.  
\end{itemize}


First, store the pointer to the time derivative function in the struct.
\begin{matlab}
%}
patches.fun=fun;
%{
\end{matlab}

Second, store the order of interpolation that is to provide the values for the inter-patch coupling conditions.
Maybe allow \verb|ordCC| of~0 and~\(-1\) to request spectral coupling??
\begin{matlab}
%}
if ~ismember(ordCC,[1:8])
    error('makePatch: ordCC out of allowed range [1:8]')
end
%{
\end{matlab}
For odd~\verb|ordCC| do interpolation based upon odd neighbouring patches as is useful for staggered grids.
\begin{matlab}
%}
patches.alt=rem(ordCC,2);
ordCC=ordCC+patches.alt;
patches.ordCC=ordCC;
%{
\end{matlab}
Might as well precompute the weightings for the interpolation of field values for coupling. (What about coupling via derivative values??  what about spectral coupling??)
\begin{matlab}
%}
if patches.alt  % eqn (7) in \cite{Cao2014a}
  patches.Cwtsr=[1
    ratio/2
    (-1+ratio^2)/8
    (-1+ratio^2)*ratio/48
    (9-10*ratio^2+ratio^4)/384
    (9-10*ratio^2+ratio^4)*ratio/3840
    (-225+259*ratio^2-35*ratio^4+ratio^6)/46080
    (-225+259*ratio^2-35*ratio^4+ratio^6)*ratio/645120 ];
else % 
  patches.Cwtsr=[ratio
    ratio^2/2
    (-1+ratio^2)*ratio/6
    (-1+ratio^2)*ratio^2/24
    (4-5*ratio^2+ratio^4)*ratio/120
    (4-5*ratio^2+ratio^4)*ratio^2/720
    (-36+49*ratio^2-14*ratio^4+ratio^6)*ratio/5040
    (-36+49*ratio^2-14*ratio^4+ratio^6)*ratio^2/40320 ];
end
patches.Cwtsr=patches.Cwtsr(1:ordCC);
patches.Cwtsl=(-1).^((1:ordCC)'-patches.alt).*patches.Cwtsr;
%{
\end{matlab}


Third, set the centre of the patches in a the macroscale grid of patches assuming periodic macroscale domain.
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
