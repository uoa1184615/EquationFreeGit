% Compute the weightings for the polynomial
% interpolation of field values for coupling.
% AJR, 7 Aug 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{patchCwts}: weights of polynomial interpolation}
\label{sec:patchCwts}



\subsection{Introduction}
Computes the weightings for the polynomial interpolation of
field values for inter-patch coupling.  Should work for any
number of dimensions as determined by the number of elements
in parameter \verb|ratio|.  Used by
\verb|configPatches|\(n\) for \(n=1,2,\ldots\)\,.
\begin{matlab}
%}
function patchCwts(ratio,ordCC,stag)
global patches
%{
\end{matlab}
\paragraph{Input} \begin{itemize}
\item \verb|ratio| row vector, one element for each axis of
the spatial domain, of either the half-width or full-width
of a patch to the spacing of the patch mid-points along that
axis direction.

\item \verb|ordCC| is the order of the polynomial
interpolation for inter-patch coupling across empty space of
the macroscale patch values to the edge-values of the
patches: must be~\(2,4,\ldots\)\,.

\item \verb|stag| is true for interpolation using only odd
neighbouring patches as for staggered grids, and false for
the usual case of all neighbour coupling---as yet only
tested in 1D.

\end{itemize}

\paragraph{Output} The \emph{global} struct \verb|patches|
has the following components added.
\begin{itemize}
\item \verb|.Cwtsr| and \verb|.Cwtsl|, when \(n\)~is the
number of elements of \verb|ratio|, are the
\(\verb|ordCC|\times n\)-array of weights for the
inter-patch interpolation onto the `right' edges and `left'
edges (respectively) with patch:macroscale ratio as
specified.
\end{itemize}


\begin{devMan}
First check, reserve storage, and define some index vectors.
\begin{matlab}
%}
assert(ordCC>0,'order of poly interp must be positive')
patches.Cwtsr=nan(ordCC,numel(ratio));
ks = (1:2:ordCC)';   
ps = (1:ordCC/2)'-1; 
%{
\end{matlab}
If staggered grid, then we need something like equation~(7)
in \cite{Cao2014a}.  But so far only tested for 1D??
\begin{matlab}
%}
if stag 
    patches.Cwtsr(ks  ,:) = [ones(size(ratio)) ...
      cumprod( (ratio.^2-ks(1:end-1).^2)/4 ,1) ...
      ./factorial(2*ps(1:end-1)) ];    
    patches.Cwtsr(ks+1,:) = [ratio/2 ...
      cumprod( (ratio.^2-ks(1:end-1).^2)/4 ,1) ...
      ./factorial(2*ps(1:end-1)+1).*ratio/2 ];
%{
\end{matlab}
For non-staggered edge-to-edge or centre-to-edge
interpolation, use these weights.
\begin{matlab}
%}
else 
    patches.Cwtsr(ks  ,:) = cumprod(ratio.^2-ps.^2,1) ...
        ./factorial(2*ps+1)./ratio;
    patches.Cwtsr(ks+1,:) = cumprod(ratio.^2-ps.^2,1) ...
        ./factorial(2*ps+2);
end
patches.Cwtsl = (-1).^((1:ordCC)'-patches.stag).*patches.Cwtsr;
%{
\end{matlab}
Fin.
\end{devMan}
%}
