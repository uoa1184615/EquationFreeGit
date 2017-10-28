%Provides the coupling across space for patches of
%simulations of a smooth lattice system such as PDE
%discretisations.
%AJR, Nov 2017
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{patchSmooth1()}}
\label{sec:patchSmooth1}

Couples patches across space so a spatially discrete system can be integrated in time via the patch or gap-tooth scheme \cite[]{Roberts06d}.
Need to pass patch-design variables to this function, so use the global struct~\verb|patches|.
\begin{matlab}
%}
function dudt=patchSmooth1(t,u)
global patches
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|u| is a vector of size \(\verb|nSubP|\cdot \verb|nPatch|\cdot \verb|nVars|\) where there are \verb|nVars| field values at each of the points in the \(\verb|nSubP|\times \verb|nPatch|\) grid.
\item \verb|t| is the current time to be passed to the user's time derivative function.
\item \verb|patches.fun| is the name of the user's function \verb|fun(t,u,x)| that computes the time derivatives on the patchy lattice. The array~\verb|u| has size \(\verb|nSubP|\cdot \verb|nPatch|\cdot \verb|nVars|\).
The time derivatives must be computed into the same sized array, but the patch edge values will be overwritten by zeros.
\item \verb|patches.x| is \(\verb|nSubP|\times \verb|nPatch|\) array of the spatial locations~\(x_{ij}\) of the microscale grid points in every patch.  
Currently it \emph{must} be a regular lattice on both macro- and micro-scales.
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is \(\verb|nSubP|\cdot \verb|nPatch|\cdot \verb|nVars|\) vector of time derivatives, but with zero on patch edges??
\end{itemize}

Try to figure out sizes of things.
An error in the reshape indicates~\verb|u| has the wrong length.
\begin{matlab}
%}
[nM,nP]=size(patches.x);
nV=round(numel(u)/numel(patches.x));
u=reshape(u,nM,nP,nV);
%{
\end{matlab}
With Dirichlet conditions on the patch edge, the half-length of a patch is \(h=dx(n_\mu-1)/2\) (or~\(-2\) for specified flux), and the ratio needed for interpolation is then \(r=h/\Delta X\).
Compute lattice sizes from inside the patches as the edge values may be \nan{}s, etc.
\begin{matlab}
%}
dx=patches.x(3,1)-patches.x(2,1);
DX=patches.x(2,2)-patches.x(2,1);
r=dx*(nM-1)/2/DX;
%{
\end{matlab}

For the moment?? assume the physical domain is macroscale periodic so that the coupling formulas are simplest.  
These index vectors point to patches and their two immediate neighbours.
\begin{matlab}
%}
j=1:nP; jp=mod(j,nP)+1; jm=mod(j-2,nP)+1;
%{
\end{matlab}
The centre of each patch, assuming odd~\verb|nM|, is at
\begin{matlab}
%}
i0=round((nM+1)/2);
%{
\end{matlab}
So compute centred differences of the mid-patch values for the macro-interpolation, of all fields.
For the moment?? assumes the domain is macro-periodic.
\begin{matlab}
%}
dmu=(u(i0,jp,:)-u(i0,jm,:))/2; % \mu\delta
ddu=(u(i0,jp,:)-2*u(i0,j,:)+u(i0,jm,:)); % \delta^2
dddmu=dmu(1,jp,:)-2*dmu(1,j,:)+dmu(1,jm,:); % \mu\delta^3
ddddu=ddu(1,jp,:)-2*ddu(1,j,:)+ddu(1,jm,:); % \delta^4
%dddddmu=dddmu(jp,:)-2*dddmu(j,:)+dddmu(jm,:); % \mu\delta^5
%ddddddu=ddddu(jp,:)-2*ddddu(j,:)+ddddu(jm,:); % \delta^6
%{
\end{matlab}
Interpolate macro-values to be Dirichlet edge values for each patch \cite[]{Roberts06d}.
Here interpolate to fourth order.
\begin{matlab}
%}
u(nM,j,:)=u(i0,j,:)+r*dmu+r^2/2*ddu ...
   +dddmu*(-1+r^2)*r/6+ddddu*(-1+r^2)*r^2/24;
u(1,j,:)=u(i0,j,:)-r*dmu+r^2/2*ddu ...
   -dddmu*(-1+r^2)*r/6+ddddu*(-1+r^2)*r^2/24;
%{
\end{matlab}
The following additions, respectively, is potentially sixth order.
\begin{verbatim}
   +dddddmu*(r/30-r^3/24+r^5/120) ...
   +ddddddu*(r^2/180-r^4/144+r^6/720);
   -dddddmu*(r/30-r^3/24+r^5/120) ...
   +ddddddu*(r^2/180-r^4/144+r^6/720);
\end{verbatim}

Ask the user for the time derivatives, overwrite edge values, the return as column vector.
\begin{matlab}
%}
dudt=patches.fun(t,u,patches.x);
dudt([1 nM],:,:)=0;
dudt=reshape(dudt,[],1);
%{
\end{matlab}
Fin.
%}
