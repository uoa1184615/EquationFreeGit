%Provides the coupling across space for patches of
%simulations of a smooth lattice system such as PDE
%discretisations.
%AJR, Nov 2017
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{patchSmooth1()}}
\label{sec:patchSmooth1}
\localtableofcontents

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
\item \verb|u| is a vector of length \(\verb|nSubP|\cdot \verb|nPatch|\cdot \verb|nVars|\) where there are \verb|nVars| field values at each of the points in the \(\verb|nSubP|\times \verb|nPatch|\) grid.
\item \verb|t| is the current time to be passed to the user's time derivative function.
\item \verb|patches| a struct set by \verb|makePatches()| with the following information.
\begin{itemize}
\item \verb|.fun| is the name of the user's function \verb|fun(t,u,x)| that computes the time derivatives on the patchy lattice. 
The array~\verb|u| has size \(\verb|nSubP|\times \verb|nPatch|\times \verb|nVars|\).
Time derivatives must be computed into the same sized array, but the patch edge values will be overwritten by zeros.
\item \verb|.x| is \(\verb|nSubP|\times \verb|nPatch|\) array of the spatial locations~\(x_{ij}\) of the microscale grid points in every patch.  
Currently it \emph{must} be a regular lattice on both macro- and micro-scales.
\item \verb|.ordCC|
\item \verb|.alt|
\item \verb|.Cwtsr| and \verb|.Cwtsl|
\end{itemize}

\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is \(\verb|nSubP|\cdot \verb|nPatch|\cdot \verb|nVars|\) vector of time derivatives, but with zero on patch edges??
\end{itemize}

Try to figure out sizes of things.
Any error arising in the reshape indicates~\verb|u| has the wrong length.
\begin{matlab}
%}
[nM,nP]=size(patches.x);
nV=round(numel(u)/numel(patches.x));
u=reshape(u,nM,nP,nV);
%{
\end{matlab}
With Dirichlet patches, the half-length of a patch is \(h=dx(n_\mu-1)/2\) (or~\(-2\) for specified flux), and the ratio needed for interpolation is then \(r=h/\Delta X\).
Compute lattice sizes from inside the patches as the edge values may be \nan{}s, etc.
\begin{matlab}
%}
dx=patches.x(3,1)-patches.x(2,1);
DX=patches.x(2,2)-patches.x(2,1);
r=dx*(nM-1)/2/DX;
%{
\end{matlab}

For the moment?? assume the physical domain is macroscale periodic so that the coupling formulas are simplest.  
Should cater for periodic, odd-mid-gap, even-mid-gap, even-mid-patch, dirichlet, neumann, ??
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
Assumes the domain is macro-periodic.
\begin{matlab}
%}
dmu=nan(patches.ordCC,nP,nV);
if patches.alt % use only odd numbered neighbours
    dmu(1,:,:)=(u(i0,jp,:)+u(i0,jm,:))/2; % \mu
    dmu(2,:,:)= u(i0,jp,:)-u(i0,jm,:); % \delta
    jp=jp(jp); jm=jm(jm); % increase shifts to \pm2
else % standard
    dmu(1,:,:)=(u(i0,jp,:)-u(i0,jm,:))/2; % \mu\delta
    dmu(2,:,:)=(u(i0,jp,:)-2*u(i0,j,:)+u(i0,jm,:)); % \delta^2
end% if
%{
\end{matlab}
Recursively take \(\delta^2\) of these to form higher order centred differences (could unroll a little to cater for two in parallel).
\begin{matlab}
%}
for k=3:patches.ordCC
    dmu(k,:,:)=dmu(k-2,jp,:)-2*dmu(k-2,j,:)+dmu(k-2,jm,:);
end
%{
\end{matlab}
Interpolate macro-values to be Dirichlet edge values for each patch \cite[]{Roberts06d}, using weights computed in \verb|makePatches()| .
Here interpolate to specified order.
\begin{matlab}
%}
u(nM,j,:)=u(i0,j,:)*(1-patches.alt) ...
    +sum(bsxfun(@times,patches.Cwtsr,dmu));
u( 1,j,:)=u(i0,j,:)*(1-patches.alt) ...
    +sum(bsxfun(@times,patches.Cwtsl,dmu));
%{
\end{matlab}

Ask the user for the time derivatives computed in the array, overwrite its edge values, then return to an integrator as column vector.
\begin{matlab}
%}
dudt=patches.fun(t,u,patches.x);
dudt([1 nM],:,:)=0;
dudt=reshape(dudt,[],1);
%{
\end{matlab}
Fin.
%}
