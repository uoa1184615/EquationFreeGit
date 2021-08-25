% mmPatchSys2() provides an interface to time integrators
% for the dynamics on moving patches coupled across space. The
% system must be a lattice system such as a PDE
% discretisation.  AJR, Aug 2021
%!TEX root = doc.tex
%{
\section{\texttt{mmPatchSys2()}: interface 2D space of
moving patches to time integrators}
\label{sec:mmPatchSys2}

\paragraph{Beware ad hoc assumptions}
In an effort to get started, I have just made some plausible generalisations from the 1D code to this 2D code.
Probably lots of details are poor??

To simulate in time with 2D patches moving in space we need
to interface a users time derivative function with time
integration routines such as \verb|ode23| or~\verb|PIRK2|.
This function \verb|mmPatchSys2()| provides an interface.
Patch edge values are determined by macroscale interpolation
of the patch-centre or edge values. Microscale heterogeneous
systems may be accurately simulated with this function via
appropriate interpolation. Communicate patch-design
variables (\cref{sec:configPatches2}) either via the global
struct~\verb|patches| or via an optional third argument
(except that this last is required for parallel computing of
\verb|spmd|).

\begin{matlab}
%}
function dudt = mmPatchSys2(t,u,patches)
if nargin<3, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector of length
$2\cdot\verb|prod(nPatch)| +\verb|prod(nSubP)| \cdot
\verb|nVars| \cdot \verb|nEnsem| \cdot \verb|prod(nPatch)|$
where there are $\verb|nVars| \cdot \verb|nEnsem|$ field
values at each of the points in the $\verb|nSubP(1)| \times
\verb|nSubP(2)| \times \verb|nPatch(1)| \times
\verb|nPatch(2)|$ grid.

\item \verb|t| is the current time to be passed to the
user's time derivative function.

\item \verb|patches| a struct set by \verb|configPatches2()|
with the following information used here.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,M,patches)| that computes the time derivatives
on the patchy lattice, where the \((I,J)\)th~patch moves at
velocity~\((M.Vx_I,M.Vy_J)\) and at current time is
displaced~\((M.Dx_I,M.Dy_J)\) from the fixed reference
positions in~\verb|.x| and~\verb|.y|\,.  The array~\verb|u|
has size $\verb|nSubP(1)| \times \verb|nSubP(2)| \times
\verb|nVars| \times \verb|nEsem| \times \verb|nPatch(1)|
\times \verb|nPatch(2)|$.  Time derivatives must be computed
into the same sized array, although herein the patch
edge-values are overwritten by zeros.

\item \verb|.x| is $\verb|nSubP(1)| \times1 \times1 \times1
\verb|nPatch(1)| \times1$ array of the spatial
locations~$x_{i}$ of the microscale $(i,j)$-grid points in
every patch.  Currently it \emph{must} be an equi-spaced
lattice on both macro- and micro-scales??

\item \verb|.y| is similarly $1 \times \verb|nSubP(2)|
\times1 \times1 \times1 \times \verb|nPatch(2)|$ array of
the spatial locations~$y_{j}$ of the microscale $(i,j)$-grid
points in every patch.  Currently it \emph{must} be an
equi-spaced lattice on both macro- and micro-scales.

\item \verb|.Xlim| ??

\end{itemize}
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is a vector\slash array of of time
derivatives, but with patch edge-values set to zero. It is
of total length $2\cdot\verb|prod(nPatch)| + \verb|prod(nSubP)| \cdot \verb|nVars|
\cdot \verb|nEnsem| \cdot \verb|prod(nPatch)|$ and the same
dimensions as~\verb|u|.
\end{itemize}



\begin{devMan}
Extract the \(2\cdot\verb|prod(nPatch)|\) displacement
values from the start of the vectors of evolving variables.
Reshape the rest as the fields~\verb|u| in a 6D-array, and
sets the edge values from macroscale interpolation of
centre-patch values. \cref{sec:patchEdgeInt2} describes
\verb|patchEdgeInt2()|.
\begin{matlab}
%}
Nx = size(patches.x,5);
Ny = size(patches.y,6);
nM = Nx*Ny;
M.Dx = reshape(u(   1:nM  ),[1 1 1 1 Nx Ny]); 
M.Dy = reshape(u(nM+1:2*nM),[1 1 1 1 Nx Ny]);
u = patchEdgeInt2(u(2*nM+1:end),patches);
%{
\end{matlab}



\paragraph{Moving mesh velocity}
Developing from standard moving meshes for \pde{}s
\cite[e.g.]{Budd2009, Huang10}, we follow
\cite{Maclean2021a}, and generalise ad hoc to~2D??  There
exists a set of macro-scale mesh points $\big(X_{IJ}(t) \C
Y_{IJ}(t) \big) := \big( X_{IJ}^0+Dx_{IJ}(t) \C
Y_{IJ}^0+Dy_{IJ}(t) \big)$ (at the centre) of each patch
with associated field values, say
$U_{IJ}(t):=\overline{u_{ijIJ}(t)}$. And remove the two
microscale dimensions from the front of the arrays, so they
are 4D arrays.
\begin{matlab}
%}
X = shiftdim( mean(patches.x,1)+M.Dx ,2);
Y = shiftdim( mean(patches.y,2)+M.Dy ,2);
U = shiftdim( mean(mean(u,1,'omitnan'),2,'omitnan') ,2);
%Uz=squeeze(U)
%{
\end{matlab}
Then for every patch~\((I,J)\) we set \(H^{pq}_{IJ}:={}\)the
\(q\)th~spatial component of the step to the next patch in
the \(p\)th~index direction, for periodic patch
indices~\((I,J)\),
\begin{matlab}
%}
I=1:Nx; Ip=[2:Nx 1]; Im=[Nx 1:Nx-1];
J=1:Ny; Jp=[2:Ny 1]; Jm=[Ny 1:Ny-1];
Hix = X(:,:,Ip,J)-X(:,:,I,J);
Hiy = Y(:,:,Ip,J)-Y(:,:,I,J);
Hjx = X(:,:,I,Jp)-X(:,:,I,J);
Hjy = Y(:,:,I,Jp)-Y(:,:,I,J);
Hix(:,:,Nx,:) = Hix(:,:,Nx,:)+diff(patches.Xlim(1:2));
Hjy(:,:,:,Ny) = Hjy(:,:,:,Ny)+diff(patches.Xlim(3:4));
%{
\end{matlab}
we discretise a moving mesh \pde\ for node
locations~$(X_{IJ},Y_{IJ})$ with field values~$U_{IJ}$ via
the second derivatives estimates ??
\begin{subequations} \label{mm2Disc}
\begin{equation}
U''_{j} := \frac2{H_j+H_{j-1}}\left[ \frac{U_{j+1} -
U_j}{H_j} - \frac{U_{j} - U_{j-1}}{H_{j-1}} \right] .
\end{equation}
First, compute first derivatives at \((I+\tfrac12,J)\) and
\((I,J+\tfrac12)\) respectively.
\begin{matlab}
%}
Ux = (U(:,:,Ip,J)-U(:,:,I,J))./Hix(:,:,I,J); %ux=squeeze(Ux)
Uy = (U(:,:,I,Jp)-U(:,:,I,J))./Hjy(:,:,I,J); %uy=squeeze(Uy)
%{
\end{matlab}
Second, compute second derivative matrix, without assuming
symmetry because the derivatives in space are not quite the
same as the derivatives in indices.  The mixed derivatives
are at \((I+\tfrac12,J+\tfrac12)\), so average to get at
patch locations.
\begin{matlab}
%}
Uxx = ( Ux(:,:,I,J)-Ux(:,:,Im,J) )*2./(Hix(:,:,I,J)+Hix(:,:,Im,J));
Uyy = ( Uy(:,:,I,J)-Uy(:,:,I,Jm) )*2./(Hjy(:,:,I,J)+Hjy(:,:,I,Jm));
Uyx = ( Uy(:,:,Ip,J)-Uy(:,:,I,J) )./Hix(:,:,I,J);
Uxy = ( Ux(:,:,I,Jp)-Ux(:,:,I,J) )./Hjy(:,:,I,J);
Uyx = (Uyx(:,:,I,J)+Uyx(:,:,Im,J)+Uyx(:,:,I,Jm)+Uyx(:,:,Im,Jm))/4;
Uxy = (Uxy(:,:,I,J)+Uxy(:,:,Im,J)+Uxy(:,:,I,Jm)+Uxy(:,:,Im,Jm))/4;
%uxx=squeeze(Uxx),uyy=squeeze(Uyy),uxy=squeeze(Uxy),uyx=squeeze(Uyx),
%{
\end{matlab}
And compute its norm over all variables and ensembles
(arbitrarily?? chose the mean square norm here, using
\verb|abs.^2| as they may be complex), shifting the variable
and ensemble dimensions out of the result to give 2D array
of values, one for each patch (use \verb|shiftdim| rather
than \verb|squeeze| as users may invoke a 1D array of 2D
patches, as in channel dispersion).
\begin{matlab}
%}
U2 = shiftdim( mean(mean( ...
         abs(Uxx).^2+abs(Uyy).^2+abs(Uxy).^2+abs(Uyx).^2 ...
      ,1),2) ,2);
Hix = shiftdim(Hix,2);  Hiy = shiftdim(Hiy,2);
Hjx = shiftdim(Hjx,2);  Hjy = shiftdim(Hjy,2);
%{
\end{matlab}
Having squeezed out all microscale information, the
global moderating coefficient in~1D??
\begin{equation}
\label{mm2Al}
\alpha := \max\left\{ 1\C \left[\frac{1}{b-a}\sum_{j} 
H_{j-1} \frac12 \left({U''_{j}}^{2/3} +
{U''_{j-1}}^{2/3}\right)\right]^3 \right\}
\end{equation}
generalises to an integral over \emph{approximate}
parallelograms in~2D?? (area approximately?? determined by
cross-product). Rather than \(\max(1,\cdot)\) surely better
to use something smooth like~\(\sqrt(1+\cdot^2)\)??
\begin{matlab}
%}
U23 = U2.^(1/3);
alpha = sum(sum( ...
    abs( Hix(Im,Jm).*Hjy(Im,Jm)-Hiy(Im,Jm).*Hjx(Im,Jm) ) ...
    .*( U23(I,J)+U23(Im,J)+U23(I,Jm)+U23(Im,Jm) )/4 ...
    ))/diff(patches.Xlim(1:2))/diff(patches.Xlim(3:4));
alpha = sqrt(1+alpha^6);
%{
\end{matlab}
Then the importance function at each patch is the 2D array 
\begin{equation}
\label{mm2R}
\rho_j := \left(1 + \frac{1}{\alpha} {U''_{j}}^2 \right)^{1/3},
\end{equation}
\begin{matlab}
%}
rho = ( 1+U2/alpha ).^(1/3);
%{
\end{matlab}
For every patch, we move all micro-grid points according to
the following velocity of the notional macro-scale node of
that patch: (Since we differentiate the importance function,
maybe best to compute it above at half-grid points of the
patches---aka a staggered scheme??)
\begin{equation} 
\label{mm2X}
V_j:= \de t{X_j} = \frac{(N-1)^2}{2\rho_j \tau } \left[
(\rho_{j+1}+\rho_j) H_j - (\rho_j + \rho_{j-1}) H_{j-1}
\right]. 
\end{equation}
Is the \verb|Nx| and \verb|Ny| correct here??
And are the derivatives appropriate since these here are 
scaled index derivatives, not actually spatial derivatives??
\begin{matlab}
%}
M.Vx = nan+M.Dx;  M.Vy = nan+M.Dy;  % allocate storage
M.Vx(:) = ( (rho(Ip,J)+rho(I,J)).*Hix(I,J) ...
           -(rho(Im,J)+rho(I,J)).*Hix(Im,J) ) ...
          ./rho(I,J) *(Nx^2/2/patches.mmTime);
M.Vy(:) = ( (rho(I,Jp)+rho(I,J)).*Hjy(I,J) ...
           -(rho(I,Jm)+rho(I,J)).*Hjy(I,Jm) ) ...
          ./rho(I,J) *(Ny^2/2/patches.mmTime);
%Vx=squeeze(M.Vx), Vy=squeeze(M.Vy), return
%{
\end{matlab}
\end{subequations}




\paragraph{Evaluate system differential equation}
Ask the user function for the time derivatives computed in
the array, overwrite its edge values with the dummy value of
zero (as \verb|ode15s| chokes on NaNs), then return to the
user\slash integrator as same sized array as input.
\begin{matlab}
%}
dudt = patches.fun(t,u,M,patches);
dudt([1 end],:,:,:,:,:) = 0;
dudt(:,[1 end],:,:,:,:) = 0;
dudt=[M.Vx(:); M.Vy(:); dudt(:)];
%{
\end{matlab}
Fin.
\end{devMan}
%}
 