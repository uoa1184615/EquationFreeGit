% mmPatchSys1() provides an interface to time integrators
% for the dynamics on moving patches coupled across space. The
% system must be a lattice system such as a PDE
% discretisation.  AJR, Aug 2021
%!TEX root = doc.tex
%{
\section{\texttt{mmPatchSys1()}: interface 1D space of
moving patches to time integrators}
\label{sec:mmPatchSys1}


To simulate in time with moving 1D spatial patches we need
to interface a user's time derivative function with time
integration routines such as \verb|ode23| or~\verb|PIRK2|.
This function \verb|mmPatchSys1()| provides an interface.
Patch edge values are determined by macroscale
interpolation of the patch-centre or edge values. 
Microscale heterogeneous systems may be accurately simulated
with this function via appropriate interpolation. 
Communicate patch-design variables
(\cref{sec:configPatches1}) either via the global
struct~\verb|patches| or via an optional third argument
(except that this last is required for parallel computing of
\verb|spmd|).

\begin{matlab}
%}
function dudt = mmPatchSys1(t,u,patches)
if nargin<3, global patches, end
%{
\end{matlab}


\paragraph{Input}
\begin{itemize}

\item \verb|u| is a vector of length $\verb|nPatch| +
\verb|nSubP| \cdot \verb|nVars| \cdot \verb|nEnsem| \cdot
\verb|nPatch|$  where there are $\verb|nVars| \cdot
\verb|nEnsem|$ field values at each of the points in the
$\verb|nSubP|\times \verb|nPatch|$ grid, and because of the
moving mesh there are an additional \verb|nPatch| patch
displacement values at its start.

\item \verb|t| is the current time to be passed to the
user's time derivative function.

\item \verb|patches| a struct set by \verb|configPatches1()|
with the following information  used here.
\begin{itemize}

\item \verb|.fun| is the name of the user's function
\verb|fun(t,u,M,patches)| that computes the time derivatives
on the patchy lattice, where the \(j\)th~patch moves at
velocity~\(M.V_j\) and at current time is
displaced~\(M.D_j\) from the fixed reference position
in~\verb|.x|\,.  The array~\verb|u| has size $\verb|nSubP|
\times \verb|nVars| \times \verb|nEnsem| \times
\verb|nPatch|$.  Time derivatives should be computed into
the same sized array, then herein the patch edge values are
overwritten by zeros.

\item \verb|.x| is $\verb|nSubP| \times1 \times1 \times
\verb|nPatch|$ array of the spatial locations~$x_{i}$ of
the microscale grid points in every patch.  Currently it
\emph{must} be an equi-spaced lattice on both macro- and
microscales ??
\end{itemize}
\end{itemize}


\paragraph{Output}
\begin{itemize}
\item \verb|dudt| is a vector of of time
derivatives, but with patch edge-values set to zero.  It is
of total length $\verb|nPatch| + \verb|nSubP| \cdot
\verb|nVars| \cdot \verb|nEnsem| \cdot \verb|nPatch|$.
\end{itemize}



\begin{devMan}
Extract the \verb|nPatch| displacement values from the start
of the vectors of evolving variables. Reshape the rest as
the fields~\verb|u| in a 4D-array, and sets the edge values
from macroscale interpolation of centre-patch values.
\cref{sec:patchEdgeInt1} describes \verb|patchEdgeInt1()|.
\begin{matlab}
%}
N = size(patches.x,4);
M.D = reshape(u(1:N),[1 1 1 N]);
u = patchEdgeInt1(u(N+1:end),patches);
%{
\end{matlab}


\paragraph{Moving mesh velocity}
Developing from standard moving meshes for \pde{}s
\cite[e.g.]{Budd2009, Huang10}, we follow
\cite{Maclean2021a}. There exists a set of macro-scale mesh
points~$X_j(t):=X_j^0+D_j(t)$ (at the centre) of each patch
with associated field values, say
$U_j(t):=\overline{u_{ij}(t)}$.
\begin{matlab}
%}
X = mean(patches.x,1)+M.D;
U = mean(u,1);
%{
\end{matlab}
Then for every patch~\(j\) we set \(H_j:=X_{j+1}-X_j\) for
periodic patch indices~\(j\)
\begin{matlab}
%}
j=1:N; jp=[2:N 1]; jm=[N 1:N-1];
H = X(:,:,:,jp)-X(:,:,:,j);  
H(N) = H(N)+diff(patches.Xlim);
%{
\end{matlab}
we discretise a moving mesh \pde\ for node locations~$X_j$
with field values~$U_j$ via the second derivative estimate
\begin{subequations} \label{mmDisc}
\begin{equation}
U''_{j} := \frac2{H_j+H_{j-1}}\left[ \frac{U_{j+1} -
U_j}{H_j} - \frac{U_{j} - U_{j-1}}{H_{j-1}} \right] ,
\end{equation}
and its norm over all variables and ensembles (arbitrarily??
chose the  mean square norm here).
\begin{matlab}
%}
U2 = ( (U(:,:,:,jp)-U(:,:,:,j))./H(:,:,:,j) ...
      -(U(:,:,:,j)-U(:,:,:,jm))./H(:,:,:,jm) ...
     )*2./(H(:,:,:,j)+H(:,:,:,jm));
U2 = squeeze( mean(mean( abs(U2).^2 ,2),3) );
H = squeeze(H);
%{
\end{matlab}
Having squeezed out all microscale information, the
coefficient
\begin{equation}
\label{mmAl}
\alpha := \max\left\{ 1\C \left[\frac{1}{b-a}\sum_{j} 
H_{j-1} \frac12 \left({U''_{j}}^{2/3} +
{U''_{j-1}}^{2/3}\right)\right]^3 \right\}
\end{equation}
Rather than \(\max(1,\cdot)\) surely better to use something 
smooth like~\(\tanh(\cdot)\)??
\begin{matlab}
%}
alpha = sum( H(jm).*( U2(j).^(1/3)+U2(jm).^(1/3) ))/2/sum(H);
alpha = max(1,alpha^3);
%{
\end{matlab}
Then the importance function
\begin{equation}
\label{mmR}
\rho_j := \left(1 + \frac{1}{\alpha} {U''_{j}}^2 \right)^{1/3},
\end{equation}
\begin{matlab}
%}
rho = ( 1+U2/alpha ).^(1/3);
%{
\end{matlab}
For every patch, we move all micro-grid points according to
the following velocity of the notional macro-scale node of
that patch:
\begin{equation} 
\label{mmX}
V_j:= \de t{X_j} = \frac{(N-1)^2}{2\rho_j \tau } \left[
(\rho_{j+1}+\rho_j) H_j - (\rho_j + \rho_{j-1}) H_{j-1}
\right]. 
\end{equation}
\begin{matlab}
%}
M.V = nan+M.D; % allocate storage
M.V(:) = ( (rho(jp)+rho(j)).*H(j) -(rho(j)+rho(jm)).*H(jm) ) ...
    ./rho(j) *((N-1)^2/2/patches.mmTime);
%{
\end{matlab}
\end{subequations}


\paragraph{Evaluate system differential equation}
Ask the user function for the advected time derivatives on
the moving patches, overwrite its edge values with the dummy
value of zero (since \verb|ode15s| chokes on NaNs), then
return to the user\slash integrator as a vector.
\begin{matlab}
%}
dudt=patches.fun(t,u,M,patches);
dudt([1 end],:,:,:) = 0;
dudt=[M.V(:); dudt(:)];
%{
\end{matlab}
Fin.
\end{devMan}
%}
 