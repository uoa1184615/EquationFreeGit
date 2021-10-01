% mmPatchSys2() provides an interface to time integrators
% for the dynamics on moving patches coupled across space. The
% system must be a lattice system such as a PDE
% discretisation.  AJR, Aug 2021 -- Sep 2021
%!TEX root = doc.tex
%{
\section{\texttt{mmPatchSys2()}: interface 2D space of
moving patches to time integrators}
\label{sec:mmPatchSys2}

\paragraph{Beware ad hoc assumptions}
In an effort to get started, I make some plausible generalisations from the 1D code to this 2D code, in the option \verb|adhoc|.
Also, I code the alternative \verb|Huang98| which aims to implement the method of \cite{Huang98}.

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
global theMethod % adhoc or Huang98
global ind % =1 for x-dirn and =2 for y-dirn testing Huang
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
\cite{Maclean2021a}, and generalise to~2D according to the algorithm of \cite{Huang98}, and also ad hoc.  
Here the patch indices~\(I,J\) play the role of mesh variables~\(\xi,\eta\) of \cite{Huang98}.  There
exists a set of macro-scale mesh points $\Xv_{IJ}(t) := \big(X_{IJ}(t) \C Y_{IJ}(t) \big) := \big( X_{IJ}^0+Dx_{IJ}(t) \C Y_{IJ}^0+Dy_{IJ}(t) \big)$ (at the centre) of each patch
with associated field values, say
$U_{IJ}(t):=\overline{u_{ijIJ}(t)}$.  
Also, remove microscale dimensions from the front of these macro-mesh arrays, so \verb|X,Y| are~2D, and \verb|U| is~4D.
\begin{matlab}
%}
X = shiftdim( mean(patches.x,1)+M.Dx ,4);
Y = shiftdim( mean(patches.y,2)+M.Dy ,4);
U = shiftdim( mean(mean(u,1,'omitnan'),2,'omitnan') ,2);
%{
\end{matlab}
Then for every patch~\((I,J)\) we set \(??:={}\)the
\(q\)th~spatial component of the step to the next patch in
the \(p\)th~index direction, for periodic patch
indices~\((I,J)\).
Throughout, use appended \verb|r,u| to denote mesh-midpoint quantities at \(I+\frac12\) and \(J+\frac12\), respectively, and use \verb|_I,_J| to respectively denote differences in the macro-mesh indices~\(I,J\) which then estimate derivatives in the mesh parameters, \(\D\xi{}\) and \(\D\eta{}\), respectively.  
\begin{matlab}
%}
I=1:Nx; Ip=[2:Nx 1]; Im=[Nx 1:Nx-1];
J=1:Ny; Jp=[2:Ny 1]; Jm=[Ny 1:Ny-1];
Xr_I = X(Ip,J)-X(I,J); % propto dX/dxi
Yr_I = Y(Ip,J)-Y(I,J); % propto dY/dxi
Xu_J = X(I,Jp)-X(I,J); % propto dX/deta
Yu_J = Y(I,Jp)-Y(I,J); % propto dY/deta
Xr_I(Nx,:) = Xr_I(Nx,:)+diff(patches.Xlim(1:2));
Yu_J(:,Ny) = Yu_J(:,Ny)+diff(patches.Xlim(3:4));
%{
\end{matlab}



\subsection{ad hoc attempt}
\begin{matlab}
%}
switch theMethod
case 'adhoc'
%{
\end{matlab}
Temporarily shift the macro-mesh info into dimensions~3 and~4:
\begin{matlab}
%}
Xr_I = shiftdim(Xr_I,-2);  Yr_I = shiftdim(Yr_I,-2);
Xu_J = shiftdim(Xu_J,-2);  Yu_J = shiftdim(Yu_J,-2);
%{
\end{matlab}
We discretise a moving mesh \pde\ for node
locations~$(X_{IJ},Y_{IJ})$ with field values~$U_{IJ}$ via
the second derivatives estimates ??
\begin{subequations} \label{mm2Disc}
\begin{equation}
U''_{j} := \frac2{H_j+H_{j-1}}\left[ \frac{U_{j+1} -
U_j}{H_j} - \frac{U_{j} - U_{j-1}}{H_{j-1}} \right] .
\end{equation}
First, compute first derivatives at \((I+\tfrac12,J)\) and
\((I,J+\tfrac12)\) respectively---these are derivatives in the macro-mesh directions, incorrectly scaled.
\begin{matlab}
%}
Ux = (U(:,:,Ip,J)-U(:,:,I,J))./Xr_I(:,:,I,J); 
Uy = (U(:,:,I,Jp)-U(:,:,I,J))./Yu_J(:,:,I,J); 
%{
\end{matlab}
Second, compute second derivative matrix, without assuming
symmetry because the derivatives in space are not quite the
same as the derivatives in indices.  The mixed derivatives
are at \((I+\tfrac12,J+\tfrac12)\), so average to get at
patch locations.
\begin{matlab}
%}
Uxx = ( Ux(:,:,I,J)-Ux(:,:,Im,J) )*2./(Xr_I(:,:,I,J)+Xr_I(:,:,Im,J));
Uyy = ( Uy(:,:,I,J)-Uy(:,:,I,Jm) )*2./(Yu_J(:,:,I,J)+Yu_J(:,:,I,Jm));
Uyx = ( Uy(:,:,Ip,J)-Uy(:,:,I,J) )./Xr_I(:,:,I,J);
Uxy = ( Ux(:,:,I,Jp)-Ux(:,:,I,J) )./Yu_J(:,:,I,J);
Uyx = (Uyx(:,:,I,J)+Uyx(:,:,Im,J)+Uyx(:,:,I,Jm)+Uyx(:,:,Im,Jm))/4;
Uxy = (Uxy(:,:,I,J)+Uxy(:,:,Im,J)+Uxy(:,:,I,Jm)+Uxy(:,:,Im,Jm))/4;
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
Xr_I = shiftdim(Xr_I,2);  Yr_I = shiftdim(Yr_I,2);
Xu_J = shiftdim(Xu_J,2);  Yu_J = shiftdim(Yu_J,2);
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
    abs( Xr_I(Im,Jm).*Yu_J(Im,Jm)-Yr_I(Im,Jm).*Xu_J(Im,Jm) ) ...
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
M.Vx = shiftdim( ...
    ((rho(Ip,J)+rho(I,J)).*Xr_I(I,J) ...
    -(rho(Im,J)+rho(I,J)).*Xr_I(Im,J) ) ...
    ./rho(I,J) *(Nx^2/2/patches.mmTime) ...
    ,-4);
M.Vy = shiftdim( ...
    ((rho(I,Jp)+rho(I,J)).*Yu_J(I,J) ...
    -(rho(I,Jm)+rho(I,J)).*Yu_J(I,Jm) ) ...
    ./rho(I,J) *(Ny^2/2/patches.mmTime) ...
    ,-4);
%{
\end{matlab}
\end{subequations}









\subsection{Huang98}
Here encode the algorithm of \cite{Huang98}.
\begin{matlab}
%}
case 'Huang98' %= theMethod
%{
\end{matlab}

The Jacobian at the \(N_x\times N_y\) mesh-points is, using centred differences,
\begin{matlab}
%}
Jac = 0.25*( (Xr_I+Xr_I(Im,J)).*(Yu_J+Yu_J(I,Jm)) ...
            -(Yr_I+Yr_I(Im,J)).*(Xu_J+Xu_J(I,Jm)) );
%{
\end{matlab}

\begin{subequations} \label{mm2Disc}
The mesh movement \pde\ is \cite[(24), with \(\gamma_2=0\)]{Huang98}
\begin{align}&
\D t\Xv =-\frac{\Xv_\xi}{\tau\sqrt{\tilde g_1}\cJ}\left\{
+\D\xi{}\left[\frac{\Xv_\eta^TG_1\Xv_\eta}{\cJ g_1}\right]
-\D\eta{}\left[\frac{\Xv_\xi^TG_1\Xv_\eta}{\cJ g_1}\right]
\right\}
\nonumber\\&
\phantom{\D t\Xv =}{}
-\frac{\Xv_\eta}{\tau\sqrt{\tilde g_2}\cJ}\left\{
-\D\xi{}\left[\frac{\Xv_\eta^TG_2\Xv_\xi}{\cJ g_2}\right]
+\D\eta{}\left[\frac{\Xv_\xi^TG_2\Xv_\xi}{\cJ g_2}\right]
\right\} ,
\label{mm2DiscdXdt}
\\&
\text{Jacobian }\cJ :=X_\xi Y_\eta -X_\eta Y_\xi\,,
\\&
g_k :=\det(G_k),
\\&
G_1 :=\sqrt{1+\|\grad U\|^2}\left[(1-\gamma_1)\cI 
+\gamma_1S(\grad\tilde \xi)\right],
\\&
G_2 :=\sqrt{1+\|\grad U\|^2}\left[(1-\gamma_1)\cI 
+\gamma_1S(\grad\tilde\eta)\right],
\\&
\text{matrix }S(\vv) :=\vv_\perp\vv_\perp^T/\|\vv\|^2
\quad\text{for }\vv_\perp:=(v_2,-v_1),
\\&\text{identity }\cI.
\end{align}
In their examples, \cite{Huang98} chose the mesh orthogonality parameter \(\gamma_1=0.1\)\,.
\begin{matlab}
%}
gamma1=0.1;
%{
\end{matlab}
The tildes appear to denote a reference mesh \cite[p.1005]{Huang98} which could be the identity map \((\tilde\xi,\tilde\eta)=(x,y)\), so here maybe \((\tilde I,\tilde J)=(\tilde\xi,\tilde\eta)=(X/H_x,Y/H_y)\).

We discretise the moving mesh \pde\ for node
locations~$(X_{IJ},Y_{IJ})$ with field values~$U_{IJ}$ via
the second derivatives estimates ??
So \cite{Huang98}'s \(\xi,i\mapsto I\), and \(\eta,j\mapsto J\), and \(\xv\mapsto \verb|Xv|\)??

\paragraph{Importance functions}
First, compute the gradients of the macroscale field formed into \(w=\sqrt{1+\|\grad\Uv\|^2}\) \cite[(30)]{Huang98}, using centred differences from patch to patch, unless we use the patches to estimate first derivatives (implicitly the interpolation).
Need to shift dimensions of macroscale mesh to cater for components of the field~\Uv.
\begin{matlab}
%}
Y_J = shiftdim( (Yu_J(I,J)+Yu_J(I,Jm))/2 ,-2);
Y_I = shiftdim( (Yr_I(I,J)+Yr_I(Im,J))/2 ,-2);
U_x = ( (Y_J(:,:,Ip,J).*U(:,:,Ip,J)-Y_J(:,:,Im,J).*U(:,:,Im,J))/2 ...
       -(Y_I(:,:,I,Jp).*U(:,:,I,Jp)-Y_I(:,:,I,Jm).*U(:,:,I,Jm))/2 ...
      )./Jac;
X_J = shiftdim( (Xu_J(I,J)+Xu_J(I,Jm))/2 ,-2);
X_I = shiftdim( (Xr_I(I,J)+Xr_I(Im,J))/2 ,-2);
U_y = ( (X_J(:,:,Ip,J).*U(:,:,Ip,J)-X_J(:,:,Im,J).*U(:,:,Im,J))/2 ...
       -(X_I(:,:,I,Jp).*U(:,:,I,Jp)-X_I(:,:,I,Jm).*U(:,:,I,Jm))/2 ...
      )./Jac;
w = sqrt(1 + sum(sum(U_x.^2+U_y.^2,1),2) );
testy(w,2+ind,'w')
%{
\end{matlab}
In order to compute~\(G_k\), it seems \(\grad\eta=(Y_\eta,-Y_\xi)/\cJ\) and \(\grad\xi=(-X_\eta,X_\xi)/\cJ\).   
Then, \(S(\grad\xi)\) has \(\vv_\perp=(X_\xi,X_\eta)/\cJ\) so \(S(\grad\xi)=\begin{bmat} X_\xi^2 &X_\xi X_\eta \\ X_\xi X_\eta &X_\eta^2 \end{bmat}/(X_\xi^2+X_\eta^2)\).
Now \cite{Huang98} has tildes on these, so they are meant to be reference coordinates?? in which case we would have \(X_\eta=0\)\,, so \(S=\begin{bmat}1 &0 \\ 0 &0 \end{bmat}\), and so \(G_1\propto \begin{bmat} 1&0\\0&1-\gamma_1 \end{bmat}\)---\emph{I do not see how this helps stop the mesh degenerating}.
\begin{matlab}
%}
G1 = w.*( (1-gamma1)*eye(2) ...
         +gamma1*[Y_I.^2 Y_I.*Y_J; Y_I.*Y_J Y_J.^2]./(Y_I.^2+Y_J.^2) );
G2 = w.*( (1-gamma1)*eye(2) ...
         +gamma1*[X_I.^2 X_I.*X_J; X_I.*X_J X_J.^2]./(X_I.^2+X_J.^2) );
testy(G1,2+ind,'G1')
testy(G2,2+ind,'G2')
%{
\end{matlab}

Apply low-pass filter \cite[(27)]{Huang98} (although unclear whether to apply the filter four times to the whole of both matrices?? or once to each of the four components of both matrices??):
\begin{matlab}
%}
for k=1:1
G1 = G1/4 ...
  +(G1(:,:,Ip,J)+G1(:,:,Im,J)+G1(:,:,I,Jp)+G1(:,:,I,Jm))/8 ...
  +(G1(:,:,Ip,Jp)+G1(:,:,Im,Jm)+G1(:,:,Im,Jp)+G1(:,:,Ip,Jm))/16;
G2 = G2/4 ...
  +(G2(:,:,Ip,J)+G2(:,:,Im,J)+G2(:,:,I,Jp)+G2(:,:,I,Jm))/8 ...
  +(G2(:,:,Ip,Jp)+G2(:,:,Im,Jm)+G2(:,:,Im,Jp)+G2(:,:,Ip,Jm))/16;
end
testy(G1,2+ind,'G1')
testy(G2,2+ind,'G2')
%{
\end{matlab}


\paragraph{Macro-mesh movement} 
These are the \(2\times2\) matrices at \(N_x\times N_y\) midpoints of the mesh-net \cite[(27)]{Huang98}:
\begin{matlab}
%}
G1r = (G1(:,:,Ip,:)+G1(:,:,I,:))/2;
G2r = (G2(:,:,Ip,:)+G2(:,:,I,:))/2;
G1u = (G1(:,:,:,Jp)+G1(:,:,:,J))/2;
G2u = (G2(:,:,:,Jp)+G2(:,:,:,J))/2;
testy(G1r,2+ind,'G1')
testy(G2r,2+ind,'G2')
testy(G1u,2+ind,'G1')
testy(G2u,2+ind,'G2')
%{
\end{matlab}
Compute \(N_x\times N_y\) determinant of matrices \cite[(27)]{Huang98}:
\begin{matlab}
%}
g1 = shiftdim( G1(1,1,:,:).*G1(2,2,:,:) ...
               -G1(1,2,:,:).*G1(2,1,:,:) ,2);
g2 = shiftdim( G2(1,1,:,:).*G2(2,2,:,:) ...
               -G2(1,2,:,:).*G2(2,1,:,:) ,2);
testy(g1,ind,'g1')
testy(g2,ind,'g2')
g1r = shiftdim( G1r(1,1,:,:).*G1r(2,2,:,:) ...
               -G1r(1,2,:,:).*G1r(2,1,:,:) ,2);
g2r = shiftdim( G2r(1,1,:,:).*G2r(2,2,:,:) ...
               -G2r(1,2,:,:).*G2r(2,1,:,:) ,2);
g1u = shiftdim( G1u(1,1,:,:).*G1u(2,2,:,:) ...
               -G1u(1,2,:,:).*G1u(2,1,:,:) ,2);
g2u = shiftdim( G2u(1,1,:,:).*G2u(2,2,:,:) ...
               -G2u(1,2,:,:).*G2u(2,1,:,:) ,2);
%{
\end{matlab}
Compute vector \verb|Xv._.| of coordinate derivatives \cite[(27)]{Huang98}---use arrays \verb|X_.| and \verb|Y_.| here as they know the macro-periodicity.  
\begin{matlab}
%}
Xr_J = 0.25*(Xu_J(I,Jm)+Xu_J(I,J)+Xu_J(Ip,Jm)+Xu_J(Ip,J));
Yr_J = 0.25*(Yu_J(I,Jm)+Yu_J(I,J)+Yu_J(Ip,Jm)+Yu_J(Ip,J));
Xvr_J = [shiftdim(Xr_J,-1);shiftdim(Yr_J,-1)];
Xvr_I = [shiftdim(Xr_I,-1);shiftdim(Yr_I,-1)];
testy(Xvr_J,1+ind,'Xvr_J')
testy(Xvr_I,1+ind,'Xvr_I')
Xu_I = 0.25*(Xr_I(Im,J)+Xr_I(I,J)+Xr_I(Im,Jp)+Xr_I(I,Jp));
Yu_I = 0.25*(Yr_I(Im,J)+Yr_I(I,J)+Yr_I(Im,Jp)+Yr_I(I,Jp));
Xvu_I = [shiftdim(Xu_I,-1);shiftdim(Yu_I,-1)];
Xvu_J = [shiftdim(Xu_J,-1);shiftdim(Yu_J,-1)];
testy(Xvu_J,1+ind,'Xvu_J')
testy(Xvu_I,1+ind,'Xvu_I')
%{
\end{matlab}
Then the two Jacobians at the \(N_x\times N_y\) midpoints of the mesh-net are \cite[(27)]{Huang98}:
\begin{matlab}
%}
Jacr = Xr_I.*Yr_J - Yr_I.*Xr_J;
Jacu = Yu_J.*Xu_I - Xu_J.*Yu_I;
testy(Jacr,ind,'Jacr')
testy(Jacu,ind,'Jacu')
%{
\end{matlab}

For vectors~\xv,\yv\ of dimension \(d\times N_x\times N_y\) and array~\(G\) of dimension \(d\times d\times N_x\times N_y\), define function to evaluate product~\(\xv^TG\yv\) of dimension \(N_x\times N_y\):
\begin{matlab}
%}
xtGy = @(x,G,y) shiftdim(sum(sum( ...
       permute(x,[1 4 2 3]).*G.*shiftdim(y,-1) ...
       )),2);
%{
\end{matlab}
The moving mesh \ode{}s~\eqref{mm2DiscdXdt} are then coded
\cite[(26)]{Huang98} as
(should \verb|sqrt(g1)| be \verb|sqrt(g1tilde)|??)
\begin{matlab}
%}
brace1 = xtGy(Xvr_J(:,I,J),G1r(:,:,I,J),Xvr_J(:,I,J))...
              ./Jacr(I,J)./g1r(I,J) ...
        -xtGy(Xvr_J(:,Im,J),G1r(:,:,Im,J),Xvr_J(:,Im,J))...
              ./Jacr(Im,J)./g1r(Im,J) ...
        -xtGy(Xvr_I(:,I,J),G1u(:,:,I,J),Xvr_J(:,I,J))...
              ./Jacu(I,J)./g1r(I,J) ...
        +xtGy(Xvr_I(:,I,Jm),G1u(:,:,I,Jm),Xvr_J(:,I,Jm))...
              ./Jacu(I,Jm)./g1r(I,Jm);
brace2 =-xtGy(Xvr_J(:,I,J),G2r(:,:,I,J),Xvr_I(:,I,J))...
              ./Jacr(I,J)./g2r(I,J) ...
        +xtGy(Xvr_J(:,Im,J),G2r(:,:,Im,J),Xvr_I(:,Im,J))...
              ./Jacr(Im,J)./g2r(Im,J) ...
        +xtGy(Xvr_I(:,I,J),G2u(:,:,I,J),Xvr_I(:,I,J))...
              ./Jacu(I,J)./g2u(I,J) ...
        -xtGy(Xvr_I(:,I,Jm),G2u(:,:,I,Jm),Xvr_I(:,I,Jm))...
              ./Jacu(I,Jm)./g2u(I,Jm);
testy(brace1,ind,'brace1')
testy(brace2,ind,'brace2')
M.Vx = shiftdim( ...
     -squeeze(X_I)./(sqrt(g1).*Jac).*brace1 ...
     -squeeze(X_J)./(sqrt(g2).*Jac).*brace2 ...
     ,-4)/patches.mmTime;
M.Vy = shiftdim( ...
     -squeeze(Y_I)./(sqrt(g1).*Jac).*brace1 ...
     -squeeze(Y_J)./(sqrt(g2).*Jac).*brace2 ...
     ,-4)/patches.mmTime;
testy(M.Vx ,4+ind,'M.Vx')
testy(M.Vy ,4+ind,'M.Vy')
%{
\end{matlab}
\end{subequations}

\begin{matlab}
%}
end% switch theMethod
%{
\end{matlab}




\subsection{Evaluate system differential equation}
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
 