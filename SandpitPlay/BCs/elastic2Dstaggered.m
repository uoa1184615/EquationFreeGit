% The aim is to simulate a heterogeneous 2D beam using an
% Equation-free Patch Scheme.  Here we code patches to use a
% microscale staggered grid of the microscale heterogeneous
% elasticity PDEs.  This function computes the time
% derivatives given that interpolation has provided the
% patch-edge values, except we code boundary conditions on
% extreme patches. AJR, 20 Aug 2022 -- 4 Feb 2023
%!TEX root = doc.tex
%{
\subsection{\texttt{elastic2Dstaggered()}: \(\D t{}\) of 2D
heterogeneous elastic patch on staggered grid}
\label{secelastic2Dstaggered}

Let's try a staggered microgrid for patches of heterogeneous
2D elasticity forming a 2D beam. \cref{figpatchgridv} draw
the microgrid in a patches: the microgrid is akin to that by
\cite{Virieux86}.
\begin{matlab}
%}
function [Ut]=elastic2Dstaggered(t,U,patches)
%{
\end{matlab}
\paragraph{Input} \begin{itemize}
\item \verb|t| is time (real scalar).
\item \verb|U| is vector of \((u,v,\dot u,\dot v)\) values
in each and every patch.
\item \verb|patches| data structure from \verb|configPatches|,
with extra physical parameters.
\end{itemize}

\paragraph{Output} \begin{itemize}
\item \verb|Ut| is corresponding vector of time derivatives
of~\verb|U|
\end{itemize}


\paragraph{Unpack the data vector} Form data vector into a
4D array: \verb|nx|~is the number of points along a patch;
\verb|Nx|~is the number of patches.
\begin{matlab}
%}
[nx,~,~,Nx] = size(patches.x);
nEnsem = patches.nEnsem;
nVars = round(numel(U)/numel(patches.x)/nEnsem);
U = reshape(U,nx,nVars,nEnsem,Nx);
%{
\end{matlab}
\input{Figs/figpatchgridv} 
Set microgrid parameters of patches like
\cref{figpatchgridv}.
\begin{matlab}
%}
ny=(nVars+2)/4;
nny=2*ny-1;
dx=diff(patches.x(2:3));
dy=dx; % now set by new sim script
i =2:nx-1;  % indices for interior x
ix=1:nx-1;
ju=1:ny-1; % indices y structures of u, u_t
jv=1:ny;   % indices y structures of v, v_t
jse=1:2:nny; jnw=2:2:nny; % indices of elasticity y structure
%{
\end{matlab}
Then split input field into the four physical fields.
\begin{matlab}
%}
u = U(:,ju       ,:,:);
v = U(:,jv+ju(end),:,:);
ut= U(:,ju+nny,:,:);
vt= U(:,jv+nny+ju(end),:,:);
%{
\end{matlab}
Unpack the physical, microscale heterogeneous, elastic
parameters.
\begin{matlab}
%}
mu  =patches.cs(:,1:nny);
lamb=patches.cs(:,2*ny:end);
visc=patches.viscosity; 
%{
\end{matlab}

Temporary trace print.
\begin{matlab}
%}
   if t<0
   usz=size(u),us=squeeze(u);
   vsz=size(v),vs=squeeze(v);
   utsz=size(ut),uts=squeeze(ut);
   vtsz=size(vt),vts=squeeze(vt);
   musz=size(mu),mus=squeeze(mu);
   lambsz=size(lamb),lambs=squeeze(lamb);
   end   
%{
\end{matlab}


\paragraph{Mathematical expression of elasticity} Strain
tensor \(\varepsilon_{ij} := \tfrac12 (u_{i,j} +u_{j,i})\)
is symmetric. Remember stress~\(\sigma_{ij}\) is the force
in the \(j\)th~direction across a surface whose normal is
the \(i\)th~direction. So across a surface whose normal
is~\nv, the stress is~\(\nv^T\sigma\),
equivalently~\(n_i\sigma_{ij}\).


Denote the 2D displacement vector by \(\uv=(u,v)\).
\footnote{Adapted from Wikipedia: linear elasticity:
Derivation of Navier--Cauchy equations} First, consider the
\(x\)-direction.  Substitute the strain-displacement
equations into the equilibrium equation in the
\(x\)-direction to get
\begin{align*}&
\sigma_{xx} =2\mu \varepsilon_{xx} +\lambda
(\varepsilon_{xx} +\varepsilon_{yy}) =(2\mu+\lambda) {\D xu}
+\lambda {\D yv}
\\&
\sigma_{xy} = 2\mu \varepsilon_{xy} =\mu \left({\D yu} +{\D
xv}\right)
\end{align*}

\paragraph{Discretise in space} Form equi-spaced microscale
grid of spacing~\(dx,dy\) and indices~\(i,j\) in the
\(x,y\)-directions respectively (\cref{figpatchgridv}). 
Field values are at various `quarter points' within the
microgrid elements (\ne, \nw, \sw, \se).
\begin{itemize}
\item locate \(u\) at \ne-points
\item locate \(v\) at \sw-points
\item evaluate diagonal strains \(\varepsilon_{xx},
\varepsilon_{yy}\) at \nw-points,  also~\(\lambda,\mu\)
\item evaluate off-diagonal strains \(\varepsilon_{xy}
=\varepsilon_{yx}\) at \se-points, also~\(\mu\).
\item evaluate \pde{}s at the appropriate \sw,\ne~points, via
gradients of stresses from the \nw,\se~points.
\end{itemize}
All the above derivatives in space then involve simple
first differences from grid points half-spacing away. They do
\emph{not} involve any averaging. Such staggered spatial
discretisation is the best chance of best behaviour.

For interpolation between patches in the \(x\)-direction,
that edge values of~\(u,v\) are in zig-zag
positions makes absolutely no difference to the patch
interpolation schemes, since all coded schemes are invariant
to translations in the \(x\)-direction just so long as the
translation is the same across all patches for each of the
edge variables.

\paragraph{Physical boundary conditions on extreme patches}
Fixed boundaries are to simply code zero displacements and
zero time derivatives (velocities) at the boundaries.  First
at the left edge of the leftmost patch:
\begin{matlab}
%}
tweak=1;% choose simple zero, or accurate one??
u(1,:,:,1) =  +u(2,:,:,1)/5*tweak; 
v(1,:,:,1) =  -v(2,:,:,1)/3*tweak; 
ut(1,:,:,1)= +ut(2,:,:,1)/5*tweak; 
vt(1,:,:,1)= -vt(2,:,:,1)/3*tweak; 
%{
\end{matlab}
Second, at the right edge of the rightmost patch: first
alternative is fixed boundary.
\begin{matlab}
%}
if ~patches.stressFreeRightBC
u(nx,:,:,Nx) =  -u(nx-1,:,:,Nx)/3*tweak; 
v(nx,:,:,Nx) =  +v(nx-1,:,:,Nx)/5*tweak; 
ut(nx,:,:,Nx)= -ut(nx-1,:,:,Nx)/3*tweak; 
vt(nx,:,:,Nx)= +vt(nx-1,:,:,Nx)/5*tweak; 
%{
\end{matlab}
otherwise code a stress free boundary at the location of the
\(x\)-direction displacements and velocities (as for the
top-bottom boundaries).
\begin{matlab}
%}
else% patches.stressFreeRightBC
%{
\end{matlab}
Need \verb|sxy| (below) to be computed zero at the
bdry---assume \(u\) fattened by zero difference in~\(y\) at
top/bottom.  Later set \verb|sxx| to zero.
\begin{matlab}
%}
uDy = diff(u(nx-1,:,:,Nx),1,2);
v(nx,:,:,Nx) = v(nx-1,:,:,Nx)-dx/dy*[0 uDy 0];
%{
\end{matlab}
For phenomenological viscosity purposes assume \(\D x{\dot
u}=\D x{\dot v}=0\) at the free boundary.
\begin{matlab}
%}
vt(nx,:,:,Nx) = vt(nx-1,:,:,Nx);
ut(nx,:,:,Nx) = ut(nx-2,:,:,Nx);
end%if ~patches.stressFreeRightBC
%{
\end{matlab}



\paragraph{Patch-interior elasticity} For the case of
\cref{figpatchgridv}, fatten~\verb|u| using zero stress on
top-bottom, although not really needed it does mean
that~\(\sigma_{xy}\) computed next is zero on the
top-bottom.
\begin{matlab}
%}
u= [nan(nx,1,nEnsem,Nx)  u  nan(nx,1,nEnsem,Nx) ];
u(ix,1   ,:,:)=u(ix,2 ,:,:)+dy/dx*diff(v(:,1 ,:,:));
u(ix,ny+1,:,:)=u(ix,ny,:,:)-dy/dx*diff(v(:,ny,:,:));
%{
\end{matlab}
\begin{itemize}
\item Compute stresses \(\sigma_{xy}\) at \se-points.
\begin{matlab}
%}
sxy = mu(ix,jse).*( diff(u(ix,:,:,:),1,2)/dy+diff(v)/dx );
%{
\end{matlab}
\item Compute stresses \(\sigma_{yy}\) at \nw-points
(omitting leftmost~\(\sigma_{yy}\)), fattening the array to
cater for ghost nodes.  Little tricky with indices of
\nw-elasticity parameters as they cross over the `grid
lines'.
\begin{matlab}
%}
syy = nan(nx,ny+1,nEnsem,Nx);
syy(ix+1,ju+1,:,:) =  ...
      (2*mu(:,jnw)+lamb(:,jnw)).*diff(v(ix+1,:,:,:),1,2)/dy ...
      +lamb(:,jnw).*diff(u(:,2:ny,:,:))/dx;
%{
\end{matlab}
\item Compute~\(\sigma_{xx}\) at \nw-points (omitting
leftmost \(i=1\)).
\begin{matlab}
%}
sxx = lamb(:,jnw).*diff(v(ix+1,:,:,:),1,2)/dy ...
     +(2*mu(:,jnw)+lamb(:,jnw)).*diff(u(:,ju+1,:,:))/dx;
%{
\end{matlab}
\end{itemize}

\paragraph{Top-bottom boundary conditions} At the top-bottom
of the beam we need \(\sigma_{xy}=\sigma_{yy}=0\). Given the
micro-grid of \cref{figpatchgridv}, the first is already
catered for. Here use \(\sigma_{yy}=0\) at \sw~points on the
top-bottom to set ghost values. 
\begin{matlab}
%}
syy(:,   1,:,:) = -syy(:,2 ,:,:);
syy(:,ny+1,:,:) = -syy(:,ny,:,:);
%{
\end{matlab}
Similarly for any stress-free right boundary condition
(remembering \verb|sxx| currently omits leftmost \(i=1\)).
\begin{matlab}
%}
if patches.stressFreeRightBC
  sxx(nx-1,:,:,Nx) = -sxx(nx-2,:,:,Nx);
end
%{
\end{matlab}
Temporary trace print.
\begin{matlab}
%}
   if t<0
   sxysz=size(sxy),sxys=squeeze(sxy);
   syysz=size(syy),syys=squeeze(syy);
   sxxsz=size(sxx),sxxs=squeeze(sxx);
   end   
%{
\end{matlab}

\paragraph{Second derivatives for viscosity} This code is
for constant viscosity---if spatially varying, then modify
each to be two first-derivatives in each direction. The
\(x\)-derivatives are well-defined second differences.
\begin{matlab}
%}
utxx = visc*diff(ut,2,1)/dx^2; 
vtxx = visc*diff(vt,2,1)/dx^2; 
%{
\end{matlab}
On \cref{figpatchgridv}, and for the purpose of viscosity,
assume \(\D y{\dot u}=\D y{\dot v}=0\) on the top-bottom (we
just need `viscosity' to be some phenomenological
dissipation).
\begin{matlab}
%}
utyy = visc*diff(ut(i,[1 ju ny-1],:,:),2,2)/dy^2;
vtyy = visc*diff(vt(i,[2 jv ny-1],:,:),2,2)/dy^2;
%{
\end{matlab}
Temporary trace print.
\begin{matlab}
%}
   if t<0
   utxxsz=size(utxx),utxxs=squeeze(utxx);
   vtxxsz=size(vtxx),vtxxs=squeeze(vtxx);
   utyysz=size(utyy),utyys=squeeze(utyy);
   vtyysz=size(vtyy),vtyys=squeeze(vtyy);
   end   
%{
\end{matlab}


\paragraph{Time derivatives} The rate of change of
displacement is the provided velocities: \(\D tu=\dot u\)
and \(\D tv=\dot v\). Time derivative array should really be
initialised to~\verb|nan|, but \verb|ode15s| chokes on any
\verb|nan|s at the patch edges.
\begin{matlab}
%}
Ut=zeros(nx,nVars,nEnsem,Nx);
Ut(i,ju        ,:,:) = ut(i,:,:,:); % dudt = ut
Ut(i,jv+ju(end),:,:) = vt(i,:,:,:); % dvdt = vt     
%{
\end{matlab}
Substitute the stresses into the dynamical \pde{}s in the
\(x\)-direction to get
\begin{align*}&
\rho \D t{\dot u} 
= {\D {x}{\sigma_{xx}}} +{\D {y}{\sigma_{yx}}} 
+\kappa(\dot u_{xx}+\dot u_{yy})
+F_{x} \,,
\end{align*}
and correspondingly for the \(y\)-direction. Here, we
include some viscous dissipation of strength
\(\kappa=\verb|visc|\), but \emph{do not} justify it as
visco-elasticity.  First code above \(x\)-direction, then
the corresponding \(y\)-direction.
\begin{matlab}
%}
Ut(i,ju+nny,:,:) ...
    = diff(sxx)/dx +diff(sxy(i,:,:,:),1,2)/dy +(utxx+utyy); 
Ut(i,jv+nny+ju(end),:,:) ...
    = diff(syy(i,:,:,:),1,2)/dy +diff(sxy)/dx +(vtxx+vtyy);   
%{
\end{matlab}
End of function.
%}
