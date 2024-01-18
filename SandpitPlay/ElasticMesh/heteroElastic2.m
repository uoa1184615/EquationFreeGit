% The aim is to simulate a heterogeneous forced 2D elastic sheet
% using an Equation-free Patch Scheme, and fixed ends.  Here
% we code patches to use a microscale staggered grid of the
% microscale heterogeneous elasticity PDEs.  This function
% computes the time derivatives given that interpolation has
% provided the patch-edge values. AJR, 20 Aug 2022 -- 25 Oct
% 2023
% Adapted from elastic2Dstaggered.m AJR, 20 Aug 2022 -- 18 Mar 2023
% JEB, Dec 2023
%{
\section{\texttt{elastic2Dstaggered()}: \(\D t{}\) of 2D
heterogeneous elastic patch on staggered grid}
\label{secelastic2Dstaggered}

Let's try a staggered microgrid for patches of heterogeneous
2D elasticity forming a 2D beam. \cref{figpatchgridv} draw
the microgrid in a patches: the microgrid is akin to that by
\cite{Virieux86}.
\begin{matlab}
%}
function [Ut]=heteroElastic2(t,U,patches)
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
[nx,~,~,~,Nx] = size(patches.x);
[~,ny,~,~,~,Ny] = size(patches.y);
U = reshape(U,nx,ny,4,1,Nx,Ny);
%{
\end{matlab}
\input{Figs/figpatchgridv} 
Set microgrid parameters of patches like
\cref{figpatchgridv}.
\begin{matlab}
%}
dy=diff(patches.y(2:3)); 
dx=diff(patches.x(2:3));
i = 2:nx-1;  % indices for interior x
j = 2:ny-1;  % indices for interior y
ix=1:nx-1; % xtra index for x
jx=1:ny-1; % xtra index for y
%{
\end{matlab}

\paragraph{Optional boundary conditions}
Apply fixed-end BCs (Dirichlet) of all four fields at
extreme left/right and top/bottom edges of extreme 
left/right and top/bottom patches
\begin{matlab}
%}
if ~patches.periodic(1)
    U( 1,:,:,:, 1,:) = 0;
    U(nx,:,:,:,Nx,:) = 0;
end%if
if ~patches.periodic(2)
    U(:, 1,:,:,:, 1) = 0;
    U(:,ny,:,:,:,Ny) = 0;
end%if
%{
\end{matlab}

Then split input field into the four physical fields.
\begin{matlab}
%}
u = U(:,:,1,:,:,:);
v = U(:,:,2,:,:,:);
ut= U(:,:,3,:,:,:);
vt= U(:,:,4,:,:,:);
%{
\end{matlab}
Unpack the physical, microscale heterogeneous, elastic
parameters.
\begin{matlab}
%}
mu  =patches.cs(:,:,1:2);
lamb=patches.cs(:,:,3:4);
visc=patches.viscosity; 
if t<0, visc=visc, end
%{
\end{matlab}

Temporary trace print.
\begin{matlab}
%}
   if t<0
   usz=size(u), vsz=size(v), utsz=size(ut)
   vtsz=size(vt), musz=size(mu), lambsz=size(lamb)
   if t<-1
   us=squeeze(u), vs=squeeze(v), uts=squeeze(ut)
   vts=squeeze(vt), mus=squeeze(mu), lambs=squeeze(lamb)
   end
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

\begin{itemize}
\item Compute stresses \(\sigma_{xy}\) at \se-points.
\begin{matlab}
%}
sxy = mu(ix,jx,2).*( diff(u(ix,:,:,:,:,:),1,2)/dy ...
                    +diff(v(:,jx,:,:,:,:))/dx ); %edit from: diff(v(:,jx+1,:,:,:,:))/dx )
%{
\end{matlab}
\item Compute stresses \(\sigma_{yy}\) at \nw-points
(omitting leftmost~\(\sigma_{yy}\)). Little tricky with
indices of \nw-elasticity parameters as they cross over the
`grid lines'.  Apply zero normal stress between the top-edge
and the next-to-top-edge positions along the beam.  Here we
set~\(v\) at the top-edge so that \(\sigma_{yy}\) computes
to zero at the next-to-top-edge.
\begin{matlab}
%}
if ~patches.periodic(2) % make syy zero at j=ny-1 => sxx is OK there
v(ix+1,ny,:,:,:,Ny) = v(ix+1,ny-1,:,:,:,Ny) ...
    -lamb(:,ny-1,1).*diff(u(:,ny-1,:,:,:,Ny))/dx ...
    *dy./(2*mu(:,ny-1,1)+lamb(:,ny-1,1));
end
syy = (2*mu(:,jx,1)+lamb(:,jx,1)).*diff(v(ix+1,:,:,:,:,:),1,2)/dy ...
      +lamb(:,jx,1).*diff(u(:,jx+1,:,:,:,:))/dx; % edited from: diff(u(:,jx,:,:,:,:))/dx;
   
%{
\end{matlab}
\item Compute~\(\sigma_{xx}\) at \nw-points (omitting
leftmost \(i=1\)).
\begin{matlab}
%}
sxx = lamb(:,jx,1).*diff(v(ix+1,:,:,:,:,:),1,2)/dy ...
     +(2*mu(:,jx,1)+lamb(:,jx,1)).*diff(u(:,jx+1,:,:,:,:))/dx; % edited from: diff(u(:,jx,:,:,:,:))/dx;
%{
\end{matlab}
\end{itemize}


\paragraph{Top-bottom boundary conditions} At the top-bottom
of the beam we need \(\sigma_{xy}=\sigma_{yy}=0\).  Here set 
zero at microscale zig-zag locations along top and bottom.
\begin{matlab}
%}
if ~patches.periodic(2)
sxy(:, 1,:,:,:, 1) = 0;
syy(:, 1,:,:,:, 1) = 0;
sxy(:,ny-1,:,:,:,Ny) = 0;
%syy(:,ny-1,:,:,:,Ny) = 0; % enforced above
end
%{
\end{matlab}
Temporary trace print.
\begin{matlab}
%}
   if t<0
   sxysz=size(sxy), syysz=size(syy), sxxsz=size(sxx)
   if t<-1
   sxys=squeeze(sxy), syys=squeeze(syy), sxxs=squeeze(sxx)
   end
   tmp=sxx(:,1,:,:,:,1); sxx(:,1,:,:,:,1)=0;
   assert(all(~isnan(sxx(:))),'sxx has some nans')
   assert(all(~isnan(sxy(:))),'sxy has some nans')
   assert(all(~isnan(syy(:))),'syy has some nans')
   sxx(:,1,:,:,:,1)=tmp;
   end   
%{
\end{matlab}

\paragraph{Second derivatives for viscosity} This code is
for constant viscosity---if spatially varying, then modify
each to be two first-derivatives in each direction. The
\(x\)-derivatives are well-defined second differences.
\begin{matlab}
%}
utxx = visc*diff(ut(:,j,:,:,:,:),2,1)/dx^2; 
vtxx = visc*diff(vt(:,j,:,:,:,:),2,1)/dx^2; 
%{
\end{matlab}
On \cref{figpatchgridv}, and for the purpose of viscosity,
assume \(\D y{\dot u}=\D y{\dot v}=0\) on the top-bottom (we
just need `viscosity' to be some phenomenological
dissipation).
\begin{matlab}
%}
if ~patches.periodic(2)
ut(i, 1,:,:,:, 1) = ut(i,  2 ,:,:,:, 1);
ut(i,ny,:,:,:,Ny) = ut(i,ny-1,:,:,:,Ny);
vt(i, 1,:,:,:, 1) = vt(i,  2 ,:,:,:, 1);
vt(i,ny,:,:,:,Ny) = vt(i,ny-1,:,:,:,Ny);
end
utyy = visc*diff(ut(i,:,:,:,:,:),2,2)/dy^2;
vtyy = visc*diff(vt(i,:,:,:,:,:),2,2)/dy^2;
%{
\end{matlab}
Temporary trace print.
\begin{matlab}
%}
   if t<0
   utxxsz=size(utxx), vtxxsz=size(vtxx)
   utyysz=size(utyy), vtyysz=size(vtyy)
   if t<-1
   utxxs=squeeze(utxx), vtxxs=squeeze(vtxx)
   utyys=squeeze(utyy), vtyys=squeeze(vtyy)
   end
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
Ut=zeros(nx,ny,4,1,Nx,Ny);
Ut(i,j,1,:) = ut(i,j,:,:); % dudt = ut
Ut(i,j,2,:) = vt(i,j,:,:); % dvdt = vt     
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
if t<0
   utt = diff(sxx(:,j,:,:,:,:))/dx +diff(sxy(i,:,:,:,:,:),1,2)/dy ...
       +(utxx+utyy);
   vtt = diff(syy(i-1,:,:,:,:,:),1,2)/dy +diff(sxy(:,j-1,:,:,:,:))/dx ...
       +(vtxx+vtyy) + patches.f(i,j,:,:,:,:);
   uttsz=size(utt), vttsz=size(vtt)
   if t<-1, utts=squeeze(utt), vtts=squeeze(vtt), end
   assert(all(~isnan(utt(:))),'utt has some nans')
   assert(all(~isnan(vtt(:))),'vtt has some nans')
end   
Ut(i,j,3,:,:,:) ...
    = diff(sxx(:,j-1,:,:,:,:))/dx +diff(sxy(i,:,:,:,:,:),1,2)/dy ...
    +(utxx+utyy); % edited from diff(sxx(:,j,:,:,:,:))/dx +diff(sxy(i,:,:,:,:,:),1,2)/dy
Ut(i,j,4,:,:,:) ...
    = diff(syy(i-1,:,:,:,:,:),1,2)/dy +diff(sxy(:,j,:,:,:,:))/dx ...
    +(vtxx+vtyy) + patches.f(i,j,:,:,:,:); % edited from +diff(sxy(:,j-1,:,:,:,:))/dx

%{
\end{matlab}
End of function.
%}