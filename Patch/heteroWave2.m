% Computes the time derivatives of forced heterogeneous
% waves (slightly damped) in 2D on patches.  AJR, Aug 2021 
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroWave2()}: heterogeneous Waves}
\label{sec:heteroWave2}

This function codes the lattice heterogeneous waves inside
the patches.  The forced wave \pde\ is
\begin{equation*}
u_t=v,\quad v_t=\grad(a\divv u)+f
\end{equation*}
for scalars~\(a(t,x,y)\) and~\(f(t,x,y)\) where~\(a\) has
microscale variations.  For 6D input arrays~\verb|u|,
\verb|x|, and \verb|y| (via edge-value interpolation of
\verb|patchSys2|, \cref{sec:patchSys2}), computes the
time derivative at each point in the interior of a patch,
output in~\verb|ut|.  The four 2D arrays of heterogeneous
interaction coefficients,~$c_{ijk}$, have previously been
stored in~\verb|patches.cs| (3D). 

Supply patch information as a third argument (required by
parallel computation), or otherwise by a global variable.
\begin{matlab}
%}
function ut = heteroWave2(t,u,patches)
  if nargin<3, global patches, end
%{
\end{matlab}
Microscale space-steps, and interior point indices.  
\begin{matlab}
%}
  dx = diff(patches.x(2:3));  % x micro-scale step
  dy = diff(patches.y(2:3));  % y micro-scale step
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
  assert(max(abs(u(:)))<9999,"u-field exploding")
%{
\end{matlab}

Form coefficients here---odd periodic extension.  To avoid
slight errors in periodicity (in full domain simulation),
first adjust any coordinates crossing \(x=\pm1\) or \(y=\pm
1\).
\begin{matlab}
%}
x=patches.x; y=patches.y;
l=find(abs(x)>1); x(l)=x(l)-sign(x(l))*2;
l=find(abs(y)>1); y(l)=y(l)-sign(y(l))*2;
%{
\end{matlab}
Then set at this time three possible forcing functions,
although only use one depending upon \verb|patches.eff|.
Forcing~\(f_1\) and~\(f_2\) are as specified by \S5.1 of
\cite{Maier2021}, whereas~\(f_3\) here is~\(f\) in
their \S5.2.
\begin{matlab}
%}
f1 = ( (abs(x)>0.4)*(20*t+230*t^2) ...
      +(abs(x)<0.4)*(100*t+2300*t^2) ).*sign(x).*sign(y);
f2 = 20*t*x.*(1-abs(x)).*y.*(1-abs(y)) ...
    +230*t^2*(sign(y).*x.*(1-abs(x))+sign(x).*y.*(1-abs(y)));
f3 = (5*t+50*t^2)*sin(pi*x).*sin(pi*y);
%{
\end{matlab}
Also set the heterogeneous interactions at this time.
\begin{matlab}
%}
ax = (patches.cs(:,:,1)+sin(2*pi*t)) ...
   .*(patches.cs(:,:,2)+sin(2*pi*t));
ay = (patches.cs(:,:,3)+sin(2*pi*t)) ...
   .*(patches.cs(:,:,4)+sin(2*pi*t));
%{
\end{matlab}

Reserve storage (using \verb|nan+u| appears quickest), and
then assign time derivatives for interior patch values due
to the heterogeneous interaction and forcing.
\begin{matlab}
%}
  ut = nan+u;  % preallocate output array
  ut(i,j,1,:) = u(i,j,2,:);
  ut(i,j,2,:) ...
  = diff(ax(:,j).*diff(u(:,j,1,:),1),1)/dx^2 ...
   +diff(ay(i,:).*diff(u(i,:,1,:),1,2),1,2)/dy^2 ...
   +(patches.eff==1)*f1(i,j,:,:) ...
   +(patches.eff==2)*f2(i,j,:,:) ...
   +(patches.eff==3)*f3(i,j,:,:) ...
   + 1e-4*(diff(u(:,j,2,:),2,1)/dx^2+diff(u(i,:,2,:),2,2)/dy^2); 
end% function
%{
\end{matlab}
In the last line above, the slight damping of~\(10^{-4}\)
causes microscale modes to decay at rate~\(e^{-28t}\), with
frequencies~\(2000\)--\(5000\), whereas macroscale modes
decay with rates roughly~\(0.0005\)--\(0.05\) with
frequencies~\(10\)--\(100\).  This slight damping term may
correspond to the weak damping of the backward Euler scheme
adopted by \cite{Maier2021} for time integration.
%}
