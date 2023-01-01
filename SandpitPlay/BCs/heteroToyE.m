% Computes the time derivatives of forced heterogeneous
% diffusion in 1D on patches.  AJR, Apr 2019 -- 20 Oct 2022
%!TEX root = doc.tex
%{
\subsection{\texttt{heteroToyE()}: forced heterogeneous toy elasticity}
\label{sec:heteroToyE}

This function codes the lattice heterogeneous toy elasticity
inside the patches.  Computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.
\begin{matlab}
%}
function uvt = heteroToyE(t,uv,patches)
  global b M vis i0 iN
%{
\end{matlab}
Separate state vector into displacement and velocity fields.
\begin{matlab}
%}
  u=uv(:,1:2,:,:); v=uv(:,3:4,:,:); % separate u and v=du/dt
%{
\end{matlab}
Compute the two different strain fields, and also a first
derivative for some optional viscosity.
\begin{matlab}
%}
  eps2 = diff(u)/(2*b);
  eps1 = [u(:,2,:,:)-u(:,1,:,:) u([2:end 1],1,:,:)-u(:,2,:,:)]/b;
  eps1(end,2,:,:)=nan; % as this value is fake
  vx1  = [v(:,2,:,:)-v(:,1,:,:) v([2:end 1],1,:,:)-v(:,2,:,:)]/b;
  vx1(end,2,:,:)=nan; % as this value is fake
%{
\end{matlab}
Set corresponding nonlinear stresses
\begin{matlab}
%}
  sig2 = eps2-M(2)*eps2.^3+eps2.^5;
  sig1 = eps1-M(1)*eps1.^3+eps1.^5;
%{
\end{matlab}
Preallocate output array, and fill in time derivatives of
displacement and velocity, from velocity and gradient of
stresses, respectively.
\begin{matlab}
%}
  uvt = nan+uv;          % preallocate output array
  i=2:size(uv,1)-1;
  % rate of change of position
  uvt(i,1:2,:,:) = v(i,:,:,:);  
  % rate of change of velocity +some artificial viscosity??
  uvt(i,3:4,:,:) = diff(sig2) ...
    +[ sig1(i,1,:,:)-sig1(i-1,2,:,:)  diff(sig1(i,:,:,:),1,2)] ... 
  +vis*[ vx1(i,1,:,:)-vx1(i-1,2,:,:)  diff(vx1(i,:,:,:),1,2) ]; 
%{
\end{matlab}
Maintain boundary value of \(u_i,\dot u_i\) by setting them
both to be constant in time, for both \(x_i=\pm b/2\).  If
\verb|i0|~is empty, then no boundary condition is set.
\begin{matlab}
%}
if ~isempty(i0), uvt(i0)=0; end
if ~isempty(iN), uvt(iN(3:4))=dLdt(t); end% vel=d/dt of end displacement
end% function
%{
\end{matlab}


\subsection{\texttt{dLdt()}: prescribed movement of length}
\begin{matlab}
%}
function Ld=dLdt(t)
Ld=-0.03*cos(t/20);
end
%{
\end{matlab}

%}
