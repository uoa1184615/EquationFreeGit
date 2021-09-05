% Computes the time derivatives of heterogeneous
% Landau--Lifshitz PDE on 1D lattice within spatial patches.
% From Leitenmaier & Runborg, arxiv.org/abs/2108.09463 and
% used by homoLanLif1D.m   AJR, Sep 2021
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{heteroLanLif1D()}: heterogeneous Landau--Lifshitz PDE}
\label{sec:heteroLanLif1D}

This function codes the lattice  heterogeneous
Landau--Lifshitz PDE \cite[(1.1)]{Leitenmaier2021} inside
patches in 1D space.  For 4D input array~\verb|M| storing
the three components of~\Mv\ (via edge-value interpolation
of \verb|patchSmooth1|, \cref{sec:patchSmooth1}), computes
the time derivative at each point in the interior of a
patch, output in~\verb|Mt|.  The column vector of
coefficients \(c_i=1+\tfrac12\sin(2\pi x_i/\epsilon)\) have
previously been stored in struct~\verb|patches.cs|.
\begin{itemize}
\item With \verb|ex5p1=0| computes the example \textsc{ex1}
\cite[p.6]{Leitenmaier2021}.
\item With \verb|ex5p1=1| computes the first 'locally
periodic' example \cite[p.27]{Leitenmaier2021}.
\end{itemize}


\begin{matlab}
%}
function Mt = heteroLanLif1D(t,M,patches)
  global alpha ex5p1
  dx = diff(patches.x(2:3));   % space step
  i = 2:size(M,1)-1;   % interior points in a patch
%{
\end{matlab}
Compute the heterogeneous \(\Hv:=\divv(a\grad\Mv)\)
\begin{matlab}
%}
  a = patches.cs ...
     +ex5p1*(0.1+0.25*sin(2*pi*(patches.x(2:end,:,:,:)-dx/2)+1.1));
  H = diff(a.*diff(M))/dx^2;
%{
\end{matlab}
At each microscale grid point, compute the cross-products
\(\Mv\times \Hv\) and \(\Mv\times(\Mv\times \Hv)\) to then 
give the time derivative \(\Mv_t=-\Mv\times \Hv -\alpha \Mv\times (\Mv\times \Hv)\) \cite[(1.1)]{Leitenmaier2021}:
\begin{matlab}
%}
  MH=nan+H; % preallocate for MxH
  MH(:,3,:,:) = M(i,1,:,:).*H(:,2,:,:)-M(i,2,:,:).*H(:,1,:,:);
  MH(:,2,:,:) = M(i,3,:,:).*H(:,1,:,:)-M(i,1,:,:).*H(:,3,:,:);
  MH(:,1,:,:) = M(i,2,:,:).*H(:,3,:,:)-M(i,3,:,:).*H(:,2,:,:);
  MMH=nan+H; % preallocate for MxMxH
  MMH(:,3,:,:)= M(i,1,:,:).*MH(:,2,:,:)-M(i,2,:,:).*MH(:,1,:,:);
  MMH(:,2,:,:)= M(i,3,:,:).*MH(:,1,:,:)-M(i,1,:,:).*MH(:,3,:,:);
  MMH(:,1,:,:)= M(i,2,:,:).*MH(:,3,:,:)-M(i,3,:,:).*MH(:,2,:,:);
  Mt = nan+M; % preallocate output array
  Mt(i,:,:,:) = -MH-alpha*MMH; 
end% function
%{
\end{matlab}
%}
