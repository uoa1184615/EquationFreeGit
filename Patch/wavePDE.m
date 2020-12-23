% Microscale discretisation of the 2D ideal wave PDE inside
% 2D patches in space.  Used by the example wave2D.m
% AJR, 4 Apr 2019 -- Nov 2020
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\subsection{\texttt{wavePDE()}: Example of simple wave PDE inside patches}
As a microscale discretisation of \(u_{tt}=\delsq(u)\), so
code \(\dot u_{ijkl}=v_{ijkl}\) and \(\dot v_{ijkl}
=\frac1{\delta x^2} (u_{i+1,j,k,l} -2u_{i,j,k,l}
+u_{i-1,j,k,l}) + \frac1{\delta y^2} (u_{i,j+1,k,l}
-2u_{i,j,k,l} +u_{i,j-1,k,l})\).
\begin{matlab}
%}
function uvt = wavePDE(t,uv,patches)
  dx = diff(patches.x(1:2));   
  dy = diff(patches.y(1:2));  % microscale spacing
  i = 2:size(uv,1)-1;  
  j = 2:size(uv,2)-1; % interior patch-points
  uvt = nan+uv;  % preallocate storage
  uvt(i,j,1,:) = uv(i,j,2,:);
  uvt(i,j,2,:) = diff(uv(:,j,1,:),2,1)/dx^2 ...
                +diff(uv(i,:,1,:),2,2)/dy^2;
end
%{
\end{matlab}
%}
