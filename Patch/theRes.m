% This functions converts a vector of values into the
% interior values of the patches, then evaluates the time
% derivative of the system at $t=1$, and returns the vector
% of patch-interior time derivatives.  Now for AutoDiff.
% AJR, 1 Feb 2023 -- 14 Jul 2026
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{theRes()}: wrapper function to zero for equilibria}
\label{sec:theRes}
This functions converts a vector of values into the interior
values of the patches, then evaluates the time derivative of
the system at time \(t=1\), and returns the vector of
patch-interior time derivatives.
\begin{matlab}
%}
function f=theRes(u)
  global patches
  switch numel(size(patches.x))
    case 4, pSys = @patchSys1;
            szv = size(patches.x);
            szii = 2:3;
    case 5, pSys = @patchSys2;
            szv = size(patches.x+patches.y);
            szii = 3:4;
    case 6, pSys = @patchSys3;
            szv = size(patches.x+patches.y+patches.z);
            szii = 4:5;
    otherwise error('number of dimensions is somehow wrong')
  end%switch
  l = length(patches.nEdge);
  szi = szv; szi(1:l) = szi(1:l)-2*patches.nEdge; % omit edges from count
  nVars = prod(size(u))/prod(szi)/patches.nEnsem;
  szv(szii) = [nVars patches.nEnsem]; % count nVars and ensem
  assert(all(szv==round(szv)),['bad non-integer array size ' num2str(szv)])
  v = nan(szv,'like',u);
  v(patches.i) = u;
  f = pSys(1,v(:),patches);
  f = f(patches.i);
end%function theRes
%{
\end{matlab}
%}
