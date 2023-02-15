% This functions converts a vector of values into the
% interior values of the patches, then evaluates the time
% derivative of the system at $t=1$, and returns the vector
% of patch-interior time derivatives. AJR, 1 Feb 2023
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
            v=nan(size(patches.x));
    case 5, pSys = @patchSys2;
            v=nan(size(patches.x+patches.y));
    case 6, pSys = @patchSys3;
            v=nan(size(patches.x+patches.y+patches.z));
    otherwise error('number of dimensions is somehow wrong')
  end%switch
  v(patches.i) = u;
  f = pSys(1,v(:),patches);
  f = f(patches.i);
end%function theRes
%{
\end{matlab}
%}
