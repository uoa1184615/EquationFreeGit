% Solve for steady state of multiscale heterogeneous diffusion
% in 2D on patches as an example application, varied from
% example of section 5.1 of Allaire & Brizzi (2005).  
% **** just exploratory AJR, 3 Feb 2026
%{
\section{\texttt{allaireDiffEquil2}: equilibrium of a 2D
multiscale heterogeneous diffusion via small patches}
\label{sec:allaireDiffEquil2}

Here we find the steady state~\(u(x,y)\) to the
heterogeneous \pde\ \cite[inspired by][\S5.1]{Allaire2005}
\begin{equation*}
u_t=\divv[a(x,y)\grad u]-1,
\end{equation*}
on square domain \([0,1]^2\) with zero-Dirichlet BCs, for
coefficient `diffusion' matrix, varying with
period~\(\epsilon\) of (their~p.13)
\begin{equation*}
a:=\frac1{[2+P\sin(2\pi x/\epsilon)]\cdot[2+P\sin(2\piy/\epsilon)]},
\quad P=1.8\,,
\end{equation*}
\cref{fig:allaireDiffEquil2} shows solutions have some 
nice microscale wiggles reflecting the heterogeneity.
\begin{figure}
\centering\caption{\label{fig:allaireDiffEquil2}%
Equilibrium of the macroscale diffusion problem of Allaire
with boundary conditions of Dirichlet zero-value except for
\(x=0\) which is Neumann (\cref{sec:allaireDiffEquil2}). 
Here the patches have a Chebyshev-like spatial distribution.
The patch size is chosen large enough to see within.}
\includegraphics[scale=0.8]{Figs/allaireDiffEquil2}
\end{figure}

Clear, and initiate globals. 
\begin{matlab}
%}
clear all
global patches
%global OurCf2eps, OurCf2eps=true %option to save plot
%{
\end{matlab}


First establish the microscale heterogeneity has
micro-period~\verb|mPeriod| on the spatial micro-grid
lattice.  Then \verb|configPatches2| replicates the
heterogeneity to fill each patch.  (These diffusion
coefficients should really recognise the half-grid-point
shifts, but let's not bother.)
\begin{matlab}
%}
mPeriod = 6
x = (0.5:mPeriod)'/mPeriod;  y=x';
a = 1./(2+1.8*sin(2*pi*x))./(2+1.8*sin(2*pi*y));
diffusivityRange = [min(a(:)) max(a(:))]
%{
\end{matlab}
Set the periodicity~\(\epsilon\), here to match Fig.1 of 
Allaire2005.
\begin{matlab}
%}
epsilon = 0.01
dx = epsilon/mPeriod
nPeriodsPatch = 1 % any integer
nSubP = nPeriodsPatch*mPeriod+2 % when edgy int
%{
\end{matlab}



\paragraph{Patch configuration} 
Loop over increasing number of patches to assess errors
\begin{matlab}
%}
nPatchs=[]; rmsErrs=[];
nPatch = 3; % to start with five 
for iPatch=1:2
    nPatch=2*nPatch-1
    nPatchs(iPatch)=nPatch; % store
%{
\end{matlab}

Choose Dirichlet boundaries in coordination
with micro-code in \cref{sec:allaireDiffForce2}
\begin{matlab}
%}
Dom.bcOffset = zeros(2);
%{
\end{matlab}
Say use \(7\times7\) patches in \((0,1)^2\), fourth order
interpolation, and either `equispace' or `chebyshev':
\begin{matlab}
%}
Dom.type='equispace';
configPatches2(@allaireDiffForce2,[0 1],Dom ...
    ,nPatch ,4 ,dx ,nSubP ,'EdgyInt',true ,'hetCoeffs',a );
%{
\end{matlab}




\paragraph{Solve for steady state} Set initial guess of
zero, with \verb|NaN| to indicate patch-edge values.
Index~\verb|i| are the indices of patch-interior points, and
the number of unknowns is then its length.
\begin{matlab}
%}
u0 = zeros(nSubP,nSubP,1,1,nPatch,nPatch);
if iPatch>1%interpolate interpolation of previous solution
    u0(:,:,:,:,1:2:end,1:2:end) = up;
    u0(:,:,:,:,2:2:end,1:2:end) = (u0(:,:,:,:,1:2:end-2,1:2:end) ...
                                  +u0(:,:,:,:,3:2:end  ,1:2:end))/2;
    u0(:,:,:,:,:,2:2:end) = (u0(:,:,:,:,:,1:2:end-2) ...
                            +u0(:,:,:,:,:,3:2:end  ))/2;
end%if iPatch
u0([1 end],:,:) = nan;  u0(:,[1 end],:) = nan;
patches.i = find(~isnan(u0));
nVariables = numel(patches.i)
%{
\end{matlab}
Solve by iteration.  Use \verb|fsolve| for simplicity and
robustness (and using \verb|optimoptions| to omit trace
information), via the generic patch system wrapper
\verb|theRes| (\cref{sec:theRes}), and give magnitudes.
Presumably the interpolation generates large residuals 
around the patch edges, but the residual would be small 
within each patch.  Whereas if just an initial-zero field 
is used then there is no such interpolation residual, only 
the uniform interior residual.
\begin{matlab}
%}
meanInitialRes=mean(abs( theRes(u0(patches.i)) ))
tic;
if 1
    uSoln = fsolve(@theRes,u0(patches.i) ...
            ,optimoptions('fsolve','Display','off')); 
else
    restart=floor(sqrt(nVariables)); maxit=restart;
    b=-theRes(zeros(nVariables,1));
    uSoln = gmres(@(u) theRes(u)+b ,b ...
            ,restart,[],maxit,[],[],u0(patches.i)); 
end
solnTime = toc
meanResidual = mean(abs( theRes(uSoln) ))
normSoln = norm(uSoln)
rmsSoln = rms(uSoln)
%{
\end{matlab}
Store the solution vector into the patches, and interpolate,
but have not bothered to set boundary values so they stay
NaN from the interpolation.
\begin{matlab}
%}
if iPatch<99, u0(patches.i) = uSoln; end
u0 = patchEdgeInt2(u0);
%{
\end{matlab}
Error of previous one compared to this one.
\begin{matlab}
%}
  if iPatch>1,
    errp = up-u0(:,:,:,:,1:2:end,1:2:end);
    rmsErrs(iPatch-1) = rms(errp(:),"omitnan")/rmsSoln
  end%if iPatch
  up = u0;%save for next iteration
%{
\end{matlab}





\paragraph{Draw solution profile} Separate patches with
NaNs, then reshape arrays to suit 2D space surface plots.
\begin{matlab}
%}
figure(iPatch), clf, colormap(0.8*hsv)
patches.x(end+1,:,:)=nan;  u0(end+1,:,:)=nan;  
patches.y(:,end+1,:)=nan;  u0(:,end+1,:)=nan;
u = reshape(permute(squeeze(u0),[1 3 2 4]) ...
    , [numel(patches.x) numel(patches.y)]);
%{
\end{matlab}
Draw the patch solution surface, with boundary-values
omitted as already~\verb|NaN| by not bothering to set them.
\begin{matlab}
%}
mesh(patches.x(:),patches.y(:),u'); view(60,55) 
xlabel('space $x$'), ylabel('space $y$'), zlabel('$u(x,y)$')
drawnow
ifOurCf2eps(mfilename) %optionally save plot
end%for iPatch
%{
\end{matlab}






\subsection{\texttt{allaireDiffForce2()}: microscale
discretisation inside patches of forced diffusion PDE}
\label{sec:allaireDiffForce2}

This function codes the lattice heterogeneous diffusion of
the \pde\ inside the patches.  For 6D input arrays~\verb|u|,
\verb|x|, and~\verb|y|, computes the time derivative at each
point in the interior of a patch, output in~\verb|ut|.   
\begin{matlab}
%}
function ut = allaireDiffForce2(t,u,patches)
  dx = diff(patches.x(2:3));  % x space step
  dy = diff(patches.y(2:3));  % y space step
  i = 2:size(u,1)-1; % x interior points in a patch
  j = 2:size(u,2)-1; % y interior points in a patch
  ut = nan+u;         % preallocate output array
%{
\end{matlab}
Set Dirichlet boundary value of zero around the square
domain.
\begin{matlab}
%}
  u( 1 ,:,:,:, 1 ,:)=0; % left edge of left patches
  u(end,:,:,:,end,:)=0; % right edge of right patches
  u(:, 1 ,:,:,:, 1 )=0; % bottom edge of bottom patches
  u(:,end,:,:,:,end)=0; % top edge of top patches
%{
\end{matlab}
Compute the time derivatives via stored forcing and
coefficients.  Easier to code by conflating the last four
dimensions into the one~\verb|,:|.
\begin{matlab}
%}
  ut(i,j,:) = diff(patches.cs(:,j).*diff(u(:,j,:)))/dx^2 ...
      + diff(patches.cs(i,:).*diff(u(i,:,:),1,2),1,2)/dy^2 ...
      - 1; 
end%function allaireDiffForce2
%{
\end{matlab}
%}
