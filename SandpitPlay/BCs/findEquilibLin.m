% findEquilibLin find the patch equilibria of a linear
% system approximated by the patch scheme.  It forms the
% linear system matrix by numerical differentiation then
% solves the system explicitly.  This is computationally
% quicker, and allows us to explore the problem's multiscale
% nature, but takes more human time to code, using fsolve()
% is less human time.  AJR, 5 Jan 2023
%{
First establish the microscale heterogeneity has
micro-period \verb|mPeriod| on the lattice, and coefficients
to match the second example of Eckhardt2210.04536 \S6.2.1.
Set the phase of the heterogeneity so that each patch centre
is a point of symmetry of the diffusivity. Then the
heterogeneity is repeated to fill each patch. 
\begin{matlab}
%}
mPeriod = 6
y = linspace(0,1,mPeriod+1)';
a = 1./(2-cos(2*pi*y(1:mPeriod)))
%a=ones(mPeriod,1) % try constant diffusivity??
global microTimePeriod; microTimePeriod=0;
%{
\end{matlab}

Set the spatial period~\(\epsilon\), via
integer~\(1/\epsilon\), and other parameters.
\begin{matlab}
%}
nPatch = 11
nPeriodsPatch = 1 % any integer
rEpsilon = 200 % up to 2000 say
dx = 1/(mPeriod*rEpsilon+1)
nSubP = nPeriodsPatch*mPeriod+2
%{
\end{matlab}

Establish the global data struct~\verb|patches| for the
microscale heterogeneous lattice diffusion system
\cref{eq:hetroDiffF} solved on domain~\([0,1]\), with
\verb|nPatch| patches,  and say fourth order interpolation to
provide the edge-values.  Setting \verb|patches.EdgyInt|
true means the edge-values come from interpolating the
opposite next-to-edge values of the patches (not the
mid-patch values).  
\begin{matlab}
%}
global patches
ordCC = 2
configPatches1(@heteroDiffF,[0 1],'equispaced',nPatch ...
    ,ordCC,dx,nSubP,'EdgyInt',true,'hetCoeffs',a);
%{
\end{matlab}

Set the forcing coefficients, either the original parabolic,
or sinusoidal.
\begin{matlab}
%}
if 0 % given forcing
  patches.f1=2*( patches.x-patches.x.^2 );
  patches.f2=2*0.5+0*patches.x;
else% simple sine forcing 
  patches.f1=sin(pi*patches.x);
  patches.f2=pi/2*sin(pi*patches.x);
end%if
%{
\end{matlab}


\paragraph{Solve for equilibrium}
At time \(t=1\), the linear system is of the form \(\vec f(\vec u)=J\vec u+\vec f_0\)
Get the constant term in the system by evaluating at zero.
\begin{matlab}
%}
u0 = 0*patches.x;
f0 = patchSys1(1,u0);
%{
\end{matlab}
Put NaNs on the patch-edges, to then find all micro-grid points interior to the patches, and hence are the variables in~\(\vec u\).
\begin{matlab}
%}
u0([1 end],:,:,:)=nan;
i=find(~isnan(u0));
nJac=length(i)
%{
\end{matlab}
Create Jacobian~\(J\) column by column: since linear we numerically differentiate with unit vectors.  The plot of singular values shows: \(N-2\) `small' values, one for each interior patch; a pair of near-large values tentatively due to the two boundary patches; and the rest are `large' values attributed to sub-patch structures. 
\begin{matlab}
%}
deltau=1;
Jac=nan(nJac);
for j=1:nJac
  uj=u0;  uj(i(j))=deltau;
  dujdt=(patchSys1(1,uj)-f0)/deltau;
  Jac(:,j)=dujdt(i);
end
figure(1),semilogy(svd(Jac),'.'),ylabel('singular values')
Jacs=Jac-Jac'; Jacs(abs(Jacs)<1e-9)=0;
figure(2),spy(Jac,'o'),hold on,spy(Jacs,'rx'),hold off
assert(rcond(Jac)>1e-9,'Jacobian seems too ill-conditioned')
%{
\end{matlab}
Solve linear system.
\begin{matlab}
%}
u0(i) = -Jac\f0(i);
ue = squeeze(u0)
%{
\end{matlab}
Check the residual.
\begin{matlab}
%}
res = patchSys1(1,u0);
normRes = norm(res(i))
%{
\end{matlab}

Fin.
%}