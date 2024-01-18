% Use elastic2Dstaggered.m to compute the
% spectra of an unforced 2D beam which is variously
% macroscale periodic, or fixed-fixed.
% AJR, 21 Jun 2023 Rewritten Oct 2023
%{
Clear and set many parameters: spectral interpolation is a
start; spectra may be sensitive to viscosity.
\begin{matlab}
%}
clear all, close all
global patches 
Evar = 0.5
nPatch = [5 1]
nx = 7
ny = 6
iOrd = 0
edgyInt = true
viscosity = 1e-3
switch "Periodic" % choose beam end BCs
case "Periodic"
        xLim=[-pi pi]
        fn = @elastic2Dstaggered
        Dom.type='periodic '
        Dom.bcOffset = nan
case "FixFix"
        xLim=[0 pi]
        fn = @elastic2Dstaggered
        Dom.type='equispace'
        Dom.bcOffset = 0
end
eV = round(-log10(viscosity));
eE = round(-log(Evar)/log(2));
basename = mfilename +"N"+num2str(nPatch(1))+num2str(nPatch(2)) ...
        +"v"+num2str(eV)+"E"+num2str(eE)
%{
\end{matlab}
Configure the cross-beam structure.  Say use equispace
patches across the beam, but remember that self-adjoint only
with either one patch, or with full-cross-beam patches.  Set
\verb|bcOffset| so that none applied along the beam, but an
offset across the beam so that setting \(\delta_y\propto
1/(n_y-1)/Ny\) is full-cross-beam.   Get separated patches
across the beam with smaller~\(\delta_y\).
\begin{matlab}
%}
beamAspectRatio = 20
Dom.type(2,:) = 'equispace'
Dom.bcOffset(2) = 0.5
yLim = [-.5 .5]*diff(xLim)/beamAspectRatio 
dy = diff(yLim)/(ny-2)/nPatch(2)
dx = dy%*rand  % optionally randomly smaller than dy
%{
\end{matlab}


\paragraph{Set the microscale heterogeneity}
The microscale heterogeneity has micro-period that fills
each patch. Small heterogeneity of~\(0.01\) helps break some
symmetries, but essentially homogeneous.  Somewhat larger
at~\(0.05\) shows some minor heterogeneous effects.  Larger
at~\(0.2\) may make interpretation more difficult.
\begin{matlab}
%}
heteroPeriod = nx-2
E  = exp( Evar*(2*rand(nx-2,ny-2,2)-1) );
nu =  0.3 +0.1*(2*rand(nx-2,ny-2,2)-1);
if 0, E=1+0*E;  nu=0.3+0*nu; end% optionally homogenous
EQuartiles  = quantile(E(:) ,0:0.25:1)
nuQuartiles = quantile(nu(:),0:0.25:1)
%{
\end{matlab}
Compute corresponding \(\lambda,\mu\)-fields, and store so 
\verb|configPatches1| \(x\)-expands micro-heterogeneity to 
fill patch as necessary.
\begin{matlab}
%}
lambda=nu.*E./(1+nu)./(1-2*nu);
mu    =    E./(1+nu)./2;
cElas = cat(3,mu,lambda);
%{
\end{matlab}


\paragraph{Configure the patches}
\begin{matlab}
%}
disp('Using the following copy of configPatches2')
disp(which('configPatches2'))
configPatches2(fn,[xLim yLim],Dom,nPatch ...
    ,iOrd,[dx dy],[nx ny],'EdgyInt',edgyInt,'hetCoeffs',cElas);
patches.viscosity = viscosity; 
patches.yLim = yLim;
%{
\end{matlab}

Set the forcing to zero (for spectra).
\begin{matlab}
%}
patches.f = 0*patches.x+0*patches.y;
x = squeeze(patches.x)
y = squeeze(patches.y)
%{
\end{matlab}

Test the function.  Get some internal traces by setting 
time to~\(-1\) or~\(-2\)
\begin{matlab}
%}
disp('Using the following copy of patchSys2')
disp(which('patchSys2'))
U0 = rand([nx,ny,4,1,nPatch]);
U0t=patchSys2(0,U0(:),patches);
figure(1),hist(U0t(:))
disp('**** Finished first test evaluation')
%{
\end{matlab}



\paragraph{Determine spectrum}
Index~\verb|i| are the indices of patch-interior
displacements, and the number of unknowns is then its
length.
\begin{matlab}
%}
U0 = zeros([nx,ny,4,1,nPatch]);
U0([1 end],:,:) = nan;  
U0(:,[1 end],:) = nan;  
i = find(~isnan(U0));
nVariables = numel(i)
%{
\end{matlab}
Compute the Jacobian. 
\begin{matlab}
%}
disp("Computing Jacobian: wait ...")
tic  
Jac=nan(nVariables);
for j=1:nVariables
    U0(i)=((1:nVariables)'==j);
    Ut=patchSys2(1,U0(:),patches);
    Jac(:,j)=Ut(i);
end 
formJacTime=toc
%{
\end{matlab}
Explore structure and symmetry of Jacobian.
\begin{matlab}
%}
figure(2),clf,spy(Jac)
l1=2*(nx-2)*(ny-2);
l2=prod(nPatch);
JJ=reshape(Jac,l1,2,l2,l1,2,l2);
Jsw=reshape(JJ(:,2,:,:,1,:),l1*l2,l1*l2);
nonSymSW=norm(Jsw-Jsw')
Jse=reshape(JJ(:,2,:,:,2,:),l1*l2,l1*l2);
nonSymSE=norm(Jse-Jse')
%{
\end{matlab}
Compute the eigenvalues.
\begin{matlab}
%}
[evecs,evals]=eig(Jac,'vector');
%[~,j]=sort(real(evals)/(viscosity+1e-7)-abs(imag(evals)),"descend");
[~,j]=sort(abs(evals));
evals=evals(j);
evecs=evecs(:,j);
leadingEvals=evals(1:2:40)
%{
\end{matlab}
Plot spectra and export
\begin{matlab}
%}
figure(1),clf
plot(real(evals),imag(evals),".")
xlabel("$\Re\lambda$"), ylabel("$\Im\lambda$")
quasiLogAxes(0.001,0.1)
%{
\end{matlab}
Optionally export spectrum
\begin{matlab}
%}
if 0
  set(gca,'Position',[.2 .2 .6 .6])
  exportgraphics(gcf,basename+".pdf",'contenttype','vector')
end%if
%{
\end{matlab}


Fin.
%}

