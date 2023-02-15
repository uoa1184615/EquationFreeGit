% Execute patch scheme on the heterogeneous miscroscale grid
% coded in function elastic2Dstaggered() with specified
% microscale boundary conditions on the beam ends.  AJR, 27
% Sep 2022 -- 4 Feb 2023
%!TEX root = doc.tex
%{
\section{\texttt{elastic2DstagHeteroSim}: simulate 2D
heterogeneous elastic patches on staggered grid}
\label{elastic2DstagHeteroSim}

Execute patch scheme on the miscroscale grid of
heterogeneous elasticity in a 2D beam as coded in the
function \verb|elastic2Dstaggered()| of
\cref{secelastic2Dstaggered}.  Here a beam of length~\(\pi\)
with either fixed or stress-free boundary conditions.  This
patch scheme provides an effective computational
homogenisation of the microscale heterogeneous beam. 
\begin{matlab}
%}
clear all
global patches 
%global OurCf2eps, OurCf2eps=true %option to save plot
%{
\end{matlab}

\paragraph{Set some patch and physical parameters}
\begin{matlab}
%}
nPatch = 5
nSubP = 7
ny = 4
xLim = [0 pi]
yLim= [-1 1]*0.06
nVars = 4*ny-2;
iOrd = 0
edgyInt = true
heteroStrength = 1 % -2--2 ??
stressFreeRightBC = true
viscosity = 2e-3; % 1e-3 sometimes has unstable micro-modes!!!
tEnd = 30
basename = mfilename;
%{
\end{matlab}
Viscosity of 1e-3 sometimes has unstable modes: only one
pair for homogeneous beam, but three pairs for heterogeneous
beam??


Set the cross-beam microscale \(y\)-coordinates and indices
in row vectors: one for \(u,\sigma_{xx},\sigma_{yy}\)-level;
and one for \(v,\sigma_{xy}\)-level.  Then show the
micro-grid spacing.
\begin{matlab}
%}
yv=linspace(yLim(1),yLim(2),ny);
yu=(yv(1:ny-1)+yv(2:ny))/2;
ju=1:ny-1; jv=1:ny;
dy=diff(yu(1:2))
%{
\end{matlab}



\paragraph{Physical elastic parameters---random} Initially
express in terms of Young's modulus and Poisson's ratio to
be passed to the micro-code via \verb|patches|. For
heterogeneous elasticity, in the microgrid, we need
\(\lambda\)~at just \nw-points, and \(\mu\)~at both
\nw,\se~points.  Might as well provide at all quarter points
as they both come from~\(E,\nu\), so we need these two
parameters at each of \(2n_y-1\) points across the beam.
Edgy interpolation requires patches to be two microgrid
points longer than any multiple of the microgrid
periodicity, whereas non-edgy requires one more than an even
multiple.  Have to be careful with heterogeneities as
\nw,\se-points are alternately either side of the notional
microgrid lines.  
\begin{matlab}
%}
xHeteroPeriod = (nSubP-1-edgyInt)/(2-edgyInt)
E  =  1*exp(heteroStrength*0.5*randn(xHeteroPeriod,2*ny-1)  );
nu = 0.3 +heteroStrength*0.1*(2*rand(xHeteroPeriod,2*ny-1)-1);
Equantiles  = quantile( E(:),[0:0.25:1])
nuQuantiles = quantile(nu(:),[0:0.25:1])
%{
\end{matlab}
Compute corresponding \(\lambda,\mu\)-fields, and store so 
\verb|configPatches1| \(x\)-expands micro-heterogeneity to 
fill patch as necessary.
\begin{matlab}
%}
lambda=nu.*E./(1+nu)./(1-2*nu);
mu    =    E./(1+nu)./2;
cElas = [mu lambda];
%{
\end{matlab}

\paragraph{Configure patches} Using \verb|hetCoeffs| to
manage form heterogeneity in the patch, and adjoin viscosity
coefficient.
\begin{matlab}
%}
Dom.type='equispace';
Dom.bcOffset = [ 0 0.75*stressFreeRightBC ];
configPatches1(@elastic2Dstaggered,xLim,'equispace',nPatch ...
    ,iOrd,dy,nSubP,'EdgyInt',edgyInt,'hetCoeffs',cElas);
patches.viscosity = viscosity;
patches.stressFreeRightBC = stressFreeRightBC;
dx=diff(patches.x(1:2))
%{
\end{matlab}


\paragraph{Set initial condition} Say choose initial
velocities zero, and displacements varying somehow.
\begin{matlab}
%}
U0 = zeros(nSubP,nVars,1,nPatch);
U0(:,ju,:,:) = 0.1*sin(patches.x/(1+stressFreeRightBC))+0.*yu;
if stressFreeRightBC
     U0(:,jv+ju(end),:,:) = 0.02*patches.x.^2.*(2.5-patches.x)+0.*yv;
else U0(:,jv+ju(end),:,:) = 0.2*sin(patches.x)+0.*yv;
end
%{
\end{matlab}
Optionally test the function
\begin{matlab}
%}
if 1 % test the new function
Ut=elastic2Dstaggered(-1,U0,patches);
Ut=reshape(Ut,nSubP,nVars,nPatch);
dudt=squeeze(Ut(:,ju,:))
dvdt=squeeze(Ut(:,jv+ju(end),:))
dutdt=squeeze(Ut(:,ju+2*ny-1,:))
dvtdt=squeeze(Ut(:,jv+2*ny-1+ju(end),:))
return
end
%{
\end{matlab}

\paragraph{Integrate in time}  
\begin{matlab}
%}
[ts,Us] = ode15s(@patchSys1, [0 tEnd], U0(:));
%{
\end{matlab}

\paragraph{Plot summary of simulation} First, subsample the
time-steps to roughly 150~steps. 
\begin{matlab}
%}
nt=numel(ts)
jt=1:ceil(nt/150):nt;
ts=ts(jt);
Us=reshape(Us(jt,:),[],nSubP,nVars,1,nPatch);
%{
\end{matlab}
Extract displacement fields.
\begin{matlab}
%}
us=Us(:,:,ju,:,:);
vs=Us(:,:,jv+ju(end),:,:);
%{
\end{matlab}
Form arrays of averages and variation across \(y\),
reshaped into 2D arrays.
\begin{matlab}
%}
uMean=reshape(mean(us,3),[],nSubP*nPatch); 
uStd=reshape(std(us,0,3),[],nSubP*nPatch); 
vMean=reshape(mean(vs,3),[],nSubP*nPatch); 
vStd=reshape(std(vs,0,3),[],nSubP*nPatch); 
%{
\end{matlab}
Plot \(xt\)-meshes coloured by the variation across~\(y\),
first the \(u\)-field of compression waves, see
\cref{figelastic2DstagHeteroSimu}.  The small space-shift in
\(x\)-location is due to the staggered micro-grid for the
patches as coded.
%\begin{figure}
%\centering
%\caption{\label{figelastic2DstagHeteroSimu}%
%mean field~\(u\) in space-time showing compression waves.  etc.}
%\input{Figs/elastic2DstagHeteroSimu}
%\end{figure}
\begin{matlab}
%}
xs = patches.x;  xs([1 end],:) = nan;
figure(1),clf
mesh(xs(:)+dx/4,ts,uMean,uStd);
ylabel('time $t$'),xlabel('space $x$'),zlabel('mean $u$')
xlim([0 pi]),view(40,55),colorbar
ifOurCf2eps([basename 'us'])%optionally save
%{
\end{matlab}
Draw the \(v\)-field showing the beam bending, see
\cref{figelastic2DstagHeteroSimv}.  The small space-shift in
\(x\)-location is due to the staggered micro-grid for the
patches as coded.
%\begin{figure}
%\centering
%\caption{\label{figelastic2DstagHeteroSimv}%
%mean field~\(v\) in space-time showing beam bending.  etc.}
%\input{Figs/elastic2DstagHeteroSimv}
%\end{figure}
\begin{matlab}
%}
figure(2),clf
mesh(xs(:)-dx/4,ts,vMean,vStd);
ylabel('time $t$'),xlabel('space $x$'),zlabel('mean $v$')
xlim([0 pi]),view(40,55),colorbar
drawnow
ifOurCf2eps([basename 'vs'])%optionally save
%{
\end{matlab}

\paragraph{Compute Jacobian and its spectrum} The aim here
is four-fold: \begin{itemize}
\item to show the patch scheme is stable for all elastic waves;
\item to have a clear separation between fast and slow waves; 
\item for compression macro-waves to match theory of
frequency \(\omega=\sqrt E k\) for integer wavenumber~\(k\);
and
\item for beam macro-waves to match beam theory of
\(mv_{tt}=EIv_{xxxx}\) for mass-density \(m:=2h\) and moment
\(I:=\int_{-h}^h y^2\,dy=\tfrac23h^3\) giving frequencies
\(\omega=\sqrt{E/3}hk\) where \(h\)~is the half-width of the
beam.
\end{itemize}

Form the Jacobian matrix, the linear operator, by numerical
construction about a zero field.  Use~\verb|i| to store the
indices of the micro-grid points that are interior to the
patches and hence are the system variables.
\begin{matlab}
%}
  u0 = zeros(nSubP,nVars,1,nPatch);
  u0([1 end],:,:,:)=nan; u0=u0(:);
  i=find(~isnan(u0));
  nVariables=length(i)
  Jac=nan(nVariables);
  for j=1:nVariables
    u0(i)=((1:nVariables)==j);
    dudt=patchSys1(0,u0);
    Jac(:,j)=dudt(i);
  end
  nJacEffZero = sum(abs(Jac(:))<1e-12)
  Jac(abs(Jac)<1e-12) = 0;
%{
\end{matlab}
Find the eigenvalues of the Jacobian, and list the 
slowest for inspection.
\begin{matlab}
%}
  [evecs,evals]=eig(Jac);  
  cabs=@(z) -real(z)*1e5+abs(imag(z));
  [~,j]=sort(cabs(diag(evals)));
  evals=diag(evals(j,j));  evecs=evecs(:,j);
  nZeroEvals=sum(abs(evals)<1e-5) 
  j=find(abs(evals)>1e-5);
  leadingNonzeroEvals=evals(j(1:4*nPatch))
%{
\end{matlab}
Plot quasi-log view of the spectrum, 
see \cref{figelastic2DstagHeteroSimEig}.
%\begin{figure}
%\centering
%\caption{\label{figelastic2DstagHeteroSimEig}%
%spectrum of the patch scheme showing slow macroscale waves, 
%and fast sub-patch microscale modes.  etc.}
%\input{Figs/elastic2DstagHeteroSimEig}
%\end{figure}
\begin{matlab}
%}
  figure(3),clf
  handle = plot(real(evals),imag(evals),'.');
  xlabel('real-part'), ylabel('imag-part')
  quasiLogAxes(handle,0.01,0.1);
  ifOurCf2eps([basename 'Eig'])%optionally save
%{
\end{matlab}

\input{elastic2Dstaggered.m}
%}

