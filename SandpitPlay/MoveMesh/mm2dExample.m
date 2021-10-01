% mm2dExample: moving patches of Burgers PDE
%!TEX root = doc.tex
%{
\section{\texttt{mm2dExample}: example of moving
patches in 2D for nonlinear diffusion}
\label{sec:mm2dExample}



The code here shows two ways to use moving patches in 2D.
Plausible generalisations from the 1D code to this 2D code is the case \verb|adhoc|.
The alternative \verb|Huang98| aims to implement the method of \cite{Huang98}.
\begin{matlab}
%}
clear all
global theMethod
if 0, theMethod = 'adhoc',
else  theMethod = 'Huang98', end
%{
\end{matlab}
However, \verb|mmPatchSys2()| has far too many ad hoc assumptions, 
so fix those before exploring predictions here.

Establish global patch data struct to interface with a
function coding the microscale.
Prefer EdgyInt as we suspect it performs better
for moving meshes.
Using \verb|nxy=3| means that there are no sub-patch modes, all modes are those of the macro-diffusion and the macro-mesh movement.  
There are \(N_x+N_y\) zero eigenvalues associated with the mesh movement.  
And there are \(N_xN_y\) slow eigenvalues of the diffusion (one of them looks like zero badly affected by round-off to be as big as \(10^{-4}\) or so).  
So far we generally see the macro-diffusion is poorly perturbed by the mesh movement in that there are some diffusion modes with imaginary part up to five.  
\begin{matlab}
%}
global patches
nxy=3 % =3 means no sub-patch dynamics
Nx=7, Ny=5
patches = configPatches2(@mmNonDiffPDE,[-3 3 -2 2], nan ...
    , [Nx Ny], 0, 0.1, nxy ,'EdgyInt',true);
patches.mmTime=0.03;
patches.Xlim=[-3 3 -2 2];
Npts = Nx*Ny;
%{
\end{matlab}
The above two amendments to \verb|patches| should eventually be
part of the configuration function.


\paragraph{Decide the moving mesh time parameter}
%Here for \(\epsilon=0.02\).
%\begin{itemize}
%\item Would be best if the moving mesh was no stiffer than
%the stiffest microscale sub-patch mode.  These would both be
%the zig-zag modes.
%\begin{itemize}
%\item Here the mesh \pde\ is \(X_t=(N^2/\tau)X_{jj}\) so its
%zig-zag mode decays with rate \(4N^2/\tau\)\,.
%\item Here the patch width is~\(h=0.2/15=1/75\), and so the
%microscale step is \(\delta=h/4=1/300\).  Hence the
%diffusion \(u_t=\epsilon u_{xx}\) has zig-zag mode decaying
%at rate \(4\epsilon/\delta^2\).
%\end{itemize}
%So, surely best to have \(4N^2/\tau \lesssim 4\epsilon
%/\delta^2\), that is, \(\tau \gtrsim N^2\delta^2 /\epsilon
%\approx 0.1\).
%
%\item But also we do not want the slowest modes of the
%moving mesh to obfuscate the system's macroscale modes---the
%macroscale zig-zag.
%\begin{itemize}
%\item The slowest moving mesh mode has wavenumber in~\(j\)
%of~\(2\pi/N\), and hence rate of decay \((N^2/\tau)
%(2\pi/N)^2 =4\pi^2/\tau\).
%\item The fastest zig-zag mode of the system \(U_t=\epsilon
%U_{xx}\) on step~\(H\) has decay rate \(4\epsilon/H^2\).
%\end{itemize}
%So best if \(4\pi^2/\tau \gtrsim 4\epsilon/H^2\), that is,
%\(\tau \lesssim \pi^2H^2 /\epsilon \approx 2\)\,.   
%
%(Computations indicate need \(\tau<0.8\)??)
%\end{itemize}




\paragraph{Spectrum of the moving patch system}
Compute the spectrum based upon the linearisation about some
state: \(u={}\)constant with \(D=0\) are equilibria;
otherwise the computation is about a 'quasi-equilibrium' on
the `fast-time'.
\begin{matlab}
%}
global ind, ind=2
evals=[];
patches.mmTime = patches.mmTime/0.95;
for iv=1:4
  patches.mmTime = 0.95*patches.mmTime;
%
u00 = 0.1
u0 = u00+sin(0*patches.x*pi/3+0*patches.y*pi/2);
u0([1 end],:,:)=nan;  u0(:,[1 end],:)=nan;
u0 = [zeros(2*Npts,1); u0(:)];
f0 = mmPatchSys2(0,u0);
normf0 = norm(f0)
%if normf0>1e-9, error('Jacobian: u0 is not equilibrium'),end
%{
\end{matlab}
But we must only use the dynamic variables, so let's find
where they are. 
\begin{matlab}
%}
i=find(~isnan( u0(:) )); 
nJac=length(i)
%{
\end{matlab}
Construct Jacobian with numerical differentiation.
\begin{matlab}
%}
deltau=1e-7;
Jac=nan(nJac);
for j=1:nJac
    uj=u0; uj(i(j))=uj(i(j))+deltau;
    fj = mmPatchSys2(0,uj);
    Jac(:,j)=(fj(i)-f0(i))/deltau;
end
%{
\end{matlab}
Compute and plot the spectrum with non-linear axis scaling
(\cref{fig:mm2dExampleSpec}).
\begin{matlab}
%}
eval=eig(Jac);
[~,k]=sort(-real(eval));
eval=eval(k);
nZero = sum(abs(real(eval))<1e-3)
nSlow = sum(-100<real(eval))-nZero
%eSlowest = eval(1:30) %(0+(1:2:nSlow))
%eFast = eval([nZero+nSlow+1 end])
evals=[evals eval];
end%iv-loop
%{
\end{matlab}

\paragraph{Plot spectrum}
Choose whether to save some plots, or not.
\begin{matlab}
%}
global OurCf2eps
OurCf2eps = false;
%{
\end{matlab}
Draw spectrum on quasi-log axes.
\begin{matlab}
%}
figure(3),clf
hp = plot(real(evals),imag(evals),'.');
xlabel('Re\lambda'), ylabel('Im\lambda')
quasiLogAxes(hp,1,1);
ifOurCf2eps([mfilename theMethod 'Spec'])
return%%%%%%%%%%%%%%
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:mm2dExampleadhocSpec}spectrum
of the \emph{adhoc} moving mesh 2D diffusion system (about \(u=0.1\)).  The
clusters are: right real, macroscale diffusion modes with some neutral mesh deformations;
left real, moving mesh and sub-patch modes.  Coloured dots are `trails' for \(\tau\)~reducing by~5\% between each, starting from blue dots.}
\includegraphics[scale=0.85]{Figs/mm2dExampleadhocSpec}
\end{figure}
\begin{figure}
\centering \caption{\label{fig:mm2dExampleHuang98Spec}spectrum
of the \emph{Huang98} moving mesh 2D diffusion system (about \(u=0.1\)). 
Currently there are badly unstable modes. The
clusters are: \dotfill\ confused.}
\includegraphics[scale=0.85]{Figs/mm2dExampleHuang98Spec}
\end{figure}


\paragraph{Simulate in time}
Set an  initial condition of a perturbed-Gaussian using
auto-replication of the spatial grid.
\begin{matlab}
%}
u0 = exp(-patches.x.^2-patches.y.^2);
u0 = 1+0*u0.*(0.9+0.0*rand(size(u0)));
D0 = zeros(2*Npts,1);
%{
\end{matlab}

Integrate in time to $t=2$ using standard functions. In
\Matlab\ \verb|ode15s| would be natural as the patch scheme
is naturally stiff, but \verb|ode23| is quicker \cite
[Fig.~4] {Maclean2020a}.  Ask for output at non-uniform
times because the diffusion slows.
\begin{matlab}
%}
disp('Simulating nonlinear diffusion h_t=(h^3)_xx+(h^3)_yy')
tic
[ts,us] = ode23(@mmPatchSys2,2*linspace(0,1).^2,[D0;u0(:)]);
cpuTime = toc
%{
\end{matlab}

\paragraph{Plots}
Extract data from time simulation.  Be wary that the patch-edge values do not change from initial, so either set to~\verb|NaN|, or set via interpolation.
\begin{matlab}
%}
nTime=length(ts);
Ds=reshape(us(:,1:2*Npts).',1,1,Nx,Ny,2,nTime);
us=reshape(us(:,2*Npts+1:end).',nxy,nxy,Nx,Ny,nTime);
us([1 end],:,:,:)=nan; us(:,[1 end],:,:)=nan; % nan edges
%{
\end{matlab}

Choose macro-mesh plot or micro-surf-patch plots. 
\begin{matlab}
%}
if 1
%{
\end{matlab}
Plot the movement of the mesh, with the field vertical, at
the centre of each patch.  
\begin{matlab}
%}
%% section marker for macro-mesh plot execution
figure(1),clf, colormap(0.8*hsv)
Us=shiftdim( mean(mean(us,1,'omitnan'),2,'omitnan') ,2);
Xs=shiftdim(mean(patches.x),4);
Ys=shiftdim(mean(patches.y),4);
for k=1:nTime
  Xk=Xs+shiftdim(Ds(:,:,:,:,1,k),2);
  Yk=Ys+shiftdim(Ds(:,:,:,:,2,k),2);
  if k==1,
    hand=mesh(Xk,Yk,Us(:,:,k));
    ylabel('space y'),xlabel('space x'),zlabel('mean field U')
    axis([patches.Xlim 0 1]), caxis([0 1])
    colorbar
    if 0, view(0,90) % vertical view
    else  view(-25,60) % 3D perspective
    end  
  else
    set(hand,'XData',Xk,'YData',Yk ...
       ,'ZData',Us(:,:,k),'CData',Us(:,:,k))
  end
  legend(['time =' num2str(ts(k),4)],'Location','north')
  if rem(k,31)==1, ifOurCf2eps([mfilename theMethod num2str(k)]), end
  pause(0.05)
end% for each time
else%if macro-mesh or micro-surf
%{
\end{matlab}
Plot the movement of the patches, with the field vertical in 
each patch.
\begin{matlab}
%}
%% section marker for patch-surf plot execution
figure(2),clf, colormap(0.8*hsv)
xs=reshape(patches.x,nxy,1,Nx,1);
ys=reshape(patches.y,1,nxy,1,Ny);
for k=1:nTime
  xk=xs+0*ys+Ds(:,:,:,:,1,k);
  yk=ys+0*xs+Ds(:,:,:,:,2,k);
  uk=reshape(permute(us(:,:,:,:,k),[1 3 2 4]),nxy*Nx,nxy*Ny);
  xk=reshape(permute(xk,[1 3 2 4]),nxy*Nx,nxy*Ny);
  yk=reshape(permute(yk,[1 3 2 4]),nxy*Nx,nxy*Ny);
  if k==1,
    hand=surf(xk,yk,uk);
    ylabel('space y'),xlabel('space x'),zlabel('field u(x,y,t)')
    axis([patches.Xlim 0 1]), caxis([0 1])
    colorbar
  else
    set(hand,'XData',xk,'YData',yk,'ZData',uk,'CData',uk)
  end
  legend(['time =' num2str(ts(k),4)],'Location','north')
%  if rem(k,31)==1, ifOurCf2eps([mfilename theMethod num2str(k)]), end
  pause(0.05)
end% for each time
%%
end%if macro-mesh or micro-surf
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:mm2dExample}field
$u(x,t)$ of the moving patch scheme applied to nonlinear diffusion~\pde.}
\begin{tabular}{@{}cc@{}}
\includegraphics[width=0.47\linewidth]{Figs/mm2dExample1}
&
\includegraphics[width=0.47\linewidth]{Figs/mm2dExample32}
\\
\includegraphics[width=0.47\linewidth]{Figs/mm2dExample63}
&
\includegraphics[width=0.47\linewidth]{Figs/mm2dExample94}
\end{tabular}
\end{figure}




Fin.
%}