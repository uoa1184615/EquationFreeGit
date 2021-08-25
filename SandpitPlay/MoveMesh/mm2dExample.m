% mm2dExample: moving patches of Burgers PDE
%!TEX root = doc.tex
%{
\section{\texttt{mm2dExample}: example of moving
patches in 2D for nonlinear diffusion}
\label{sec:mm2dExample}



The code here shows one way to use moving patches in 2D.
However, \verb|mmPatchSys2()| has far too many ad hoc assumptions, 
so fix those before exploring predictions here.
%The simulation seems perfectly happy for the patches to move
%so that they overlap in the shock! and then separate again
%as the shock decays.

Establish global patch data struct to interface with a
function coding a nonlinear `diffusion' \pde: to be solved
on $6\times4$-periodic domain, with $9\times7$ patches,
spectral interpolation~($0$) couples the patches, each
patch of half-size ratio~$0.4$?? (relatively large for
visualisation), and with $5\times5$ points forming each
patch.  \cite{Roberts2011a} established that this scheme is
consistent with the \pde\ (as the patch spacing decreases).
Prefer EdgyInt as we suspect it performs better
for moving meshes.
\begin{matlab}
%}
clear all
global patches
nxy=5
Nx=9, Ny=7
patches = configPatches2(@mmNonDiffPDE,[-3 3 -2 2], nan ...
    , [Nx Ny], 0, 0.2, nxy ,'EdgyInt',true);
patches.mmTime=1;
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


\paragraph{Simulate in time}
Set an  initial condition of a perturbed-Gaussian using
auto-replication of the spatial grid.
\begin{matlab}
%}
u0 = exp(-patches.x.^2-patches.y.^2);
u0 = u0.*(0.9+0.0*rand(size(u0))) +0.001;
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
Choose whether to save some plots, or not.
\begin{matlab}
%}
global OurCf2eps
OurCf2eps = true;
%{
\end{matlab}

Plot the movement of the mesh, and the field vertical, at the centre of each patch.
\begin{matlab}
%}
nTime=length(ts);
Ds=reshape(us(:,1:2*Npts).',Nx,Ny,2,nTime);
us=reshape(us(:,2*Npts+1:end).',nxy,nxy,Nx,Ny,nTime);
Us=shiftdim( mean(mean(us,1),2) ,2);
%% section marker for plot execution
figure(1),clf, colormap(0.8*hsv)
Xs=shiftdim(mean(patches.x),4);
Ys=shiftdim(mean(patches.y),4);
for k=1:nTime
  Xk=Xs+Ds(:,:,1,k);
  Yk=Ys+Ds(:,:,2,k);
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
  if rem(k,31)==1, ifOurCf2eps([mfilename num2str(k)]), end
  pause(0.05)
end
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




\paragraph{Spectrum of the moving patch system}
Compute the spectrum based upon the linearisation about some
state: \(u={}\)constant with \(D=0\) are equilibria;
otherwise the computation is about a 'quasi-equilibrium' on
the `fast-time'.
\begin{matlab}
%}
u00 = 0.1
u0 = u00+0.1*exp(-patches.x.^2-patches.y.^2);
u0([1 end],:,:,:)=nan;  u0(:,[1 end],:,:)=nan;
u0 = [zeros(2*Npts,1); u0(:)];
f0 = mmPatchSys2(0,u0);
normf0=norm(f0)
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
(\cref{fig:mm2dExample}).
\begin{matlab}
%}
eval=eig(Jac);
k=find(abs(imag(eval))<1e-6);
eval(k)=real(eval(k));
[~,k]=sort(-real(eval));
eval=eval(k);
nZero = sum(abs(real(eval))<1e-6)
nSlow = sum(-3*u00^2*300<real(eval))-nZero
eSlow = eval(nZero+(1:2:nSlow))
eFast = eval([nZero+nSlow+1 end])
figure(3),clf
plot(asinh(real(eval)),asinh(imag(eval)),'.')
xlabel('Re\lambda'), ylabel('Im\lambda')
ticks=[1;2;5]*10.^(0:6);
ticks=sort([0;ticks(:);-ticks(:)]);
set(gca,'Xtick',asinh(ticks) ...
    ,'XtickLabel',cellstr(num2str(ticks,4)) ...
    ,'XTickLabelRotation',30)
set(gca,'Ytick',asinh(ticks) ...
    ,'YtickLabel',cellstr(num2str(ticks,4)))
grid
ifOurCf2eps([mfilename 'Spec'])
%{
\end{matlab}
\begin{figure}
\centering \caption{\label{fig:mm2dExample}spectrum
of the moving mesh 2D diffusion system (about \(u=0.1\)).  The
clusters are: right real, macroscale diffusion modes with some neutral mesh deformations;
left real, moving mesh and sub-patch modes.}
\includegraphics[scale=0.85]{Figs/mm2dExampleSpec}
\end{figure}

Fin.
%}