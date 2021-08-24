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
patch of half-size ratio~$0.4$ (relatively large for
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
patches = configPatches2(@mmNonDiffPDE,[-3 3 -2 2], nan ...
    , [9 7], 0, 0.4, nxy ,'EdgyInt',true);
patches.mmTime=1;
patches.Xlim=[-3 3 -2 2];
%{
\end{matlab}
The above two amendments to \verb|patches| should eventually be
part of the configuration function.

If we use more patches, then the algorithm goes berserk after some time??

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
u0 = u0.*(0.9+0.1*rand(size(u0)));
Nx = size(patches.x,5)
Ny = size(patches.y,6)
Npts = Nx*Ny;
D0 = zeros(2*Npts,1);
%{
\end{matlab}

Integrate in time to $t=4$?? using standard functions. In
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
OurCf2eps = false;
%{
\end{matlab}

Plot the movement of the mesh, the centre of each patch, as
a function of time: spatial domain horizontal, and time
vertical.
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
  pause(0.05)
end
%{
\end{matlab}

Fin.
%}