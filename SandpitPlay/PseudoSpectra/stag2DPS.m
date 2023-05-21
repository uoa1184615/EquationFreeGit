% compute pseudospectra of stag2D to for ideal wave system
% on staggered patches of staggered micro-grid in 2D. 
% Minimal case takes a minute.  Pseudospectrum seems OK.
% AJR, May 2023
global p N r dX micro lint
% p = patch micro-grid paramater, size is (4p-1)^2 p=1,2,...
p=2
% N = even number of patches in each dirn of periodic domain
% for spectral interpolation need N/2 to be odd
N=6
% r = fraction of distance a patch edge is to centre of next patch
r=0.1
% dX = space step of the inter-patch distance of patch centres
dX=2*pi/N
X=(0.5:N)*dX; Y=X; % mid-point patches
% micro = index into locations of the PDEs on the micro-grid

% lint = 1 for linear interp, 0 for spectral
lint=0

% dx = space step of the micro-grid
dx=r*dX/(2*p-1)
% locations on axes of micro-grid points in all patches 
[x,y]=ndgrid(dx*(1-2*p:2*p-1),X);
x=x+y; y=x;
% Maybe use double letters for 2D grids
% mesh of all locations, used for ICs
[xx,XX,yy,YY]=ndgrid(dx*(1-2*p:2*p-1),X,dx*(1-2*p:2*p-1),X);
xx=xx+XX; yy=yy+YY;

% micro-grid indexes inside a pxp patch where
n=4*p-1 % patch microgrid is nxn
iMid=2*p; % mid-point index
iie=2:2:n-1; % interior even points
iio=3:2:n-2; % interior odd points (empty if p=1)
%io=1:2:n; % odd points including edges.

% index of patch types: IH,IH are H-patches; IQ,IH are
% U-patches; IH,IQ are V-patches; and IQ,IQ are empty patches
IH=1:2:N; IQ=2:2:N;

% u(i,I,j,J) = (i,j)the microgrid value in (I,J)th patch
u=nan(n,N,n,N);
% Must set some values to find pattern of data in the 4D array
h0=@(x,y) 1; %sin(x+y);
u0=@(x,y) 0; %sin(x);
v0=@(x,y) 0; %-cos(y);
% say h (e,H,e,H) (e,H,o,Q) (o,Q,e,H)
u(iie,IH,iie,IH)=h0(xx(iie,IH,iie,IH),yy(iie,IH,iie,IH)); 
u(iie,IH,iio,IQ)=h0(xx(iie,IH,iio,IQ),yy(iie,IH,iio,IQ));
u(iio,IQ,iie,IH)=h0(xx(iio,IQ,iie,IH),yy(iio,IQ,iie,IH));
% u field (o,H,e,H) (o,H,o,Q) (e,Q,e,H)
u(iio,IH,iie,IH)=u0(xx(iio,IH,iie,IH),yy(iio,IH,iie,IH)); 
u(iio,IH,iio,IQ)=u0(xx(iio,IH,iio,IQ),yy(iio,IH,iio,IQ));
u(iie,IQ,iie,IH)=u0(xx(iie,IQ,iie,IH),yy(iie,IQ,iie,IH));
% v field (e,H,o,H) (e,H,e,Q) (o,Q,o,H)
u(iie,IH,iio,IH)=v0(xx(iie,IH,iio,IH),yy(iie,IH,iio,IH)); 
u(iie,IH,iie,IQ)=v0(xx(iie,IH,iie,IQ),yy(iie,IH,iie,IQ));
u(iio,IQ,iio,IH)=v0(xx(iio,IQ,iio,IH),yy(iio,IQ,iio,IH));
% finally get the gather/scatter indices of micro-grid PDE points 
micro=find(~isnan(u));
nMicro=length(micro)
huv=u(micro); 


%% Compute a time derivative as quasi-check
%ut=nan(n,N,n,N);
%ut(micro)=stag2Ddudt(0,huv);


%% Compute eigenvalues then plot in complex plane
delta=1 % perturb the surface by this much
huv=0*huv; % about zero, assumed equilibrium
Jac=nan(nMicro);
for j=1:nMicro
    Jac(:,j)=stag2Ddudt(0,huv+delta*((1:nMicro)==j)')/delta;
end
figure(1)
multiscalePseudoSpectra(Jac,1,1,200)
figfn=mfilename;
set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 14 14] ...
      ,'renderer','Painters')
print('-depsc2',figfn)
matlab2tikz([figfn '.tex'],'showInfo',false ...
    ,'noSize',true,'parseStrings',false,'showWarnings',false ...
    ,'extraCode',['\tikzsetnextfilename{' figfn '}'] ...
    ,'extraAxisOptions','\extraAxisOptions' ...
    ,'checkForUpdates',false)






function dhuvdt=stag2Ddudt(t,huv)
% Compute time derivatives of ideal wave system on staggered
% patches of staggered micro-grid, grid0.  AJR, 7 Aug 2018
% I/O: t = time, ignored as autonomous
% huv = h,u,v values on staggered patches of staggered micro-grid
% dhuvdt = array of time derivatives of huv values
% Parameters (coupling has other global params)
% micro = index into locations of the PDEs on the micro-grid
global micro lint

% u(i,I,j,J) = (i,j)the microgrid value in (I,J)th patch
% Specify the function that couples patches by
% providing patch edge values through  whatever
% mechanism
if lint % case 1 do linear interp
[u,IH,IQ,iio,iie,dx]=stag2Dcouple(huv);
else % case 0 do spectral interp
[u,IH,IQ,iio,iie,dx]=stag2Dspectral(huv);
end

%% Compute time derivatives in 3 different patch types
% dx = space step of the micro-grid
dx2=2*dx; dy2=dx2;
ut=nan(size(u));
% h_t=-u_x-v_y at (e,H,e,H) (e,H,o,Q) (o,Q,e,H)
ut(iie,IH,iie,IH)=-(u(iie+1,IH,iie,IH)-u(iie-1,IH,iie,IH))/dx2 ...
                  -(u(iie,IH,iie+1,IH)-u(iie,IH,iie-1,IH))/dy2 ; 
ut(iie,IH,iio,IQ)=-(u(iie+1,IH,iio,IQ)-u(iie-1,IH,iio,IQ))/dx2 ...
                  -(u(iie,IH,iio+1,IQ)-u(iie,IH,iio-1,IQ))/dy2 ;
ut(iio,IQ,iie,IH)=-(u(iio+1,IQ,iie,IH)-u(iio-1,IQ,iie,IH))/dx2 ...
                  -(u(iio,IQ,iie+1,IH)-u(iio,IQ,iie-1,IH))/dy2 ;
% u_t=-h_x at (o,H,e,H) (o,H,o,Q) (e,Q,e,H)
ut(iio,IH,iie,IH)=-(u(iio+1,IH,iie,IH)-u(iio-1,IH,iie,IH))/dx2 ; 
ut(iio,IH,iio,IQ)=-(u(iio+1,IH,iio,IQ)-u(iio-1,IH,iio,IQ))/dx2 ;
ut(iie,IQ,iie,IH)=-(u(iie+1,IQ,iie,IH)-u(iie-1,IQ,iie,IH))/dx2 ;
% v_t=-h_y at (e,H,o,H) (e,H,e,Q) (o,Q,o,H)
ut(iie,IH,iio,IH)=-(u(iie,IH,iio+1,IH)-u(iie,IH,iio-1,IH))/dy2 ; 
ut(iie,IH,iie,IQ)=-(u(iie,IH,iie+1,IQ)-u(iie,IH,iie-1,IQ))/dy2 ;
ut(iio,IQ,iio,IH)=-(u(iio,IQ,iio+1,IH)-u(iio,IQ,iio-1,IH))/dy2 ;
% extract and reshape
dhuvdt=ut(micro);
end





function [u,IH,IQ,iio,iie,dx]=stag2Dspectral(huv)
% Compute patch-edge values on staggered patches of
% staggered micro-grid-A.  AJR, 5 Sep 2018
% I/O: huv = vector of h,u,v values on staggered
%     patches of staggered micro-grid
% u = 4D array of huv values, with patch-edge values interpolated
% IH = H-patch indices
% IQ = U,V-patch indices, in conjunction with IH
% iio = indices of internal lattice, odd-points
% iie = indices of internal lattice, even-points
% dx = microscale lattice spacing
% Parameters
% p = patch micro-grid paramater, size is (4p-1)^2 p=1,2,...
% N = even number of patches in each dirn of periodic domain
% r = fraction of distance a patch edge is to centre of next patch
% dX = space step of the inter-patch distance of patch centres
% micro = index into locations of the PDEs on the micro-grid
global p N r dX micro

% dx = space step of the micro-grid
dx=r*dX/(2*p-1);

% micro-grid indexes inside a nxn patch where
n=4*p-1; % patch microgrid is nxn
iMid=2*p; % mid-point index
iie=2:2:n-1; % interior even points
iio=3:2:n-2; % interior odd points (empty when p=1)

% index of patch types: IH,IH are H-patches; IQ,IH are
% U-patches; IH,IQ are V-patches; and IQ,IQ are empty patches
IH=1:2:N; IQ=2:2:N;

% u(i,I,j,J) = (i,j)the microgrid value in (I,J)th patch
u=nan(n,N,n,N);
u(micro)=huv;

%% Interpolate spectrally to set necessary patch edge values
% Omits the corner points on all patches
% N/2, 2D Fourier coeffs of the three types of patches
cHk=fft2(squeeze(u(iMid,IH,iMid,IH))); 
cUk=fft2(squeeze(u(iMid,IQ,iMid,IH))); 
cVk=fft2(squeeze(u(iMid,IH,iMid,IQ))); 
% do we need some dXs in here? no, they cancel
kMax=(N/2-1)/2; %max wavenumber for interpolation, N/2 must be odd
% shift one patch to the right/up
k1=pi/(N/2)*(rem(kMax+(0:N/2-1),N/2)-kMax);
kr=r*k1; % shift patch-half-width to the right/up
krp=(1+r)*k1; % shift one + patch-half-width to the right/up
krm=(1-r)*k1; % shift one - patch-half-width to the right/up

for j=iie
  ks=(j-iMid)*2/(n-1)*kr;
  % first interpolate U and V values to H-patches
  u(n,IH,j,IH)=ifft2(cUk.*exp(1i*( -krm'+ks )));
  u(1,IH,j,IH)=ifft2(cUk.*exp(1i*( -krp'+ks )));
  u(j,IH,n,IH)=ifft2(cVk.*exp(1i*( -krm +ks')));
  u(j,IH,1,IH)=ifft2(cVk.*exp(1i*( -krp +ks')));
  % second interpolate H to U-patches
  u(n,IQ,j,IH)=ifft2(cHk.*exp(1i*( +krp'+ks )));
  u(1,IQ,j,IH)=ifft2(cHk.*exp(1i*( +krm'+ks )));
  % third interpolate H to V-patches
  u(j,IH,n,IQ)=ifft2(cHk.*exp(1i*( +krp +ks')));
  u(j,IH,1,IQ)=ifft2(cHk.*exp(1i*( +krm +ks')));
end
for j=1:2:n%iio
  ks=(j-iMid)*2/(n-1)*kr;
  % second interpolate V to U-patches
  u(j,IQ,n,IH)=ifft2(cVk.*exp(1i*( k1'+ks'-krm )));
  u(j,IQ,1,IH)=ifft2(cVk.*exp(1i*( k1'+ks'-krp )));
  % third interpolate U to V-patches
  u(n,IH,j,IQ)=ifft2(cUk.*exp(1i*( k1 +ks -krm')));
  u(1,IH,j,IQ)=ifft2(cUk.*exp(1i*( k1 +ks -krp')));
  % also do top/bot u of V-patches
  u(j,IH,n,IQ)=ifft2(cUk.*exp(1i*( krp +ks' -k1')));
  u(j,IH,1,IQ)=ifft2(cUk.*exp(1i*( krm +ks' -k1')));
  % and left/right v of U-patches
  u(n,IQ,j,IH)=ifft2(cVk.*exp(1i*( krp' +ks -k1 )));
  u(1,IQ,j,IH)=ifft2(cVk.*exp(1i*( krm' +ks -k1 )));
end
end

