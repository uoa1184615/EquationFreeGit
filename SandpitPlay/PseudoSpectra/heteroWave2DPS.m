%{ 
compute pseudospectra of ww2D to for damped linear ideal
wave system on staggered patches of staggered micro-grid in
2D. Minimal case takes a minute.  Here with microscale
heterogeneity of the bed, modelled within Cao2014b eqn (1). 
And with microscale heterogeneous drag (linear). AJR, May
2023
%}
clear all
global p N r dX micro lint b c
% p = patch micro-grid paramater, size is (4p-1)^2 p=1,2,...
% p=2 seems minimal heterogeneity, p=3 smoother variation
% (p,N)=(3,10) has 6075 variables, takes tens seconds
% (p,N)=(2,10) has 1875 variables, takes couple of seconds
p=2
% N = even number of patches in each dirn of periodic domain
% for spectral interpolation need N/2 to be odd
N=6
% r = fraction of distance a patch edge is to centre of next patch
r=1e-3*N   % so microperiod is 1/1000 of domain width
% overall magnitude of drag
cDrag = 1e-3
% dX = space step of the inter-patch distance of patch centres
dX=2*pi/N
X=(0.5:N)*dX; Y=X; % mid-point patches
% micro = index into locations of the PDEs on the micro-grid

% lint = 1 for linear interp, 0 for spectral
lint=0

n=4*p-1 % patch microgrid is nxn
% dx = space step of the micro-grid
dx=r*dX/(2*p-1)
xi=dx*(1-2*p:2*p-1)'; % local subgrid-variable
% locations on axes of micro-grid points in all patches 
[x,y]=ndgrid(xi,X);
x=x+y; y=shiftdim(x,-2);
% Maybe use double letters for 2D grids
% mesh of all locations, used for ICs
[xx,XX,yy,YY]=ndgrid(xi,X,xi,X);
xx=xx+XX; yy=yy+YY;

% heterogeneous microscale bottom shape is the following
bMean=0
microPeriod=r*dX
microAmp=0.1
b=-bMean+microAmp*cos(xi*2*pi/microPeriod).*cos(xi'*2*pi/microPeriod)
% heterogeneous microscale drag is random
c=2*cDrag*rand(n)
b=reshape(b,n,1,n);%reshape the bottom
c=reshape(c,n,1,n);%reshape the drag
% account for the optional macroscale modulation of the bed microscale
b=b+0*x+0*y;
sizeb=size(b)

% micro-grid indexes inside a pxp patch where
iMid=2*p; % mid-point index
iie=2:2:n-1; % interior even points
iio=3:2:n-2; % interior odd points (empty if p=1)
io=1:2:n; % odd points including edges.

% index of patch types: IH,IH are H-patches; IQ,IH are
% U-patches; IH,IQ are V-patches; and IQ,IQ are empty patches
IH=1:2:N; IQ=2:2:N;

% u(i,I,j,J) = (i,j)the microgrid value in (I,J)th patch
u=nan(n,N,n,N);
% Must set some values to find pattern of data in the 4D array
% relative to the equilibrium
h0=@(x,y) 0*x+0*y; %sin(x+y);
u0=@(x,y) 0; %sin(x);
v0=@(x,y) 0; %-cos(y);
% say h (e,H,e,H) (e,H,o,Q) (o,Q,e,H)
u(iie,IH,iie,IH)=h0(xx(iie,IH,iie,IH),yy(iie,IH,iie,IH))-b(iie,IH,iie,IH); 
u(iie,IH,iio,IQ)=h0(xx(iie,IH,iio,IQ),yy(iie,IH,iio,IQ))-b(iie,IH,iio,IQ);
u(iio,IQ,iie,IH)=h0(xx(iio,IQ,iie,IH),yy(iio,IQ,iie,IH))-b(iio,IQ,iie,IH);
u(iio,IQ,iio,IQ)=h0(xx(iio,IQ,iio,IQ),yy(iio,IQ,iio,IQ))-b(iio,IQ,iio,IQ);
% u field (o,H,e,H) (o,H,o,Q) (e,Q,e,H)
u(iio,IH,iie,IH)=u0(xx(iio,IH,iie,IH),yy(iio,IH,iie,IH)); 
u(iio,IH,iio,IQ)=u0(xx(iio,IH,iio,IQ),yy(iio,IH,iio,IQ));
u(iie,IQ,iie,IH)=u0(xx(iie,IQ,iie,IH),yy(iie,IQ,iie,IH));
u(iie,IQ,iio,IQ)=u0(xx(iie,IQ,iio,IQ),yy(iie,IQ,iio,IQ));
% v field (e,H,o,H) (e,H,e,Q) (o,Q,o,H)
u(iie,IH,iio,IH)=v0(xx(iie,IH,iio,IH),yy(iie,IH,iio,IH)); 
u(iie,IH,iie,IQ)=v0(xx(iie,IH,iie,IQ),yy(iie,IH,iie,IQ));
u(iio,IQ,iio,IH)=v0(xx(iio,IQ,iio,IH),yy(iio,IQ,iio,IH));
u(iio,IQ,iie,IQ)=v0(xx(iio,IQ,iie,IQ),yy(iio,IQ,iie,IQ));
% finally get the gather/scatter indices of micro-grid PDE points 
micro=find(~isnan(u));
nMicro=length(micro)
cf34Nsqnm2sq=0.75*N^2*(n-2)^2
huv=u(micro); 


%% Compute a time derivative as first check
ut=nan(n,N,n,N);
ut(micro)=hetero2Ddudt(0,huv);
normTest=norm(ut(micro))
if normTest>1e-10
nhtHH=norm(reshape(ut(iie,IH,iie,IH),[],1))
nhtHQ=norm(reshape(ut(iie,IH,iio,IQ),[],1))
nhtQH=norm(reshape(ut(iio,IQ,iie,IH),[],1))
nhtQQ=norm(reshape(ut(iio,IQ,iio,IQ),[],1))
nutHH=norm(reshape(ut(iio,IH,iie,IH),[],1))
nutHQ=norm(reshape(ut(iio,IH,iio,IQ),[],1))
nutQH=norm(reshape(ut(iie,IQ,iie,IH),[],1))
nutQQ=norm(reshape(ut(iie,IQ,iio,IQ),[],1))
nvtHH=norm(reshape(ut(iie,IH,iio,IH),[],1))
nvtHQ=norm(reshape(ut(iie,IH,iie,IQ),[],1))
nvtQH=norm(reshape(ut(iio,IQ,iio,IH),[],1))
nvtQQ=norm(reshape(ut(iio,IQ,iie,IQ),[],1))
end
assert(norm(ut(micro))<1e-10,'**** not an equilibrium')


%% Compute eigenvalues then plot in complex plane
% relative to the equilibrium identified above
delta=1 % perturb the surface by this much
Jac=nan(nMicro);
for j=1:nMicro
    Jac(:,j)=hetero2Ddudt(0,huv+delta*((1:nMicro)==j)')/delta;
end

figure(1)
figSc = inf % <- set figure scaling
multiscalePseudoSpectra(Jac,figSc,1,200)

figName = [mfilename num2str(round(abs(log10(figSc)),1)*10)]
set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 14 14] ...
      ,'renderer','Painters')
print('-depsc2',figName)
matlab2tikz([figName '.tex'],'showInfo',false ...
    ,'noSize',true,'parseStrings',false,'showWarnings',false ...
    ,'extraCode',['\tikzsetnextfilename{' figName '}'] ...
    ,'extraAxisOptions','\extraAxisOptions' ...
    ,'checkForUpdates',false)








function dhuvdt=hetero2Ddudt(t,huv)
% Compute time derivatives of ideal wave system on staggered
% patches of staggered micro-grid, grid0.  This version adds
% a spatially varying borrom b(x,y) via array b, and as
% modelled within Cao2014b eqn (1).  And also has a heterogeneous drag c.
% To interpolate to the edges via cross-patch values we have
% to include the 'ghost' patch to make four patches in every
% macro-cell.
% AJR, 7 Aug 2018 -- 13 Feb 2023
% I/O: t = time, ignored as autonomous
% huv = h,u,v values on staggered patches of staggered micro-grid
% dhuvdt = array of time derivatives of huv values
% Parameters (coupling has other global params)
% micro = index into locations of the PDEs on the micro-grid
% b = array of heterogeneous spatially varying bottom shape,
% defined at just all the h points.  Make the same variation
% in every patch of the same type, so 3D array with
% singleton for the 2nd dim.
global micro lint b c 

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
ut(iio,IQ,iio,IQ)=-(u(iio+1,IQ,iio,IQ)-u(iio-1,IQ,iio,IQ))/dx2 ...
                  -(u(iio,IQ,iio+1,IQ)-u(iio,IQ,iio-1,IQ))/dy2 ;
% u_t=-(h+b)_x at (o,H,e,H) (o,H,o,Q) (e,Q,e,H)
ut(iio,IH,iie,IH)=-(u(iio+1,IH,iie,IH)-u(iio-1,IH,iie,IH) ...
                   +b(iio+1,IH,iie,IH)-b(iio-1,IH,iie,IH))/dx2  ...
                   -c(iio,1,iie).*u(iio,IH,iie,IH); 
ut(iio,IH,iio,IQ)=-(u(iio+1,IH,iio,IQ)-u(iio-1,IH,iio,IQ) ...
                   +b(iio+1,IH,iio,IQ)-b(iio-1,IH,iio,IQ))/dx2  ...
                   -c(iio,1,iio).*u(iio,IH,iio,IQ);
ut(iie,IQ,iie,IH)=-(u(iie+1,IQ,iie,IH)-u(iie-1,IQ,iie,IH) ...
                   +b(iie+1,IQ,iie,IH)-b(iie-1,IQ,iie,IH))/dx2  ...
                   -c(iie,1,iie).*u(iie,IQ,iie,IH);
ut(iie,IQ,iio,IQ)=-(u(iie+1,IQ,iio,IQ)-u(iie-1,IQ,iio,IQ) ...
                   +b(iie+1,IQ,iio,IQ)-b(iie-1,IQ,iio,IQ))/dx2  ...
                   -c(iie,1,iio).*u(iie,IQ,iio,IQ);
% v_t=-(h+b)_y at (e,H,o,H) (e,H,e,Q) (o,Q,o,H)
ut(iie,IH,iio,IH)=-(u(iie,IH,iio+1,IH)-u(iie,IH,iio-1,IH) ...
                   +b(iie,IH,iio+1,IH)-b(iie,IH,iio-1,IH))/dy2  ...
                   -c(iie,1,iio).*u(iie,IH,iio,IH); 
ut(iie,IH,iie,IQ)=-(u(iie,IH,iie+1,IQ)-u(iie,IH,iie-1,IQ) ...
                   +b(iie,IH,iie+1,IQ)-b(iie,IH,iie-1,IQ))/dy2  ...
                   -c(iie,1,iie).*u(iie,IH,iie,IQ);
ut(iio,IQ,iio,IH)=-(u(iio,IQ,iio+1,IH)-u(iio,IQ,iio-1,IH) ...
                   +b(iio,IQ,iio+1,IH)-b(iio,IQ,iio-1,IH))/dy2  ...
                   -c(iio,1,iio).*u(iio,IQ,iio,IH);
ut(iio,IQ,iie,IQ)=-(u(iio,IQ,iie+1,IQ)-u(iio,IQ,iie-1,IQ) ...
                   +b(iio,IQ,iie+1,IQ)-b(iio,IQ,iie-1,IQ))/dy2  ...
                   -c(iio,1,iie).*u(iio,IQ,iie,IQ);
% extract and reshape
dhuvdt=ut(micro);
end%function






function [u,IH,IQ,iio,iie,dx]=stag2Dspectral(huv)
% Compute patch-edge values on staggered patches of
% staggered micro-grid-A.  Changed to interpolate the
% central-cross values, instead of the centre-point value. 
% See macroStagPatches.pdf
% AJR, 5 Sep 2018  -- 13 Feb 2023
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

% Interpolate spectrally to set all patch edge values.
% Interpolate central-cross values, including NaNs in the
% un-used microgrid points
% First interpolate in x-direction to edges
cHk=fft(u(iMid,IH,2:n-1,:),[],2);
cQk=fft(u(iMid,IQ,2:n-1,:),[],2);
% do we need some dXs in here? no, they cancel
kMax=(N/2-1)/2; %max wavenumber for interpolation, N/2 must be odd
% shift one patch to the right/up, in index 2
k1=pi/(N/2)*(rem(kMax+(0:N/2-1),N/2)-kMax);
krp=(1+r)*k1; % shift one + patch-half-width to the right/up
krm=(1-r)*k1; % shift one - patch-half-width to the right/up
u(n,IQ,2:n-1,:)=ifft(cHk.*exp(+1i*krp),[],2, 'symmetric');
u(1,IQ,2:n-1,:)=ifft(cHk.*exp(+1i*krm),[],2, 'symmetric');
u(n,IH,2:n-1,:)=ifft(cQk.*exp(-1i*krm),[],2, 'symmetric');
u(1,IH,2:n-1,:)=ifft(cQk.*exp(-1i*krp),[],2, 'symmetric');


% Second interpolate in y-direction to edges
cHk=fft(u(:,:,iMid,IH),[],4);
cQk=fft(u(:,:,iMid,IQ),[],4);
krp=shiftdim(krp,-2);
krm=shiftdim(krm,-2);
u(:,:,n,IQ)=ifft(cHk.*exp(+1i*krp),[],4, 'symmetric');
u(:,:,1,IQ)=ifft(cHk.*exp(+1i*krm),[],4, 'symmetric');
u(:,:,n,IH)=ifft(cQk.*exp(-1i*krm),[],4, 'symmetric');
u(:,:,1,IH)=ifft(cQk.*exp(-1i*krp),[],4, 'symmetric');

end%function
