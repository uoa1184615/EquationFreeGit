%{
Plots relative error in the rate of macroscale evolution
when computing via projective integration od MidPoint rule.  Uses ODE
dx/dt=alpha*x with macroscale step Delta and microburst
delta.  Then aD=alpha*Delta and dD=delta/Delta.  For small
delta find aD \approx \pm\sqrt{-6*err}

Unresolved issue is the stability of fast modes when using
these projective integration.  Need to include some plot
from computation on dy/dt=-beta*y.  AJR Jan 2019

%}
aD=linspace(-1,1,20)*0.3;
dD=linspace(0,1,41)';

%%
figure(2)
serr=log(1+aD.*(1-dD).*exp(aD.*dD).*(1+aD/2.*(1-3*dD)))./aD-(1-dD);
serr(abs(imag(serr))>1e-8)=nan;
%serr(abs(serr)>1)=nan;
cs=[1;3]*10.^(-4:-1);
h=contour(aD,dD,real(serr),[-cs(:);0;cs(:)]);
clabel(h)
xlabel('\alpha\Delta'),ylabel('\delta/\Delta')
title('relative error in macroscale rate of PIRK2')

%%
figure(3), clf()
dD=reshape(linspace(0,0.3,21),[],1);
bD=exp( reshape(linspace(0,4,23),1,1,[]) );
ratio= exp(-aD-bD.*dD).*(1-(1-dD).*bD.*exp(-bD.*dD).*(1-bD/2.*(1-3*dD)));
f=@(r) (1-r).*(0.1-r).*(0.25-r).*(0.5-r);
serr=log(1+aD.*(1-dD).*exp(aD.*dD).*(1+aD/2.*(1-3*dD)))./aD-(1-dD);
serr(abs(imag(serr))>1e-8)=nan;
[faces,verts,colors] = isosurface(aD,dD,bD,f(ratio),0 ...
    ,0*log10(-real(serr))+0*bD+ratio); 
patch('Vertices',verts,'Faces',faces,'FaceVertexCData',colors,...
    'FaceColor','interp','EdgeColor','interp')
view(45,25), alpha 0.4
xlabel('\alpha\Delta'),ylabel('\delta/\Delta'),zlabel('\beta\Delta')
c=colorbar; c.Label.String='ratio';%'log10(macro-error)';
