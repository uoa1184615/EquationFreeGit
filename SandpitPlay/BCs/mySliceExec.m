% an example script to test slices()
x = -2:.2:2%; x(7)=nan
y = -1:.3:2
z = -2:.3:2
[X,Y,Z] = ndgrid(x,y,z);
V = X.*exp(-X.^2-0.5*Y.^2-Z.^2);
ix = [7,11,15]  
iy = length(y)
iz = [6 10]
clf, colormap(0.8*jet)
if 1, slices(x,y,z,V,ix,iy,iz)
else  slices(X,Y,Z,V,ix,iy,iz,@contour4)
end
axis equal, colorbar, view(-30,20)
xlabel('$x$'), ylabel('$y$'), zlabel('$z$')
if ~exist('OCTAVE_VERSION','builtin')
    exportgraphics(gcf,mfilename+".pdf")
end
