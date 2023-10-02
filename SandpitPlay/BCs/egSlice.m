[X,Y,Z] = meshgrid(-2:.2:2);%,-2:.3:2,-2:.4:2);
V = X.*exp(-X.^2-Y.^2-Z.^2);
xslice = [-1.2,0.8,2];   
yslice = [];
zslice = 0;
slice(X,Y,Z,V,xslice,yslice,zslice)
%xlabel('$x$'), ylabel('$y$'), zlabel('$z$')
