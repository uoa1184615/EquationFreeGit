nPat=7
dx=0.2
nx=10
p=configPatches1(@sin,[0 9],'chebyshev',nPat,0,dx,nx,'EdgyInt',true,'nEdge',2);
x=squeeze(p.x)
