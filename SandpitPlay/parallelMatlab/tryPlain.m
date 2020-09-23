function z=tryPlain(a,b,c)
%{
for testGlobalParallel: because this function is executed
within spmd...end, then the global variable pat is not
accessible, and so here is empty valued
%}
nArgs= nargin
if nargin<3, global c, end
cVal=c
global pat 
warning('**** tryPlain gpat=pat')
gpat=pat
warning('**** tryPlain pat.y=pat.y-labindex')
pat.y=pat.y-labindex
warning('**** tryPlain z=labindex*pat.y')
z=labindex*pat.y 
warning('**** ending tryPlain function')
end 