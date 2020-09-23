function z=trySPMD
% for testGlobalParallel
global pat 
warning('**** beginning trySPMD function')
%pa=getfield(pat{labindex},'a')
%warning('**** done a simple assignment')
spmd 
warning('**** result of pat.y=pat.y-labindex')
pat.y=pat.y-labindex
nnn=nan(size(pat.y),codistributor1d(2))
warning('**** result of z=labindex*pat.y')
z=labindex*pat.y 
mpiprofile viewer
end
%warning('**** a simple assignment outside spmd fails')
%pata=pat.a
warning('**** ending trySPMD function')
end 