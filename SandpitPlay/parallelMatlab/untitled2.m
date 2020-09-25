clear all
gcp('nocreate')
q=distributed(zeros(6))
spmd
  r=nan([4 size(q)]);
  for i=1:4, r(i,:,:)=gather((q+i)*labindex); end
end
classr=class(r)
r=r{1}
%delete(gcp)
