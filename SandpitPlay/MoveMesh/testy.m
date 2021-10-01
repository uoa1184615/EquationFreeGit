function testy(A,Is,lab)
return
for i=Is
  stdA = max(reshape(std(A,0,i),1,[]));
  if stdA>1e-10
    warning(lab)
    badArray = A
    stdA =stdA
    error(['failed invariance in index ' num2str(i)])
  end
end
end%function
