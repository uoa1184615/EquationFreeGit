%{
An Issue with SPMD
https://au.mathworks.com/matlabcentral/answers/566349-an-issue-with-spmd?s_tid=answers_rc1-1_p1_BOTH
%}
parpool('local', 4)
spmd
    x = magic(labindex); % x is type 'double'
end
class(x) % gets 'Composite'
size(x) % gets [1,4]
size(x{2}) % gets [2,2] - the element from worker with labindex==2
spmd
    class(x) % on the worker, we see the 'double'
    size(x) % gets [labindex,labindex] on each worker
end