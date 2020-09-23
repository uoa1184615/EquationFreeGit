% testing global parallel, AJR Sep 2020
%{
Let Lab denote the index of a pcu in the parallel pool.

An ordinary variable assigned inside spmd acquires a set of
values, one for each Lab, whether global or not.  Outside
spmd can access values via {Lab}.  Called a "composite".

For a "distributed" variable: outside spmd composed and
returns full array of values;  inside spmd composed of only
part of array, and returns that part for each Lab.

However, upon "delete(ppool)" the composites and the
distributed all lose their values because the parallel pool
has closed/shut-down.

Appears can use distributed array in a field of a global
struct provided the field is codistributed, and one only
refers to it inside "spmd ... end" pairs.

Further, because the global struct in its entirety is copied
for each worker, once you include a distributed variable to
the struct, then the entire struct becomes 'composite', and
so one can only thereafter refer to all fields in the struct
from within "spmd ... end" pairs.  If outside spmd, then one
has to subsequently use functions like "getfield(pat{1},'a')".

From the Matlab documentation: The body of an spmd statement
cannot contain global or persistent variable declarations.
The reason is that these variables are not synchronized
between workers. You can use global or persistent variables
within functions, but their value is only visible to the
worker that creates them. Seems to  =>  if using global then
have to start and stop spmd environment inside each
function, and stop-start around each function call.

From Matlab Doc---Nested spmd Statements:  The body of an
spmd statement cannot directly contain another spmd.
However, it can call a function that contains another spmd
statement. The inner spmd statement does not run in parallel
in another parallel pool, but runs serially in a single
thread on the worker running its containing function.

Appears setting distributed field outside spmd DNW.
%}
clear all
global pat c
c=rand
ppool=gcp
addAttachedFiles(ppool,{'tryPlain.m'})
%spmd, % cannot declare global inside spmd blocks (but maybe in functions)
pat.a=rand
%pat.y=11:18 %meaningless assignation just to test
% but make sure clear such vars to ensure different parallelisms proceed?
% Appears setting distributed outside spmd DNW
%pat.z=distributed(11:16)
warning('**** initialising pat.y')
spmd, mpiprofile on,pat.y=codistributed(1:6,codistributor1d(2)), end

%trySPMD

%spmd
tryPlain(1,2,3)
%end
