%PIG implements Projective Integration scheme with any user-specified
%solver for the coarse variables, and a user-specified microsolver.
%JM, September 2018.
%!TEX root = ../Doc/equationFreeDoc.tex

%{
\subsection{\texttt{PIG()}: Projective Integration with a General macrosolver}
\label{sec:PIG}

This is an approximate Projective Integration scheme in which the macrosolver is given by any user-specified scheme. Unlike the \verb|PIRK| functions, \verb|PIG()| cannot estimate the slow vector field at the times expected by any user-specified scheme, but instead provides an estimate of the slow vector field at a slightly different time, after an application of the microsolver. Consequently \verb|PIG()| will incur an additional global error term proportional to the burst length of the microscale simulator. For that reason, \verb|PIG()| should be used with very stiff problems, in which the burst length of the microsolver can be short, or with the `constraint defined manifold' based microsolver provided by \verb|cdmc()|, that attempts to project the variables onto the slow manifold without affecting the time.

\begin{matlab}
%}
function [t, x, tms, xms, svf] = PIG(solver,bT,macro,IC)
%{
\end{matlab}

The inputs and outputs are necessarily a little different to the two \verb|PIRK2| functions.
\paragraph{Input}
\begin{itemize}
\item \verb|solver()|, a function that produces output from the user-specified code for micro-scale simulation. Usage:\\
\verb|[tout,xout] = solver(t_in,x_in,tSim)|
Inputs:
\verb|t_in|, the initialisation time; \verb|x_in| \(\in\RR^n\),
 the initial state; \verb|tSim|, the time to simulate for. \\
 Outputs:
 \verb|tout|, the vector of solution times, and \verb|xout|, the corresponding states.
\item \verb|bT|, a scalar, the minimum amount of time thought needed for integration of the microsolver to relax the fast variables to the slow manifold.
  \end{itemize}
  The remaining inputs to \verb|PIG()| set the solver and parameters used for macroscale simulation.
  \begin{itemize}
  \item \verb|macro|, a struct holding information about the macrosolver.
  \begin{itemize}
  \item \verb|macro.solver()|, the numerical method that the user wants to apply on a slow time scale. The solver should be formatted as a standard numerical method in Matlab/Octave that is called as
  \verb|[t_out,x_out] = solver(f,tspan,IC)| for an ordinary differential equation \(\frac{dx}{dt}=f(t,x)\), vector of input times \verb|tspan| and initial condition \verb|IC|. The function \verb|f(t,x)| is not an input for \verb|PIG()| but will instead be estimated by PI.
\item \verb|macro.tspan|, a vector of times at which the user requests
output, of which the first element is always the initial time. If
macro.solver can adaptively select time steps (e.g. \verb|ode45|), then
tspan can consist of an initial and final time only.
\end{itemize}

\item \verb|IC|, an \(n\)-vector of initial values at the time \verb|tspan(1)|.
\end{itemize}

\paragraph{Output}
Standard usage is to output only macrosolver information, with the following usage:\\

\verb|x = PIG(micro,macro,IC)|
\begin{itemize}
\item  \verb|x|, a cell array. \verb|x{1}|, an \(\ell\)-vector of times at which PI produced output. \verb|x{2}|, an \(\ell \times n\) array of the approximate solution vector. Each row corresponds to an element of \verb|x{1}|.
\end{itemize}
It is also possible to return the microsolver applications called by the PI method in executing the user-defined macrosolver. Much of this microscale data will not be an accurate solution of the system, but rather will consist of simulations used to relax the fast variables close to the slow manifold in the process of executing a single macroscale time step.


\verb|[x,xm] = PIRK4(micro,tspan,IC)|
\begin{itemize}
\item  \verb|xm|, a cell array containing the output of all applications of the microsolver. \verb|xm{1}| is an \(L\) dimensional column vector containing times; \verb|xm{2}| is an \(L\times n\) array of the corresponding microsolver output.  \end{itemize}
\verb|[x,xm,dx] = PIRK4(micro,tspan,IC)|
\begin{itemize}
\item  \verb|dx|, a cell array containing the PI estimates of the slow vector field. \verb|dx{1}| is an \(\ell\) dimensional column vector containing all times at which the PI scheme extrapolated along microsolver data to form a macrostep. \verb|dx{2}| is an \(\ell\times n\) array containing the estimated slow vector field.
\end{itemize}

\begin{body}
The main body of the function now follows.

Get microscale and macroscale information from inputs, and compute the number of time steps at which output is expected.
\begin{matlab}
%}
tspan = macro.tspan;
csolve = macro.solver;
nT=length(tspan)-1;
sIC = length(IC);
%{
\end{matlab}

Get the number of expected outputs and set logical indices to flag what data should be saved.
\begin{matlab}
%}
nArgs=nargout();
saveMicro = (nArgs>1);
saveSvf = (nArgs>2);
%{
\end{matlab}



Run a first application of the microsolver on the initial conditions. This is done in addition to the microsolver in the main loop, because the initial conditions are often far from the attracting slow manifold.
\begin{matlab}
%}
IC = reshape(IC,[],1);
[relax_t,relax_IC] = solver(tspan(1),IC,bT);
%{
\end{matlab}

Update the initial time.
\begin{matlab}
%}
tspan(1) = tspan(1)+bT;
%{
\end{matlab}
Allocate cell arrays for times and states for any of the outputs requested by the user. If saving information, then
 record the first application of the microsolver. Note that it is unknown a priori how many applications of the microsolver
 will be required; this code may be run more efficiently if the correct number is used in place of \verb|N+1| as the dimension
 of the cell arrays.
\begin{matlab}
%}
if saveMicro
    tms=cell(nT+1,1); xms=cell(nT+1,1);
    n=1;
    tms{n} = reshape(relax_t,[],1);
    xms{n} = relax_IC;
    
    if saveSvf
        svf.t = cell(nT+1,1);
        svf.dx = cell(nT+1,1);
    end
end
%{
\end{matlab}

The idea of \verb|PIG()| is to use the output from the microsolver to approximate an unknown function \verb|ff(t,x)|,
 that describes the slow dynamics. This approximation is then used in the user-defined `coarse solver' \verb|csolve()|.
 The approximation is described in


\begin{matlab}
%}
function [dx]=genProjection(tt,xx)
%{
\end{matlab}
Run a microsolver from the given initial conditions.
        \begin{matlab}
%}
    [t_tmp,x_micro_tmp] = solver(tt,reshape(xx,[],1),bT);
%{
\end{matlab}
Compute the standard PI approximation of the slow vector field.
\begin{matlab}
%}
	del = t_tmp(end)-t_tmp(end-1);
	dx = (x_micro_tmp(end,:)-x_micro_tmp(end-1,:))'/(del);
%{
\end{matlab}
Save the microscale data, and the PI slow vector field, if requested.
\begin{matlab}
%}
	if saveMicro
        n=n+1;
        tms{n} = [reshape(t_tmp,[],1); nan];
        xms{n} = [x_micro_tmp; nan(1,sIC)];
        if saveSvf
            svf.t{n-1} = tt;
            svf.dx{n-1} = dx;
        end
	end
%{
\end{matlab}
End \verb|genProjection()|.
\begin{matlab}
%}
end
%{
\end{matlab}


Define the approximate slow vector field according to PI.
\begin{matlab}
%}
ff=@(t,x) genProjection(t,x);
%{
\end{matlab}


Do Projective Integration of \verb|ff()| with the user-specified solver.
\begin{matlab}
%}
[t,x]=feval(csolve,ff,tspan,relax_IC(end,:)');
%{
\end{matlab}

Write over \verb|x(1,:)| and \verb|t(1)|, which the user expect to be \verb|IC| and \verb|tspan(1)| respectively, with the given initial conditions.

\begin{matlab}
%}
x(1,:) = IC';
t(1) = tspan(1);
%{
\end{matlab}

Output the macroscale steps:
\begin{matlab}
%}

%{
\end{matlab}

For each additional requested output, concatenate all the cells of time and state data into two arrays. Then, return the two arrays
 in a cell.
\begin{matlab}
%}
if saveMicro
    tms = cell2mat(tms);
    xms = cell2mat(xms);
    if saveSvf
        svf.t = cell2mat(svf.t);
        svf.dx = cell2mat(svf.dx);
    end
end
%{
\end{matlab}


This concludes \verb|PIG()|.
\begin{matlab}
%}
end
%{
\end{matlab}
\end{body}
%}
