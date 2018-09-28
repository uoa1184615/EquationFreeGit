%PIRK2 is an implementation of Projective Integration with the
%macrosolver given by an explicit second order Runge-Kutta scheme and a user-specified microsolver. JM,
%July 2018.
%!TEX root = ../equationFreeDoc.tex
%{
\subsection{\texttt{PIRK2()}}
\label{sec:PIRK2}

This is a Projective Integration scheme in which the macrosolver is given by a second order Runge-Kutta method.

\begin{matlab}
%}
function varargout = PIRK2(micro,tspan,IC)
%{
\end{matlab}

\paragraph{Input}
\begin{itemize}
\item \verb|micro| is a struct holding all information needed to run the microsolver:
\begin{itemize}
\item \verb|micro.bT|, a scalar, the minimum amount of time thought needed for integration of the microsolver to relax the fast variables to the slow manifold. 
\item \verb|micro.solver()|, a function that produces output from the user-specified code for micro-scale simulation. Usage:\\
\verb|[tout,xout] = solver(t_in,x_in,bT)|\\
Inputs:
\verb|t_in|, the initialisation time; \verb|x_in| \(\in\RR^n\),
 the initial state; \verb|bT|, the total time to simulate for. \\
 Outputs:
 \verb|tout|, the vector of solution times, and \verb|xout|, the corresponding states.
  \end{itemize}
  \end{itemize}
  The remaining inputs to \verb|PIRK2()| relate to the macroscale solution.
  \begin{itemize}
\item \verb|tspan| is an \(\ell\)-vector of times at which the user requests output, of which the first element is always the initial time. \verb|PIRK2()| does not use adaptive time stepping; the time steps used will be approximately the gap between elements of \verb|tspan|.

\item \verb|IC| is an \(n\)-vector of initial values at the time \verb|tspan(1)|.
\end{itemize}

\paragraph{Output}
Standard use is to output only estimates of the state at the times given in \verb|tspan|, with the following usage:\\

\verb|x = PIRK2(micro,tspan,IC)|
\begin{itemize}
\item  \verb|x|, an \(\ell \times n\) array of the approximate solution vector. Each row corresponds to an element of \verb|tspan|.
\end{itemize}
However, particularly while first implementing the scheme, output of the underlying PI machinery may be helpful. \verb|PIRK2()| can provide one, two or three additional outputs. These outputs all take the form of a cell array, where the first cell in the array is a vector of microscale times and the second cell is a vector or matrix of information obtained by \verb|PIRK2()| at those times.\\


\verb|[x,xm] = PIRK2(micro,tspan,IC)|
\begin{itemize}
\item  \verb|xm|, a cell array containing the output of the applications of the microsolver that immediately follow a PI macrostep. \verb|xm{1}| is an \(L\) dimensional column vector containing times; \verb|xm{2}| is an \(L\times n\) array of the corresponding microsolver output. The data in \verb|xm| is an accurate simulation of the state and can be used alongside \verb|x| when visualising the solution. \end{itemize}
\verb|[x,xm,xr] = PIRK2(micro,tspan,IC)|
\begin{itemize}
\item \verb|xr|, a cell array containing the `remaining' applications of the microsolver required by the PI method during the calculation of the macrostep. The applications of the microsolver in \verb|xr| do not have the same physical interpretation as those in \verb|xm|; they are required in order to estimate the slow vector field during the calculation of the Runge-Kutta increments, and will not in general resemble the true dynamics.
\end{itemize}
\verb|[x,xm,xr,dx] = PIRK2(micro,tspan,IC)|
\begin{itemize}
\item  \verb|dx|, a cell array containing the PI estimates of the slow vector field. \verb|dx{1}| is an \(\ell\) dimensional column vector containing all times at which the PI scheme extrapolated along microsolver data to form a macrostep. \verb|dx{2}| is an \(\ell\times n\) array containing the estimated slow vector field.
\end{itemize}

The main body of the function now follows.\\


Get microscale information from inputs.
\begin{matlab}
%}
solver = micro.solver;
bT = micro.bT;
%{
\end{matlab}

Compute the number of time steps and create storage for output.

\begin{matlab}
%}
N=length(tspan)-1; 
x=zeros(N+1,length(IC)); 
%{
\end{matlab}

Get the number of expected outputs and set logical indices to flag what data should be saved.
\begin{matlab}
%}
nargs=nargout(); 
save_micro = (nargs>1); 
save_full_micro = (nargs>2); 
save_dx = (nargs>3); 
%{
\end{matlab}


Run a first application of the microsolver on the initial conditions. This is done in addition to the microsolver in the main loop, because the initial conditions are often far from the attracting slow manifold.
\begin{matlab}
%}
IC = reshape(IC,[],1); 
[relax_t,relax_IC] = solver(tspan(1),IC,bT);
%{
\end{matlab}

Use the end point of the microsolver as the initial conditions.
\begin{matlab}
%}
tspan(1) = tspan(1)+bT; 
x(1,:)=relax_IC(end,:); 
%{
\end{matlab}

Allocate cell arrays for times and states for any of the outputs requested by the user. If saving information, then record the first application of the microsolver.
\begin{matlab}
%}
if save_micro 
    t_micro=cell(1,N+1); x_micro=cell(1,N+1);
    t_micro{1} = reshape(relax_t,[],1);
    x_micro{1} = relax_IC;
    
    if save_full_micro 
        t_r_micro = cell(1,N+1);
        x_r_micro = cell(1,N+1);
        
        if save_dx 
            dx_t = cell(1,N+1); 
            dx_hist = cell(1,N+1); 
        end
    end
end
%{
\end{matlab}


Main loop:
\begin{matlab}
%}
for n=1:N
    t=tspan(n);
%{
\end{matlab}

    Run the first application of the microsolver; since this application directly follows from the initial conditions, or from the latest macrostep, this microscale information is physically meaningful as a simulation of the system.\\
     Extract the size of the final time step. 
   \begin{matlab}
%}
 [t1,xm1] = solver(t,x(n,:)',bT);
    del=t1(end)-t1(end-1);
%{
\end{matlab}

      Find the needed time step to reach time \verb|tspan(n+1)| and form a first estimate \verb|dx1| of the slow vector field.
\begin{matlab}
%}
    Dt=tspan(n+1)-t-bT;  
    
    dx1 = (xm1(end,:)-xm1(end-1,:))/del; 
%{
\end{matlab}

  Project along \verb|dx1| to form an intermediate approximation of \verb|x|; run another application of the microsolver and form a second estimate of the slow vector field.  
\begin{matlab}
%}
    xint = xm1(end,:) + (Dt-bT)*dx1;
    
    [t2,xm2] = solver(t+Dt,xint',bT);
    
    del=t2(end)-t2(end-1);
    
    dx2 = (xm2(end,:)-xm2(end-1,:))/del; 
%{
\end{matlab}

    Use the weighted average of the estimates of the slow vector field to take a macrostep.
\begin{matlab}
%}
    x(n+1,:) = xm1(end,:) + Dt*(dx1+dx2)/2; 

%{
\end{matlab}

If saving trusted microscale data, populate the cell arrays for the current loop iterate with the time steps and output of the first application of the microsolver.    
\begin{matlab}
%}
    if save_micro 
        t_micro{n+1} = reshape(t1,[],1);
        x_micro{n+1} = xm1;
%{
\end{matlab}

If saving all microscale data, repeat for the remaining applications of the microsolver.         
\begin{matlab}
%}
        if save_full_micro 
            t_r_micro{n+1} = reshape(t2,[],1);
            x_r_micro{n+1} = xm2;
%{
\end{matlab}

If saving PI estimates of the slow vector field, populate the corresponding cells with times and estimates.
\begin{matlab}
%}
            if save_dx 
                dx_t{n} = [t1(end); t2(end)];
                dx_hist{n} = [dx1; dx2];
            end
        end
    end
    %{
\end{matlab}
Terminate the main loop:
\begin{matlab}
%}
end
%{
\end{matlab}


Write over \verb|x(1,:)|, which the user expects to correspond to the input \verb|tspan(1)|, with the given initial condition.

\begin{matlab}
%}
x(1,:) = IC'; 
%{
\end{matlab}

Output the macroscale steps:
\begin{matlab}
%}
varargout{1} = x;
%{
\end{matlab}

For each additional requested output, concatenate all the cells of time and state data into two arrays. Then, return the two arrays in a cell.
\begin{matlab}
%}
if save_micro
    t_micro = cat(1,t_micro{:});
    x_micro = cat(1,x_micro{:});
    varargout{2} = {t_micro x_micro}; 
    
    if save_full_micro
        t_r_micro = cat(1,t_r_micro{:});
        x_r_micro = cat(1,x_r_micro{:});
        varargout{3} = {t_r_micro x_r_micro}; 
        
        if save_dx
            dx_t = cat(1,dx_t{:});
            dx_hist = cat(1,dx_hist{:});
            varargout{4} = {dx_t dx_hist};
        end
    end
end
%{
\end{matlab}

This concludes \verb|PIRK2()|.
\begin{matlab}
%}
end
%{
\end{matlab}
%}


