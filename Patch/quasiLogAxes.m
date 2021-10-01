% quasiLogAxes() transforms the current 2D plot to a
% quasi-log axes (via asinh).  AJR, 25 Sep 2021
%!TEX root = doc.tex
%{
\section{\texttt{quasiLogAxes()}: transforms plot to quasi-log axes}
\label{sec:quasiLogAxes}
%\localtableofcontents

This function rescales and labels the axes of the current
2D~plot.  The aim is to effectively show the complex
spectrum of multiscale systems such as the patch scheme. 
The eigenvalues are over a wide rab=nge of magnitudes, but
are signed.  So we use a nonlinear asinh transformation of
the axes, and then label the axes with reasonable ticks.

Herein \verb|x,y| denotes the original data scale, and
\verb|h,v| denotes nonlinearly transformed quantities.
\begin{matlab}
%}
function quasiLogAxes(handle,xScale,yScale)
%{
\end{matlab}
\paragraph{Input} 
\begin{itemize}
\item \verb|handle|: handle to your plot to transform, for
example, obtained by \verb|handle=plot(...)|
\item \verb|xScale| (optional, default~1): let \(x\) denote
every horizontal coordinate, then transform the plot-data
with the \(\text{asinh}()\) function so that
\begin{itemize}
\item for \(|x|\lesssim x_{\text{scale}}\) the horizontal
axis scaling is approximately linear, whereas
\item for \(|x|\gtrsim x_{\text{scale}}\) the horizontal
axis scaling is approximately signed-logarithmic.
\end{itemize}
\item \verb|yScale| (optional, default~1): corresponds to
\verb|xScale| for the vertical axis scaling.
\end{itemize}

\paragraph{Example}
If invoked with no arguments, then execute an example.
\begin{matlab}
%}
if nargin==0  
   % first generate your data
   n=99;  fast=(rand(n,1)<0.8);
   z = -rand(n,1).*(1+1e3*fast)+1i*randn(n,1).*(5+1e2*fast);
   % second plot data and transform axes
   handle = plot(real(z),imag(z),'o');
   xlabel('real-part'), ylabel('imag-part')
   quasiLogAxes(handle,1,10);
   return
end% example
%{
\end{matlab}


Default values for scaling.
\begin{matlab}
%}
if nargin<3, yScale=1; end
if nargin<2, xScale=1; end
%{
\end{matlab}

\paragraph{Output} None, just the transformed plot.

\begin{devMan}
Scale the plot data in the 2D~plot \verb|handle|. Give an
error if it appears that the plot-data has already been
transformed.
\begin{matlab}
%}
for k=1:length(handle)
    assert(~strcmp(handle(k).UserData,'quasiLogAxes'), ...
       'Replot graph, as it appears plot data is already transformed')
    handle(k).XData = asinh(handle(k).XData/xScale);
    handle(k).YData = asinh(handle(k).YData/yScale);
    handle(k).UserData = 'quasiLogAxes';
end%for
%{
\end{matlab}
Get limits of nonlinearly transformed data, and reset with
4\%~padding around all margins---crude but serviceable.
But if the range is too small, then set to plus/minus one.
\begin{matlab}
%}
axis tight; lims=axis; dl=0.04*diff(lims); 
dl(abs(dl)<4e-5) = 1;
lims = lims+[-dl(1) +dl(1) -dl(3) +dl(3)];
axis(lims)
%{
\end{matlab}

\paragraph{Horizontal scaling}
Get the order of magnitude of the horizontal data.
\begin{matlab}
%}
hmax=max(abs(lims(1:2)));
hmag=floor(log10(xScale*sinh(hmax)));
%{
\end{matlab}
The following decision affecting the number of ticks plotted
more-or-less assumes that the plot should include relatively
near-zero data, or the data crosses zero.  It assumes that
because, if not, then surely a user should use a standard
log-plot.
\begin{matlab}
%}
if hmag<3, ticks=[1;2;5];
  elseif hmag<7, ticks=[1;3];
  else ticks=[1]; 
  end
%{
\end{matlab}
Form a range of ticks, geometrically spaced, trim off the
small values that would be too dense near zero (omit those
within 6\% of \verb|hmax|).
\begin{matlab}
%}
ticks=ticks*10.^(hmag+(-7:0));
j=find(ticks>xScale*sinh(0.06*hmax));
ticks=sort([0;ticks(j);-ticks(j)]);
%{
\end{matlab}
Set the ticks in place according to the transformation.
Getting the `parent' gives the axes of the plot~(gca).
\begin{matlab}
%}
theaxes = get(handle(1),'parent');
set(theaxes,'Xtick',asinh(ticks/xScale) ...
    ,'XtickLabel',cellstr(num2str(ticks,4)) ...
    ,'XTickLabelRotation',30)
%{
\end{matlab}

\paragraph{Vertical scaling}
Do the same for the vertical axis.
\begin{matlab}
%}
vmax=max(abs(lims(3:4)));
vmag=floor(log10(yScale*sinh(vmax)));
if vmag<3, ticks=[1;2;5];
  elseif vmag<7, ticks=[1;3];
  else ticks=[1]; 
  end
ticks=ticks*10.^(vmag+(-7:0));
j=find(ticks>yScale*sinh(0.06*vmax));
ticks=sort([0;ticks(j);-ticks(j)]);
set(theaxes,'Ytick',asinh(ticks/yScale) ...
    ,'YtickLabel',cellstr(num2str(ticks,4)))
%{
\end{matlab}

Turn the grid on by default.
\begin{matlab}
%}
grid on
%{
\end{matlab}
\end{devMan}
%}
