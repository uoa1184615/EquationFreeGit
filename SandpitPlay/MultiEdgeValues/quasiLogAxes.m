% quasiLogAxes() transforms selected axes of the given plot
% to a quasi-log axes (via asinh).  AJR, 25 Sep 2021 -- 17
% Apr 2023
%!TEX root = doc.tex
%{
\section{\texttt{quasiLogAxes()}: transforms some axes of a
plot to quasi-log}
\label{sec:quasiLogAxes}

This function rescales some coordinates and labels the axes
of the given 2D or 3D~plot.  The original aim was to
effectively show the complex spectrum of multiscale systems
such as the patch scheme. The eigenvalues are over a wide
range of magnitudes, but are signed.  So we use a nonlinear
asinh transformation of the axes, and then label the axes
with reasonable ticks. The nonlinear rescaling is useful in
other scenarios also.

\begin{matlab}
%}
function quasiLogAxes(handle,xScale,yScale,zScale,cScale)
%{
\end{matlab}
\paragraph{Input} 
\begin{itemize}
\item \verb|handle|: handle to your plot to transform, for
example, obtained by \verb|handle=plot(...)|

\item \verb|xScale| (optional, default~inf): if inf, then no
transformation is done in this coordinate.  Otherwise, with
\(x\) denoting every horizontal coordinate, then transforms
the plot-data with the \(\text{asinh}()\) function so that
\begin{itemize}
\item for \(|x|\lesssim x_{\text{scale}}\) the x-axis scaling 
is approximately linear, whereas
\item for \(|x|\gtrsim x_{\text{scale}}\) the x-axis scaling 
is approximately signed-logarithmic.
\end{itemize}

\item \verb|yScale| (optional, default~inf): corresponds to
\verb|xScale| for the second axis scaling.

\item \verb|zScale| (optional, default~inf): corresponds to
\verb|xScale| for a third axis scaling if it exists.

\item \verb|cScale| (optional, default~inf): corresponds to
\verb|cScale| for a colormap, and colorbar scaling if it exists.

\item axis limits (optional): if the axis limits of the plot
do not 'fit' the plot data, then we assume you have set the
axis limits, in which case your limits are used (each direction 
considered separately).
\end{itemize}

\paragraph{Output} None, just the transformed plot.


\paragraph{Example}
If invoked with no arguments, then execute an example.
\begin{matlab}
%}
if nargin==0  
   % generate some data
   n=99;  fast=(rand(n,1)<0.8);
   z = -rand(n,1).*(1+1e3*fast)+1i*randn(n,1).*(5+1e2*fast);
   % plot data and transform axes
   handle = plot(real(z),imag(z),'.');
   xlabel('real-part'), ylabel('imag-part')
   quasiLogAxes(handle,1,10);
   return
end% example
%{
\end{matlab}


Default values for scaling, \verb|inf| denotes no
transformation of that axis.
\begin{matlab}
%}
if nargin<5, cScale=inf; end
if nargin<4, zScale=inf; end
if nargin<3, yScale=inf; end
if nargin<2, xScale=inf; end
%{
\end{matlab}


\begin{devMan}
Get current limits of the plot so we can attempt to detect
if a user has set some limits that we should keep.  And also
get the pointer to the axes and to the figure of the plot. 
\begin{matlab}
%}
xlim0=xlim; ylim0=ylim; zlim0=zlim; clim0=caxis;
theAxes = get(handle(1),'parent');
theFig = get(theAxes,'parent');
%{
\end{matlab}

Find overall factors so the data is nonlinearly mapped to
order oneish---so that then pgfplots et al.\ do not think
there is an overall scaling factor on the axes.
\begin{matlab}
%}
xFac=1e-99; yFac=xFac; zFac=xFac; cFac=xFac;
for k=1:length(handle)
    if ~isinf(xScale)
    temp = asinh(handle(k).XData/xScale);
    xFac = max(xFac, max(abs(temp(:)),[],'omitnan') );
    end
    if ~isinf(yScale)
    temp = asinh(handle(k).YData/yScale);
    yFac = max(yFac, max(abs(temp(:)),[],'omitnan') );
    end
    if ~isinf(zScale)
    temp = asinh(handle(k).ZData/zScale);
    zFac = max(zFac, max(abs(temp(:)),[],'omitnan') );
    end
    if ~isinf(cScale)
    temp = asinh(handle(k).CData/cScale);
    cFac = max(cFac, max(abs(temp(:)),[],'omitnan') );
    end
end%for
xFac=9/xFac; yFac=9/yFac; zFac=9/zFac; cFac=9/cFac;
%{
\end{matlab}

Scale the plot data in the plot \verb|handle|. Give an
error if it appears that the plot-data has already been
transformed.   Color data has to be transformed first
because usually there is automatic flow from z-data to c-data.
\begin{matlab}
%}
for k=1:length(handle)
    assert(~strcmp(handle(k).UserData,'quasiLogAxes'), ...
       'Replot graph---it appears plot data is already transformed')
    if ~isinf(cScale)
    handle(k).CData = cFac*asinh(handle(k).CData/cScale); 
    end
    if ~isinf(xScale)
    handle(k).XData = xFac*asinh(handle(k).XData/xScale); 
    end
    if ~isinf(yScale)
    handle(k).YData = yFac*asinh(handle(k).YData/yScale); 
    end
    if ~isinf(zScale)
    handle(k).ZData = zFac*asinh(handle(k).ZData/zScale); 
    end
    handle(k).UserData = 'quasiLogAxes';
end%for
if ~isinf(xScale), xlim0=xFac*asinh(xlim0/xScale); end
if ~isinf(yScale), ylim0=yFac*asinh(ylim0/yScale); end
if ~isinf(zScale), zlim0=zFac*asinh(zlim0/zScale); end
if ~isinf(cScale), clim0=cFac*asinh(clim0/cScale); end
%{
\end{matlab}
Get limits of nonlinearly transformed data, and reset with
4\%~padding around all margins---crude but serviceable.
\begin{matlab}
%}
axis tight; 
xlim1=xlim+0.04*diff(xlim)*[-1 1];
ylim1=ylim+0.04*diff(ylim)*[-1 1];
zlim1=zlim+0.04*diff(zlim)*[-1 1];
clim1=caxis+ 0*diff(caxis)*[-1 1];
%{
\end{matlab}
But if the scaled range is too different from the original,
then restore the original.  Then set the scaled limits.
\begin{matlab}
%}
if diff(xlim1)<0.5*diff(xlim0) | diff(xlim1)>2*diff(xlim0)
    xlim1=xlim0; end
if diff(ylim1)<0.5*diff(ylim0) | diff(ylim1)>2*diff(ylim0)
    ylim1=ylim0; end
if diff(zlim1)<0.5*diff(zlim0) | diff(zlim1)>2*diff(zlim0)
    zlim1=zlim0; end
if diff(clim1)<0.5*diff(clim0) | diff(clim1)>2*diff(clim0)
    clim1=clim0; end
xlim(xlim1); ylim(ylim1); zlim(zlim1); caxis(clim1);
%{
\end{matlab}

\paragraph{Tick marks on the axes}
\begin{matlab}
%}
if ~isinf(xScale)
    tickingQuasiLogAxes(theAxes,'X',xlim1,xScale,xFac)
end%if
if ~isinf(yScale)
    tickingQuasiLogAxes(theAxes,'Y',ylim1,yScale,yFac)
end%if
if ~isinf(zScale)
    tickingQuasiLogAxes(theAxes,'Z',zlim1,zScale,zFac)
end%if
%{
\end{matlab}
But for color, only if we can find a colorbar. 
\begin{matlab}
%}
if ~isinf(cScale)
  for p=1:numel(theFig.Children)
    ca = theFig.Children(p);
    if class(ca) == "matlab.graphics.illustration.ColorBar"
      tickingQuasiLogAxes(ca,'C',clim1,cScale,cFac)
      break  
    end
  end
end%if
%{
\end{matlab}
Turn the grid on by default.
\begin{matlab}
%}
grid on
end%function
%{
\end{matlab}


\subsection{\texttt{tickingQuasiLogAxes()}: typeset ticks
and labels on an axis}
\begin{matlab}
%}
function tickingQuasiLogAxes(ca,Q,qlim1,qScale,qFac)
%{
\end{matlab}
\paragraph{Input} 
\begin{itemize}
\item \verb|ca|: pointer to axes/colorbar dataset.
\item \verb|Q|: character, either \texttt{X,Y,Z,C}.
\item \verb|qlim1|: the scaled limits of the axis.
\item \verb|qScale|: the scaling parameter for the axis.
\item \verb|qFac|: the scaling factor for the axis.
\end{itemize}

\paragraph{Output} None, just the ticked and labelled axes.

Get the order of magnitude of the horizontal data.
\begin{matlab}
%}
    qmax=max(abs(qlim1));
    qmag=floor(log10(qScale*sinh(qmax/qFac)));
%{
\end{matlab}
Form a range of ticks, geometrically spaced, trim off the
small values that would be too dense near zero (omit those
within 6\% of \verb|qmax|).
\begin{matlab}
%}
    ticks=10.^(qmag+(-7:0));
    j=find(ticks>qScale*sinh(0.06*qmax/qFac));
    nj=length(j);
    if nj<3,     ticks=[1;2;5]*ticks(j);
    elseif nj<5, ticks=[1;3]*ticks(j);
    else         ticks=ticks(j);
    end
    ticks=sort([0;ticks(:);-ticks(:)]);
%{
\end{matlab}
Set the ticks in place according to the transformation.
\begin{matlab}
%}
    if Q=='C', p='s'; Q=''; else p=''; end
    set(ca,[Q 'Tick' p],qFac*asinh(ticks/qScale) ...
          ,[Q 'TickLabel' p],cellstr(num2str(ticks,4)))
    if Q=='X', set(ca,[Q 'TickLabelRotation'],40), end
end%function qScaling
%{
\end{matlab}
\end{devMan}
%}
