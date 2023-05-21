% find pseudospectra of heterogeneous diffusion in 2D on
% patches for curiosity.  Code for either edgy or 
% centred interpolation, of any order.  AJR, May 2023
% AJR, May 2023
%!TEX root = ../Doc/eqnFreeDevMan.tex
%{
\section{\texttt{heteroDiffPseudoSpectra2}: pseudospectrum 
of computational homogenisation of a 2D diffusion}
\label{sec:heteroDiffPseudoSpectra2}



This section extends to 2D the 1D code discussed in
\cref{sec:homoDiffEdgy1}. First set random heterogeneous
diffusivities of random period in each of the two
directions. Crudely normalise by the harmonic mean so the
decay time scale is roughly one. 
\begin{matlab}
%}
mPeriod = [3 3];
cHetr = exp(1*randn([mPeriod 2]));
cHetr = cHetr*mean(1./cHetr(:)) 
%{
\end{matlab}




\subsection{Compute Jacobian and its spectrum}
Let's explore the Jacobian dynamics for a range of orders of
interpolation, all for the same patch design and
heterogeneity.  Except here use a small ratio as we do not
plot.
\begin{matlab}
%}
edgyInt = true
nSubP = (2-edgyInt)*mPeriod+1+edgyInt
nPatch = [3 3] +2*edgyInt
ratio = 0.2
nLeadEvals=prod(nPatch)+max(nPatch);
leadingEvals=[];
%{
\end{matlab}

Evaluate eigenvalues for spectral as the base case for
polynomial interpolation of order \(2,4,\ldots\).
\begin{matlab}
%}
maxords = 0;
for ord = maxords:2:maxords
    ordInterp = ord    
%{
\end{matlab} 
Configure with same parameters, then because they are reset
by this configuration, restore coupling.
\begin{matlab}
%}
    configPatches2(@heteroDiff2,[-pi pi -pi pi],nan,nPatch ...
        ,ord,ratio,nSubP,'EdgyInt',edgyInt,'hetCoeffs',cHetr);
%{
\end{matlab}
Find which elements of the 6D array are interior micro-grid
points and hence correspond to dynamical variables.
\begin{matlab}
%}
    u0 = zeros([nSubP,1,1,nPatch]);
    u0([1 end],:,:) = nan;
    u0(:,[1 end],:) = nan;
    i = find(~isnan(u0));
%{
\end{matlab}
Construct the Jacobian of the scheme as the matrix of the
linear transformation, obtained by transforming the standard
unit vectors.
\begin{matlab}
%}
    nJac = length(i)
    Jac = nan(nJac);
    sizeJacobian = size(Jac)
    for j = 1:nJac
      u0(i) = (j==(1:nJac));
      dudt = patchSys2(0,u0);
      Jac(:,j) = dudt(i);
    end
%{
\end{matlab}
Test for symmetry, with error if we know it should be
symmetric.
\begin{matlab}
%}
    notSymmetric = norm(Jac-Jac')    
    if notSymmetric>1e-7, disp("failed symmetry"), end 
    Jac(abs(Jac)<1e-12) = 0;
%{
\end{matlab}
Find all the eigenvalues (as \verb|eigs| is unreliable).
\begin{matlab}
%}
    evals = eig(Jac);
    biggestImag = max(abs(imag(evals)));
    if biggestImag>0, biggestImag=biggestImag, end
%{
\end{matlab}
Sort eigenvalues on their real-part with most positive
first, and most negative last. Store the leading eigenvalues
in \verb|egs|, and write out when computed all orders.
The number of zero eigenvalues, \verb|nZeroEv|, gives
the number of decoupled systems in this patch configuration.
\begin{matlab}
%}
    [~,k] = sort(real(evals),'Descend');
    evals = evals(k);
    if ord==0, nZeroEv=sum(abs(evals(:))<1e-5), end
    leadingEvals=[leadingEvals evals(nZeroEv*(1:nLeadEvals))]
%{
\end{matlab}


\paragraph{Pseuso-spectrum of the Jacobian}
\begin{matlab}
%}
figure(1)
multiscalePseudoSpectra(Jac,1,1,200)
figfn=mfilename;
set(gcf,'PaperUnits','centimeters' ...
      ,'PaperPosition',[0 0 14 14] ...
      ,'renderer','Painters')
print('-depsc2',figfn)
matlab2tikz([figfn '.tex'],'showInfo',false ...
    ,'noSize',true,'parseStrings',false,'showWarnings',false ...
    ,'extraCode',['\tikzsetnextfilename{' figfn '}'] ...
    ,'extraAxisOptions','\extraAxisOptions' ...
    ,'checkForUpdates',false)
%{
\end{matlab}

End of the for-loop over orders of interpolation.
\begin{matlab}
%}
end
%{
\end{matlab}

End of the main script.


%\input{../Patch/heteroDiff2.m}

%}