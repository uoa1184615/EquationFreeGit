
<h3>Equation-Free function toolbox for Matlab/Octave</h3>

<h4>Users</h4> 
<a href="https://github.com/uoa1184615/EquationFreeGit">
Click to download the toolbox via Github</a>, and then place
the toolbox's folder in a path searched by Matlab/Octave. 

<p>This <em>equation-free toolbox</em> empowers the
computer-assisted analysis of complex, multiscale systems.
Its aim is to enable you to use microscopic simulators to
perform system level tasks and analysis, because microscale
simulations are often the best available description of a
system. The methodology bypasses the derivation of
macroscopic evolution equations by computing only short
bursts of of the microscale simulator, and often only
computing on small patches of the spatial domain. This suite
of functions empowers users to start implementing such
methods in their own applications. The document
<a href="https://github.com/uoa1184615/EquationFreeGit/raw/master/eqnFreeUserMan-newest.pdf">
<tt>eqnFreeUserMan-newest.pdf</tt></a> describes how to use the
main functions of interest to you.

<p><a href="https://au.mathworks.com/matlabcentral/fileexchange/73632-equation-free-toolbox"> 
View Equation-Free Toolbox on Matlab's File Exchange</a>

<p><strong>Please contact us</strong> via
<a href="http://www.maths.adelaide.edu.au/anthony.roberts/">
http://www.maths.adelaide.edu.au/anthony.roberts/ </a>



<p><img src="https://github.com/uoa1184615/EquationFreeGit/raw/master/Patch/Figs/nonlinearDiffusion.png" style="max-width:100%; width=300"/>

<h4>Sparse simulation on only small patches of space</h4>

The above graph illustrates an `equation-free' computation on
only small well-separated patches of the spatial domain. The
micro-scale simulations within each patch, here a nonlinear
diffusive system, are craftily coupled to neighbouring
patches and thus interact to provide accurate macro-scale
predictions over the whole spatial domain. We have proved
that the patches may be tiny, and still the scheme makes
accurate macro-scale predictions. Thus the computational
savings may be enormous, especially when combined with
projective integration



<h4>Projective Integration skips over time</h4>

<img src="https://github.com/uoa1184615/EquationFreeGit/raw/master/ProjInt/Figs/egPIMM2.png" width=250
align="left" style="margin: 0px 10px 5px 0px;"/>
Simulation&nbsp;over&nbsp;time is a complementary dynamic problem. The
`equation-free' approach is to simulate for only short
bursts of time, and then to extrapolate over un-simulated
time into the future, or into the past, or perform system
level analysis. This graph plots one example where the gaps
in time show the un-computed times between bursts of
computation.




<h4>Contributors</h4> 
This project aims to collectively develop a
Matlab/Octave toolbox of <em>equation-free</em> algorithms.
Initially the algorithms are to be straightforward.  The
plan is to subsequently develop more and more capability.
<a href="https://github.com/uoa1184615/EquationFreeGit/raw/master/Doc/eqnFreeDevMan.pdf">
<tt>Doc/eqnFreeDevMan.pdf</tt></a> fully details the
function code and many examples.  The Appendix outlines how
to contribute to the project.

<p>Paul Petersik is also developing some equation-free
projective integration, but in python, see
<a href="https://github.com/pjpetersik/eqnfree">
https://github.com/pjpetersik/eqnfree </a>

