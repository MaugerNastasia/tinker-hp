# Tinker-HP: High-Performance Massively Parallel Evolution of Tinker on CPUs & GPUs
-----------------------------------------------------------------------------------------------------------------------------------------------
A phase advance update (GPUs) has been pushed to GitHub (28/10/2020) to support COVID-19 research: the Tinker-HP (multi)-GPUs plateform is now available: https://github.com/TinkerTools/tinker-hp/tree/master/GPU 

!!! Check out the JCTC paper (Open Access) about the native GPUs version of Tinker-HP (03/23/2021) : 
Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems. O. Adjoua, L. Lagardère, L.-H. Jolly, Arnaud Durocher, Z. Wang, T. Very, I. Dupays, T. Jaffrelot Inizan, F. Célerse, P. Ren, J. Ponder, J-P. Piquemal, J. Chem. Theory. Comput., 2021, XX, XX, online (Open Access) https://doi.org/10.1021/acs.jctc.0c01164 
-----------------------------------------------------------------------------------------------------------------------------------------------
<H2><B>Versions</B></H2>

<B>Update 03/02/2021 : PLUMED Support for GPUs !</B>

<B>Update 02/11/2021 : PLUMED Support for version 1.2 (CPUs) !</B>

<B>Update 11/24/2020 : all versions have been pushed to GitHub</B>

Current Github version: 1.1v (enhanced AVX512 vectorized CPUs version), 1.2 (CPUs) + 1.2/phase advance (multi)-GPUs

Current Development version: 1.3 (CPUs + multi-GPUs)


All releases of the Tinker-HP code are now being performed on Github. For news, benchmarks and additional tutorials, please visit the Tinker-HP website
http://tinker-hp.ip2ct.upmc.fr/   (new website design to appear)

Tinker-HP is a CPUs and GPUs based, multi-precision, MPI massively parallel package dedicated to long polarizable molecular dynamics simulations and to polarizable QM/MM. Tinker-HP is an evolution of the popular Tinker package that conserves it simplicity of use but brings new 
capabilities allowing performing very long molecular dynamics simulations on modern supercomputers that use thousands of cores. 
The Tinker-HP approach offers various strategies using domain decomposition techniques for periodic boundary conditions in the 
framework of the (n)log(n) Smooth Particle Mesh Ewald. Tinker-HP proposes a high performance scalable computing environment for 
polarizable (AMOEBA, Amberpol...) and classical (Amber, Charmm, OPLS...) force fields giving access to large systems up to millions of atoms. It can be used on supercomputers as well as on lab clusters. Tinker-HP supports Intel (AVX5212 enhanced version) and AMD CPUs platforms as well as NVIDIA GPUs (1080, 2080, 3090, P100, V100, A100). 

Tinker-HP is available free of charge for ALL Academic Institutions, National Laboratories and supercomputer centers through the global Tinker license (https://dasher.wustl.edu/tinker/downloads/license.pdf).
Non-academic entities (e.g., companies, for profit organizations) should contact the managing universities (see license).

If you want support, support is only provided to registered users:

i) <B>Please fill in the form at:</B>
https://dasher.wustl.edu/tinker/downloads/license.pdf

and sent it back to TinkerHP_Download@ip2ct.upmc.fr to be added to the user support.

ii) <B>Please cite:</B>

- Tinker-HP: a Massively Parallel Molecular Dynamics Package for Multiscale Simulations of Large Complex Systems 
with Advanced Polarizable Force Fields.
L. Lagardère, L.-H. Jolly, F. Lipparini, F. Aviat, B. Stamm, Z. F. Jing, M. Harger, H. Torabifard, G. A. Cisneros, 
M. J. Schnieders, N. Gresh, Y. Maday, P. Ren, J. W. Ponder, J.-P. Piquemal, Chem. Sci., 2018, 9, 956-972 (Open Access)
https://doi.org/10.1039/C7SC04531J

- if you are using the AVX512 vectorized version (1.1v) dedicated to Intel's CPUs (Skylake, CascadeLake etc...), please also cite:
Raising the Performance of the Tinker-HP Molecular Modeling Package [Article v1.0].
L. H. Jolly, A. Duran, L. Lagardère, J. W. Ponder, P. Y. Ren, J.-P. Piquemal, LiveCoMS, 2019, 1 (2), 10409  (Open Access)
 https://doi.org/10.33011/livecoms.1.2.10409
 
- if you are using the GPUs version, please also cite:
Tinker-HP : Accelerating Molecular Dynamics Simulations of Large Complex Systems with Advanced Point Dipole Polarizable Force Fields using GPUs and Multi-GPUs systems
Olivier Adjoua, Louis Lagardère, Luc-Henri Jolly, Arnaud Durocher, Thibaut Very, Isabelle Dupays, Zhi Wang, Théo Jaffrelot Inizan, Frédéric Célerse, Pengyu Ren, Jay W. Ponder, Jean-Philip Piquemal, J. Chem. Theory. Comput., 2021, XX, XX (Open Access) https://doi.org/10.1021/acs.jctc.0c01164

iii) <b>Tinkertools</b>
Tinker-HP is part of the Tinker distribution and uses the same tools as Tinker. These tools can be found here : https://github.com/TinkerTools/tinker
if you use the Tinkertools please cite :

Tinker 8: Software Tools for Molecular Design.
J. A. Rackers, Z. Wang, C. Lu, M. L. Maury, L. Lagardère, M. J. Schnieders, J.-P. Piquemal, P. Ren, J. W. Ponder,  J. Chem. Theory. Comput., 2018, 14 (10), 5273–5289 DOI: http://dx.doi.org/10.1021/acs.jctc.8b00529 ; PMC free text : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6335969/

iii) <B>Support:</B>

We provide support to registered users only (see i) ).

Email: TinkerHP_Support@ip2ct.upmc.fr


