
The Fukushima Inverse Problem
=============================

This code and all associated files are the supplementary material to the paper
M. Martinez-Camara, I. Dokmani\'{c}, J. Ranieri, R. Scheibler, M. Vetterli, and A. Stohl,
__The Fukushima inverse problem__, ICASSP 2013, Vancouver, 2013.

The original paper is available [online](http://infoscience.epfl.ch/record/182697/files/icassp.pdf?version=2)
and attacks the problem of estimating the magnitude and timing of radioactive release
at the Fukushima Dai-Ichi power plant during the accident in March 2011.

Abstract
--------

Knowing what amount of radioactive material was released from Fukushima in
March 2011 is crucial to understand the scope of the consequences. Moreover, it
could be used in forward simulations to obtain accurate maps of deposition. But
these data are often not publicly available, or are of questionable quality. We
propose to estimate the emission waveforms by solving an inverse problem.
Previous approaches rely on a detailed expert guess of how the releases
appeared, and they produce a solution strongly biased by this guess.  If we
plant a nonexistent peak in the guess, the solution also exhibits a nonexistent
peak. We propose a method based on sparse regularization that solves the
Fukushima inverse problem blindly. Together with the atmospheric dispersion
models and worldwide radioactivity measurements our method correctly
reconstructs the times of major events during the accident, and gives plausible
estimates of the released quantities of Xenon.

Contact
-------

Do not hesitate to contact us if you have any questions or would
like some help running this code!

Marta Martinez-Camara
EPFL-IC-LCAV
BC Building, Station 14
1015 Lausanne, Switzerland

marta.martinez-camara@epfl.ch


Reproduce
---------

### Requirements

* MATLAB (We used version 2012b)
* CVX    (Available from: http://cvxr.com/cvx/)


### Settings

* If necessary, the location of the CVX installation directory
  can be changed directly in makeFigures.m

* All the parameters of the algorithms and simulations are described
  in a detailed fashion in the comments in `makeFigures.m`.
  If you need more detail, do not hesitate to contact us.


### Re-create figures from the paper

To recreate all four figures from the paper, follow these steps.

1. Open MATLAB
2. Navigate to directory folder containing this README.txt
3. Run `makeFigures`
4. Wait

Note: We could unfortunately not release the exact position of the
CTBTO stations. As a consequence, figure 2C is omitted.


### Software files

* matrixCleaning.m: 
  Reduce the transport matrix until a given condition number is obtained.

* reconstructFromSyntheticData_Xe.m: 
  Simulate reconstruction error for synthetic data for different reconstruction
  algorithm.

* reconstructSourceL1Pos.m: 
  Reconstruct the emission source using L1 minimization with positivity
  constraints.

* opt_routines: 
  Directory containing all the optimization routines we use.


### Data files

* Data/matrixGFSXe.mat   : the transport matrix

* Data/measXe.mat        : the measurement vector

* Data/aPrioriSource.mat : the a priori source used in [1]


Errata and caveats
-------------------

* The scale of the y-axis of Figure 1. as published in the ICASSP proceedings
  is wrong. The version of the paper in this repository contains the correct
  figure. This code also produces the correct figure.

* The exact lambda and epsilon parameters used to produce figure 1. were
  unfortunately not recorded when creating the figure for the paper.
  This code contains parameters producing an equally illustrative figure.

* In figure 4. of the paper published in the ICASSP proceedings, the emissions
  for 0m-50m and 300m-1000m were mistakenly shifted to the right by 5 ticks.
  Their magnitude was however correct.
  The version of the paper in this repository contains the correct figure. This
  code also produces the correct figure.

* The seeds of the random number generator experiments were unfortunately
  not recorded when running the synthetic experiments.


Reference
---------

[1] A. Stohl, P. Seibert, G. Wotawa, D. Arnold, J. F. Burkhart, S. Eckhardt, C.
    Tapia, A. Vargas, and T. J. Yasunari, “Xenon- 133 and caesium-137 releases into
    the atmosphere from the Fukushima Dai-ichi nuclear power plant: determination
    of the source term, atmospheric dispersion, and deposition,” Atmos. Chem.
    Phys., vol. 12, no. 5, pp. 2313–2343, 2012.


License
-------

2013 (c) M. Martinez-Camara, I. Dokmani\'{c}, J. Ranieri, R. Scheibler, M. Vetterli, and A. Stohl,
All the code is published under a CC-BY-SA 3.0 License
For details about the license, refer to http://creativecommons.org/licenses/by-sa/3.0/
  * For attribution of non-commercial reuse of this work, a similar notice to this one is sufficient
  * For attribution of commercial reuse of this work, please contact us.


