# acoustic-tweezer-and-dynamic-stability-of-Mie-objects
A minimal reproducible package related to retrieval algorithm with stability constraints. In the codes, the partial derivatives of acoustic radiation forces are provided to evaluate the Lyapunov Stability.

Major features:
- Developed for retrieving the transducer amplitudes and phases of a transducer array for the stable trapping and manipulation of Mie objects.
- Implements a nonlinear inverse framework in which the acoustic radiation force and torque are expressed explicitly as analytical functions of the transducer parameters.
- Enables the direct retrieval of transducer amplitudes and phases once the desired acoustic radiation force and torque acting on the object are prescribed.
- Integrates the inversion algorithm with the Lyapunov stability theorem through the partial derivatives of the acoustic radiation forces.
- Ensures that the retrieved transducer parameters not only generate the prescribed acoustic radiation force and torque but also satisfy the translational dynamic equilibrium stability conditions.

Theoretical background:
- The partial wave expansion method [1];
- The translation addition theorem [2];
- The conformal mapping technique [3];
- The Lyapunov stability theorem [4].

## Initialization

The program support parallel computation (Refer to user manual in ```./doc/user_manual``` for details). 

The source codes are found in the ```./src``` folder. To start the Soundiation GUI, add the folder path:

``` matlab
addpath('<folderpath>/Soundiation-Acoustophoresis-main/src');
```

Then, open the the Soundiation GUI by typing in Command Window (MATLAB):

``` matlab
main_interface;
```

## Requirements

Soundiation has been tested with MATLAB2010a and above and should run on most personal laptops and desktop machines.

## Documentation

- File ```arXiv_2202.04526``` is a description preprint for the software.

- ```./docs``` folder contains:
  -  A numerical model (by COMSOL Multiphysics 5.5) to validate the calculation results, if needed;
  -  A user manual.

- ```./src``` folder contains:
  -  All source codes (".m") and a GUI framework ("main_interface.fig") for the software.

- ```./data file (example)``` folder contains: 
  -  An example of the user designed particle geometry ("particle_data.stl");
  -  An example of the predicted dynamic data (namely "Myfilename.txt" by default).

   **Note**: By default, the above data files are automatically saved in the ```./src``` folder.


## Functionality

Major functionalities includes:
- Design a non-spherical particle and output a "particle_data.stl" file;
- Prediction of the acoustic radiation force and torque on non-spherical particles;
- Prediction of the dynamics (translational and rotational motions) of non-spherical particles above an user-specified transducer array (the dynamic data is saved in a "Myfilename.txt" by default).


## Contact
Tianquan Tang

- Email address: tianquan@connect.hku.hk; ttqtianquan@gmail.com.

- Researchgate: https://www.researchgate.net/profile/Tianquan-Tang.


## Major references and paper

- References

[1] E. G. Williams, Fourier acoustics: sound radiation and nearfeld acoustical holography, Academic Press, 1999, Chapter 6.

[2] P. A. Martin, Multiple scattering: interaction of time-harmonic waves with N obstacles, Cambridge University Press, 2006.

[3] T. Tang, L. Huang, An efficient semi-analytical procedure to calculate acoustic radiation force and torque for axisymmetric irregular bodies, Journal of Sound and Vibration 532 (2022) 117012. DOI: https://doi.org/10.1016/j.jsv.2022.117012

[4] A. P. Seyranian and A. A. Mailybaev, Multiparameter stability theory with mechanical applications, World Scientifc, 2003.


- Paper

[1] T. Tang, M. Wu, L. Huang, Partial derivatives of acoustic radiation force and dynamic equilibrium stability of acoustophoresis, Physical Review E (2026)

