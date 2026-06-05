This GitHub Repo contains code used for analysing the Lorentzian output of flexural resonators in order to estimate the 
density and viscosity of fluids. 

The code currently utilises three analytical models derived from literature, being
the least-squares algorithm - applied to both the real and imaginary admittance/impedance spectra,
then Reittinger's frequency tracking model [1] and Zhang's Hydrodynamic model [2].

Utilising a calibration resistor file, and calibration files for the fork in air and fluids, the code
easily uses whichever models desired to produce density and viscosity measurements to a high accuracy.
Simply by placing in the desired values into SampleCode.py, the fork's output is analysed.

Parallel capacitance is currently analysed using a method from Reittinger's patent [1], by taking a 
weighted average of the imaginary admittance spectrum. There are other methods, that haven't been implemented,
potentially the first thing I would improve with this code.

The code can definitely be improved for efficiency and accuracy, as well as adding more models to estimate 
density and viscosity, as well as to estimate the parallel capacitance. Not currently included in the code
is implementation of dielectric constant measurement, the change parallel capacitance between the resonator
in vacuum and the resonator in fluid is directly proportional to the dielectric constant [3] so these resonators
could be used as a sensor for dielectric constant as well.

Description of models used can be found in the pdf attached of my master's thesis:
[Master_s_thesis-5.pdf](https://github.com/user-attachments/files/28636348/Master_s_thesis-5.pdf)

## Contact
For questions, bug reports, or collaboration inquiries, please open a GitHub Issue or contact j.bond6@lancaster.ac.uk

## References
[1] P. W. Reittinger, System and method for determining producibility of a formation using flexural mechanical resonator measurements, 
US Patent 7,844,401, 2010
[2] M. Zhang et al. A Hydrodynamic Model for Measuring Fluid Density and Viscosity by Using Quartz Tuning Forks, 
Sensors (Basel, Switzerland), vol. 20, no. 1, p. 198, 2019
[3] L. F. Matsiev, Application of flexural mechanical resonators to high throughput liquid characterization, IEEE Ultrasonics Symposium. Proceedings. An International Symposium (Cat. No.00CH37121), IEEE, 2000
