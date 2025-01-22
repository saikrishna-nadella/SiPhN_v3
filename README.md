# SiPhN
**SiPhN** (pronounced as "siphon") is a program to analyse Single Phase Natural circulation. It can perform steady state, transient and stability analysis of natural circulation loops with a heater (surface heating and/or internal heat generation) and a cooler with any working fluid essentially in single phase. Analysis is possible in various orientations of heater and cooler.

This package is based on my Master's thesis work. All the formulation and solution techniques implemented in the SiPhN are based on my Master's thesis (included in the repository).

It is a very easy-to-use package with User Manual included. Installation instructions are also included.

Currently, two versions are available: SiPhN_v2 (v2.1.2) and SiPhN_v3 (v3.2.0).

## SiPhN_v2 (v2.1.2)

It consists of 5 modules (listed below) and an utility (Thermo-Physical Property Calculator)

1. SSA_SHF
2. SSA_IHG
3. LSA_SHF
4. LSA_IHG
5. TA

The modules 2 and 4 (having formulation to analyze natural circulation loops with internal heat generation in working fluid) are exclusively available in this version. For all other modules, it is recommended to use SiPhN_v3.

This version is currently not available publicly. If required for your research, drop an email at saikrishna.nadella9@gmail.com.

## SiPhN_v3 (v3.2.0)

It consists of 3 modules (listed below) and a utility (Thermo-Physical Property Calculator)

1. SSA_SHF
2. LSA_SHF
3. TA
   
This is the latest version and has more features for the above 3 modules compared to those modules in SiPhN_v2 (refer User manual for v3.0 for details of additional features). However, for steady state and linear stability analysis of natural circulation loops with internal heat generation in working fluid, use modules 2 and 4 of SiPhN_v2.

This package is shared here under GPL 3.0 license (see LICENSE). I hope this package helps you, if you are looking to do some research on natural circulation loops. If this package is used in your research work, please aknowledge it by following citation

"Saikrishna Nadella, Development of computer code for stability analysis of molten salt natural circulation loop with and without internal heat generation, M. Tech Thesis, December 2017, Homi Bhabha National Institute, Mumbai." 


## Publications based on this code
1. Saikrishna Nadella, Abhishek Kumar Srivastava, Naresh Kumar Maheshwari, A semi-analytical model for linear stability analysis of rectangular natural circulation loops, Chemical Engineering Science, Volume 192, 2018, Pages 892-905, ISSN 0009-2509, https://doi.org/10.1016/j.ces.2018.08.034.
2. A.K. Srivastava, N. Saikrishna, N.K. Maheshwari, Steady state performance of molten salt natural circulation loop with different orientations of heater and cooler, Applied Thermal Engineering, Volume 218, 2023, 119318, ISSN 1359-4311, https://doi.org/10.1016/j.applthermaleng.2022.119318.
3. Srivastava, A. K., Saikrishna, N., and Maheshwari, N. K. (September 15, 2023). "Linear and Nonlinear Stability Analysis of Molten Salt Natural Circulation Loop." ASME. ASME J of Nuclear Rad Sci. October 2023; 9(4): 041401. https://doi.org/10.1115/1.4063040.

## Acknowledgement:

I sincerely thank Dr. Abhishek Kumar Srivastava, PhD (https://www.researchgate.net/profile/Abhishek-Srivastava-10) for his valuable feedback for incorporating various features and improving user-friendliness of the code.
