# Atlung Method for Intercalant Diffusion --- AMID
API for analyzing diffusion measurement data using the AMID.

## Installation Notes
It is recommended to install Python via anaconda: [download here](https://www.anaconda.com/products/distribution)
Conda recommends creating a virtual environment for install dependencies: [venv doc](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
The environment.yml file supplied in this repository can be used.
If you have a system Python install and prefer to use pip and venv, the requirements.txt file in this 
repository can be used instead.

## How to use the template Jupyter Notebook
Press shift+enter to run each block of code

### Specifying source and destination file paths
The file paths for folder structures will be different on various operating systems (Windows, Mac, Linux). If the file path is not reading your data correctly, check if the slashes representing folders should be forward (/) or back (\) slashes. 

### Input parameters for fitting
Particle radius: this should be the primary particle radius for polycrystalline active materials, in units of cm. This can be measured from SEM images with ImageJ, see publications below for detailed examples. For single crystal materials, the particle radius can be measured by PSA or SEM.

ftol: the tolerance criteria for the fits. If fits are poor, try making smaller. If a RuntimeWarning is raised, try making larger. You may still get a good fit despite a RuntimeWarning being raised. Defaults to 5e-14.

D_bounds: bounds for Dc. Needs to be a list with the lower and upper bound. You can search the literature for your material's expected diffusivity to find a range. Generally for layered NMC materials, we would recommend the bounds for Dc be [1e-16, 1e-9]. If the outputted Dc calculated by this analysis program is at the limit of one of these ranges, then you can re-run the analysis with a new expanded range.

shape: particle geometry. Currently supports "sphere" and "plane". Defaults to "sphere".

nalpha: the number of expansion coefficients to use in the Atlung expression. Default to 150.

nQ : the number of Q values (on a log scale) to use in the theorical Atlung curve, τ vs Q , for comparing fit quality. Default is 2000. There is normally no reason to change this.

save: whether or not to save figures comparing theory to fitted values for each voltage interval. Default is True.

label: an additional label that be used for saving figures and data. For example, perhaps one wishes to use different particle sizes, then "label=r1" and "label=r2" could be used in 2 separate calls to fit_atlung().

## Sample protocols

Protocol for NMC/Li half cells in a voltage window of 3.0-4.3 V:

Formation: charge and discharge current density of C/20 (10 mA/g, 8 µA/cm2) for the first cycle and then charged to 4.3 V at C/40. 

Discharge direction AMID: After formation, charge cells to 4.3V at C/40, then discharge in 0.1 V intervals from 4.3V until 3.6 V at current densities of 2C, 1C, C/2.5, C/5, C/10, C/20, C/40, C/80, and C/160 with a 15 min open circuit relaxation between each current step. Then a final discharge interval from 3.6-3.0V with the same C-rates from 2C to C/160.

Charge directon AMID: After formation, discharge cells to 3.0V at C/40, then charge in AMID voltage intervals of 3.0-3.7V, and then subsequent 0.1V intervals to 4.3V. In each interval the same C-rates should be used as in the discharge protocol (2C to C/160) with 15 mins OCV periods in between each charge step.


## Sample data

Example protocol voltage vs time and voltage vs capacity plots for an NMC622/Li, low-IR designed diffusion cell at 30oC:
![protocol_vis_EZ_C-263_01V_30C_220705](https://user-images.githubusercontent.com/95938840/184234331-6ab90da3-53b9-40f4-9957-2a71876fb896.jpg)
Example capacity vs rate plot:
![cap-rate_EZ_C-263_01V_30C_220705](https://user-images.githubusercontent.com/95938840/184234346-fae86e91-b525-428d-9fe0-58e7edd3056d.jpg)
Example Atlung sphere fittings for various voltage intervals:
![EZ_C-263_01V_30C_220705_Atlung-sphere_4 250](https://user-images.githubusercontent.com/95938840/184234360-67aae4d7-2904-4cb8-954a-8f0f35be8cb4.jpg)
![EZ_C-263_01V_30C_220705_Atlung-sphere_4 150](https://user-images.githubusercontent.com/95938840/184234677-479bfe95-9837-4bb5-bec3-c2a071c50ddd.jpg)
![EZ_C-263_01V_30C_220705_Atlung-sphere_4 050](https://user-images.githubusercontent.com/95938840/184234683-e7859ab3-265b-4311-bf3f-1ad5b4051373.jpg)
![EZ_C-263_01V_30C_220705_Atlung-sphere_3 950](https://user-images.githubusercontent.com/95938840/184234687-45f27d5d-6bff-481d-aaf1-d256d06ef65e.jpg)
![EZ_C-263_01V_30C_220705_Atlung-sphere_3 850](https://user-images.githubusercontent.com/95938840/184234695-2d271a93-41b3-46b2-9552-e6cd667a1177.jpg)
![EZ_C-263_01V_30C_220705_Atlung-sphere_3 750](https://user-images.githubusercontent.com/95938840/184234698-3fddfa12-e45b-4d6b-8894-32edf5ecafe2.jpg)
![EZ_C-263_01V_30C_220705_Atlung-sphere_3 650](https://user-images.githubusercontent.com/95938840/184234707-63517af1-dea5-411a-9ee4-51b29f1f1725.jpg)
![EZ_C-263_01V_30C_220705_Atlung-sphere_3 300](https://user-images.githubusercontent.com/95938840/184234712-24a5f28f-d4e2-4948-b266-cad154df5ea3.jpg)

Example figures of merit Diffusivity vs Voltage output plots. For a good quality of fit the following metrics are desirable: low fit error, high capacity span, low IR drop.
![D-V_EZ_C-263_01V_30C_220705_](https://user-images.githubusercontent.com/95938840/184234373-9b55173c-d49a-4c2a-8be4-a066d32d5050.jpg)


# Citation
When using this analysis tool in your work, please cite and give credit to at least one of the following papers: 

Liu, A.; Phattharasupakun, N.; Cormier, M. M. E.; Zsoldos, E.; Zhang, N.; Lyle, E.; Arab, P.; Sawangphruk, M.; Dahn, J. R. Factors That Affect Capacity in the Low Voltage Kinetic Hindrance Region of Ni-Rich Positive Electrode Materials and Diffusion Measurements from a Reinvented Approach. J. Electrochem. Soc. 2021, 168 (7), 070503. https://doi.org/10.1149/1945-7111/ac0d69.

Phattharasupakun, N.; Cormier, M. M. E.; Lyle, E.; Zsoldos, E.; Liu, A.; Geng, C.; Liu, Y.; Li, H.; Sawangphruk, M.; Dahn, J. R. Correlating Cation Mixing with Li Kinetics: Electrochemical and Li Diffusion Measurements on Li-Deficient LiNiO 2 and Li-Excess LiNi 0.5 Mn 0.5 O 2. J. Electrochem. Soc. 2021, 168 (9), 090535. https://doi.org/10.1149/1945-7111/ac24ba.	

Phattharasupakun, N.; Cormier, M. M. E.; Liu, Y.; Geng, C.; Zsoldos, E.; Hamam, I.; Liu, A.; Johnson, M. B.; Sawangphruk, M.; Dahn, J. R. A Baseline Kinetic Study of Co-Free Layered Li 1+x (Ni 0.5 Mn 0.5 ) 1−x O 2 Positive Electrode Materials for Lithium-Ion Batteries. J. Electrochem. Soc. 2021, 168 (11), 110502. https://doi.org/10.1149/1945-7111/ac3157.


# Troubleshooting Help
For any questions, you can email us at: mball@dal.ca, eniko.zsoldos@dal.ca, marc.cormier@dal.ca
