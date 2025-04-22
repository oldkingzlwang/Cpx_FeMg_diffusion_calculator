# Growth-Then-Diffusion Modeling of Fe-Mg in Clinopyroxene

Supplementary data and code to "[Cooling rate of clinopyroxene reveals the thickness and effusion volume of Chang’E-5 basaltic flow units](https://www.sciencedirect.com/science/article/pii/S0019103522004985)", which was published in *Icarus* in 2023.



## Overview

This repository contains MATLAB scripts and functions for modeling Fe–Mg interdiffusion in clinopyroxene crystals using a **growth-then-diffusion** approach. The approach uses the observed profile of the slow-diffusing element (Ca) as the initial condition for modeling Fe-Mg diffusion, which is more realistic than the traditional step-function assumption.

**This repository serves as the supplementary code and data to the article Wang et al. (2023) published in *Icarus***, and reproduces key modeling steps and figures based on the clinopyroxene profiles reported in that study.



## Included Files

- `FeMg_diffusion_figure.m`
  The **main script**. This is the script users should run directly to generate the five diffusion profile figures corresponding to Figure 3 of *Wang et al. (2023, Icarus)*. It automatically calls the necessary supporting functions and reads the input data from the spreadsheet.

- `FeMg_Cpx.m`
   A standalone function that numerically solves the radial diffusion equation for Fe–Mg in clinopyroxene using an implicit finite-difference method under a linear cooling path.

- `fit_to_Ca.m`
  A standalone function that estimates the effective diffusion duration necessary to reproduce the observed Ca profile using Fe–Mg diffusion parameters. This result is used to initialize the Fe–Mg diffusion calculation with a "growth-then-diffusion" condition.

- `CE5_cpx.xlsx` 
  
  Supplementary EPMA dataset from Wang et al. (2023). This Excel file contains five sheets (`cpx1`–`cpx5`), each corresponding to an EPMA profile:
  
  - **Columns A–M**: Original EPMA measurements
  - **Columns O–AB**: 6-oxygen normalized data
  - **Columns AD–AG**: Core and rim data for Ca profile
  - **Columns AH–AL**: Core and rim data for Mg# profile
  
- `README.md` – This file

Please see the detailed description for the usages in the respective `.m` files.



## Theoretical Background

### Traditional Diffusion Chronometry

Diffusion chronometry typically assumes a step-function initial condition based on core and rim compositions:
$$
C(r, 0) = \begin{cases} C_0 & \text{for } r \leq r_{\text{core}} \\ C_1 & \text{for } r > r_{\text{core}} \end{cases}
$$
This simplification neglects the effects of crystal growth during cooling, which can bias timescale estimates.

### Growth-Then-Diffusion Concept

Instead, this method models the Ca profile (a slow diffuser) first, assuming it largely preserves the growth history. This profile is then rescaled and used as the initial condition for Fe–Mg diffusion modeling, which is a fast diffuser and more sensitive to thermal relaxation. This strategy has been developed in Morgan & Blake (2006) and Brugman et al. (2022), which were adopted in this study.

### Diffusion Equation

The model solves the 1D radial diffusion equation (Fick’s second law):
$$
\frac{\partial C}{\partial t} = D(T) \left( \frac{\partial^2 C}{\partial r^2} + \frac{2}{r} \frac{\partial C}{\partial r} \right)
$$
with Neumann boundary conditions at the center and edge:
$$
\frac{\partial C}{\partial r} = 0 \quad \text{at } r = 0 \text{ and } r = A
$$


### Temperature Dependence

Temperature evolves linearly from T<sub>1</sub> to T<sub>2</sub> across the total time interval (i.e., **linear cooling**), and diffusion coefficient D(T) follows the Arrhenius relationship:
$$
D(T) = D_0 \cdot \exp\left(-\frac{E_a}{R \cdot T} \right)
$$
Parameters used (according to Müller et al. (2013) for Fe-Mg interdiffusion in clinopyroxene):

- D<sub>0</sub>=2.77×10<sup>5</sup> m<sup>2</sup>/s (pre-exponential factor)
- Ea=320.7 kJ/mol (activation energy)
- T is temperature in Kelvin



## Example Usage

```matlab
cpx = readtable('CE5_cpx.xlsx','Sheet','cpx1','VariableNamingRule','preserve');
T1 = 1033; T2 = 950;

[best_CR, ~, ~] = fit_to_Ca(cpx, T1, T2, 1);
plot_diffusion(cpx, T1, T2, best_CR, 0.0001, 0.0003, 0.0005);
```



## Additional Notes (important!)

1. **Figure Generation**
   Running the `FeMg_diffusion_figure.m` script will automatically generate five output figures. These correspond directly to the five sub-panels (a–e) in Figure 3 of Wang et al. (2023, *Icarus*), which illustrates the modeled Fe–Mg diffusion profiles in clinopyroxene.
2. **Alignment of Initial Condition Profiles**
   The function `rescale_and_shift_to_match` is used within `FeMg_diffusion_figure.m` to align the Ca-derived initial profile with the Fe–Mg diffusion model. However, due to limitations in the current implementation of this function, the resulting profile may not always be perfectly centered at the core–rim boundary. In some cases (e.g., in `cpx4.png` and `cpx5.png`), the initial condition may appear offset from what is shown in the article's Figure 3d and 3e. Manual adjustment may be needed to center the initial condition symmetrically before use in further diffusion modeling.
3. **Equivalent Diffusion Time Subtraction Method**
   This code uses an indirect approach—the “equivalent diffusion time subtraction method”—to calculate the Fe–Mg diffusion duration and associated cooling rate. This method subtracts the time required to form the Ca-derived initial profile from the total cooling time. While convenient, it is not equivalent to directly computing the Fe–Mg diffusion using a known initial condition, and may introduce small discrepancies. A more rigorous, direct method should be implemented in future improvements.
4. **Manual Selection of Cooling Rates**
   This version of the code does not include automated inversion or fitting of Fe–Mg profiles to observed data. Instead, it reproduces the results of Figure 3 in Wang et al. (2023) based on predefined values of cooling rate. These rates were manually selected by visual comparison during the original study, without optimization algorithms. As such, the reported values are approximate and reflect qualitative rather than quantitatively precise fits.



## References

- Brugman K., Till C. B., Bose M. 2022. Common assumptions and methods yield overestimated diffusive timescales, as exemplified in a Yellowstone post-caldera lava. *Contributions to Mineralogy and Petrology* 177(6): 63. https://doi.org/10.1007/s00410-022-01926-5.
- Morgan D. J., Blake S.. 2006. Magmatic residence times of zoned phenocrysts: introduction and application of the binary element diffusion modelling (BEDM) technique. *Contributions to Mineralogy and Petrology* 151: 58-70. https://doi.org/10.1007/s00410-005-0045-4.
- Müller T., Dohmen R., Becker H. W., et al. 2013. Fe–Mg interdiffusion rates in clinopyroxene: experimental data and implications for Fe–Mg exchange geothermometers. *Contributions to Mineralogy and Petrology* 166: 1563-1576. https://doi.org/10.1007/s00410-013-0941-y.
- Wang Z. L., Tian W., Li H. J., et al. 2023. Cooling rate of clinopyroxene reveals the thickness and effusion volume of Chang'E-5 basaltic flow units. *Icarus* 394: 115406. https://doi.org/10.1016/j.icarus.2022.115406.