# IICR_selection
R and python code related to the study of Boitard et al (2022)

The R codes provided in directory *scripts/* allow to compute and plot the exact **Inverse Instantaneous Coalescence Rate (IICR) for different models with variable classes of Ne along the genome**, which can be used to approximate models with **linked selection**. More precisely:
- Models where each class is assumed **panmictic and stationary** (Figures 1-2 & S1-S4 of the study) are considered in *figures_pan_statsel.R*
- Models where each class is asusmed to evolve under a **stationary n-island model** (Figure 4) are considered in *figures_statstruc_statsel.R*
- Models where each class is asusmed to evolve under a **non-stationary n-island model** (Figures 5 & S7) are considered in *figures_varstruc_statsel.R*
- Models where each class is asusmed to evolve under a **non-stationary panmictic model** (i.e. with population size changes, Figures 3, 6 & S5) are considered in *figures_pan_varsel.R*.

In addition, empirical IICRs (i.e. based on simulated coalescence times) for a single selective sweep model (Figure 6) can be obtained using the script *iicr_selection.sh*. The script *psmc_tests.py* simulates genomic sequences for models with variable Ne along the genome and prints the PSMC command that was usd to estimate the IICR from these data (Figure 7).



