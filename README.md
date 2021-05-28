# IICR_selection
R and python code related to the study of Boitard et al (under review)

The R codes provided in directory *scripts/* allow to compute and plot the exact **Inverse Instantaneous Coalescence Rate (IICR) for different models with variable classes of Ne along the genome**, which can be used to approximate models with **linked selection**. More precisely:
- Models where each class is assumed **panmictic and stationary** (Figures 1-3 of the study) are considered in *figures_pan_statsel.R*
- Models where each class is asusmed to evolve under a **stationary n-island model** (Figure 4) are considered in *figures_statstruc_statsel.R*
- Models where each class is asusmed to evolve under a **non-stationary n-island model** (Figure 5) are considered in *figures_varstruc_statsel.R*
- Models where each class is asusmed to evolve under a **non-stationary panmictic model** (i.e. with population size changes, Figure 6) are considered in *figures_pan_varsel.R*.

In addition, empirical IICRs (i.e. based on simulated coalescence times) for a single selective sweep model (Figure 6) can be obtained using the script *iicr_selection.sh*. 

