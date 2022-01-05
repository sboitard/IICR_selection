## Figure 6: simulates coalescence times in 200 independent 15Mb sweep regions and computes the empirical IICR  based on these values

# simulates and stores coalescence times
python iicr_selection.py schrider_1000_02_200rep 200 15000000 0.0004 0.0004 1000 0.2

# plots the IICR
R -f iicr_selection_loop.R

