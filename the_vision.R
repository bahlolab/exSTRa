# example of how the package could be used

library(strexpansion)

strdatabase <- strdb_read("/Users/tankard/Documents/Research/repeats/disease_repeats/repeat_disorders.xlsx") # class strdb
strdatabase
# or 
# strdatabase = read.strs.ucsc("simpleRepeat.txt.gz")

strcounts <- strs_read(file = "/Users/tankard/Documents/Research/repeats/read_simulation/str_simulations/summarising_simulations/simulation_summary_03.txt", database = strdatabase) # class strdata
strcounts

strcount.perm <- str_chisq_permutation_test(strcounts) # class strdataperm (maybe a subclass of strcount???)

plot(strcount.perm) # plot all the disease p-values with confidence intervals in one plot

plot(strcount.perm, multi = TRUE, auto.layout = TRUE) # plot each disease in a separate plot, with the layout made automatically

