# example of how the package could be used

library(strexpansion)

strdatabase = read.strs.xlsx("repeat_disorders.xlsx") # class strdb
# or 
# strdatabase = read.strs.ucsc("simpleRepeat.txt.gz")

strcounts <- read.strs(data = "read_counts.txt", database = strdatabase) # class strdata

strcount.perm <- str_permutation_testing(strcounts) # class strdataperm (maybe a subclass of strcount???)

plot(strcount.perm) # plot all the disease p-values with confidence intervals in one plot

plot(strcount.perm, multi = TRUE, auto.layout = TRUE) # plot each disease in a separate plot, with the layout made automatically

