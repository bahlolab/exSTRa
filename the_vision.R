# example of how the package could be used

## ---- strexpansion_prepare
library(strexpansion)

strdatabase <- strdb_read("/Users/tankard/Documents/Research/repeats/disease_repeats/repeat_disorders.xlsx") # class strdb
#strdatabase
# or 
# strdatabase = read.strs.ucsc("simpleRepeat.txt.gz")

strcounts <- strs_read(file = "/Users/tankard/Documents/Research/repeats/read_simulation/str_simulations/summarising_simulations/simulation_summary_09.txt", 
                       database = strdatabase, 
                       groups.regex = c(control = "normal", case = "expanded")
                       ) # class strdata

strcounts.byalignment <- strs_read(file = "/Users/tankard/Documents/Research/repeats/read_simulation/str_simulations/summarising_simulations/simulation_summary_03.txt", 
  database = strdatabase, 
  groups.regex = c(control = "normal", case = "expanded")
) # class strdata
#strcounts

shortplotlabels <- c("initial_normal_11" = "n11", "initial_expanded_03" = "e03")
set_plotnames(strcounts, shortplotlabels)
set_plotnames(strcounts.byalignment, shortplotlabels)

## ---- permutation_testing
strcounts.perm <- str_chisq_permutation_test(strcounts,
                                            cols = c("up_00", "up_01", "up_11", "up_02", "up_12", "up_22"), 
                                            keep.cols = c("up_01", "up_11", "up_02", "up_12"),
                                            #loci = strloci(data),
                                            B = 100000, 
                                            require.nozero = TRUE) # class strdataperm (maybe a subclass of strcount???)

plot(strcounts.perm) # plot all the disease p-values with confidence intervals in one plot

plot(strcounts.perm, multi = TRUE, auto.layout = TRUE, cex = 1.3) # plot each disease in a separate plot, with the layout made automatically

plot(strcounts.perm, multi = TRUE, auto.layout = TRUE, cex = 1.3, read.counts = NULL) 

## ---- loglinear_testing
# A different kind of test
strcounts.ll <- str_loglin_test(strcounts, 
                                cols = c("up_00", "up_01", "up_11", "up_02", "up_12", "up_22"), 
                                )
plot(strcounts.ll)

## ---- loglinear_testing_other exploration
strcounts.ll$df
class(strcounts.ll)



plot(strcounts.ll, multi = TRUE, auto.layout = TRUE, cex = 1.3)

## ---- compare_tests
image.dir <- "../strexpansion_testing/images/"
height <- 8
width = 13
dir.create(image.dir)

pdf(paste0(image.dir, "perm_test_QQ_repeat_expansion_loci.pdf"), width = width, height = height)
plot(strcounts.perm, multi = TRUE, auto.layout = TRUE, cex = 1.3, diseases = strloci(strcounts.perm))
dev.off()

pdf(paste0(image.dir, "loglin_test_QQ_repeat_expansion_loci.pdf"), width = width, height = height)
plot(strcounts.ll, multi = TRUE, auto.layout = TRUE, cex = 1.3, diseases = strloci(strcounts.perm))
dev.off()

pdf(paste0(image.dir, "perm_test_QQ_repeat_expansion_loci_single.pdf"))
plot(strcount.perm)
dev.off()

pdf(paste0(image.dir, "loglin_test_QQ_repeat_expansion_loci_single.pdf"))
plot(strcounts.ll)
dev.off()


# plot each test result vs the other
# check we can do what we think we can do: 
if(sum(rownames(strcounts.ll$p.value) != rownames(strcounts.perm$p.value)) != 0) {
  stop("Row names not equal")
}
if(sum(colnames(strcounts.ll$p.value) != colnames(strcounts.perm$p.value)) != 0) {
  stop("Col names not equal")
}
# all good :)
plot(-log10(as.vector(as.matrix(strcounts.ll$p.value))), -log10(as.vector(as.matrix(strcounts.perm$p.value))))
abline(a = 0, b = 1, col = "red")
# needs to be just of expanded samples, and just normals

make.comparable <- function(X, sampleregex) {
  Y <- as.matrix(X$p.value)
  Y <- Y[, grepl(sampleregex, colnames(Y))]
  -log10(as.vector(Y))
}


for(sample.type in c("expanded", "normal")) {
  par(mar = c(5, 4, 4, 2) + 0.1)
  pdf(paste0(image.dir, "comparison_", sample.type, "_samples.pdf"))
  x <- make.comparable(strcounts.ll, sample.type)
  y <- make.comparable(strcounts.perm, sample.type)
  r <- cor(x, y)
  plot(x, y, 
       main = sprintf("Comparision of tests for %s samples R = %.3f", sample.type, r), 
       xlab = "Log-linear model test -log10(p-value)", ylab = "Permutation test -log10(p-value)")
  abline(a = 0, b = 1, col = "red")
  dev.off()
}


## ---- loglinear_testing_inccase

# log-linear test including cases in control pop
strcounts.ll.inccase <- str_loglin_test(strcounts, 
                                cols = c("up_00", "up_01", "up_11", "up_02", "up_12", "up_22"), 
                                include.test.case = TRUE
)

## ---- loglinear_testing_inccase_plot_command

pdf(paste0(image.dir, "inccase_loglin_test_QQ_repeat_expansion_loci_single.pdf"))
## ---- loglinear_includingcase_plot
plot(strcounts.ll.inccase, main = "Q-Q all loci including case in nulls")
## ---- end
dev.off()

## ---- boxplot_strdata
pdf(paste0(image.dir, "boxplots_with_expected.pdf"))
for(locus in strloci(strcounts)) {
  boxplot(strcounts, locus, coverage = 50 * 129/149, read.length = 129, cases.known = TRUE)
}
dev.off()

pdf(paste0(image.dir, "boxplots_with_expected_by_alignment.pdf"))
for(locus in strloci(strcounts.byalignment)) {
  boxplot(strcounts.byalignment, locus, coverage = 50, read.length = 149, cases.known = TRUE)
}
dev.off()

pdf(paste0(image.dir, "boxplots_example_of_variation.pdf"))
x <- data.frame( group = rep(1:20, each = 20), 
  counts = rpois(20 * 20, 20))
boxplot(counts ~ group, x)
dev.off()
