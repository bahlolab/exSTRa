# An example of exSTRa usage, for known STR expansion disorder loci

## ---- strexpansion_prepare
library(exSTRa)

# read ub database
str_database <- exstra_db_read("data/repeat_disorders.xlsx") # exstra_db object

str_score <- exstra_score_read (
  file = "data/repeat_scores_at_known_loci.txt",
  database = strdatabase, # alternatively, this may be the direct file path, if no import options are required (TODO)
  groups.regex = c(case = "", control = "^control") # here
)
  


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

## ---- Crazy trimming of rep_in_read plots
sample.name <- 'SCA2-1'
locus.name <- 'SCA2'
for(trimming in 0:60) {
  strrir.trim <- trim.rep_in_read_data(strrir, trimming)
  pdf(sprintf("vision_docs/images/trimming_abc_%02d.pdf", trimming))
  with(strrir.trim$data[locus.name ], 
    plot(a + .5 * c, sqrt(3/2) * c, 
      col = "white", #ifelse(sample == sample.name, "red", "black"),
      xlim = c(0, 131 - 2 * trimming),
      ylim = c(0, sqrt(3/2) * (131 - 2 * trimming)),
      main = paste("Trim", trimming)
    )
  ) 
  with(strrir.trim$data[locus.name ], 
    text(a + .5 * c, sqrt(3/2) * c, 
      col = ifelse(sample == sample.name, "red", "black"),
      labels = sample, 
      cex = 0.5
    )
  )
  dev.off()
}

#### ---- STR score ----
plot(str_score)

loc_scores <- str_score$data[locus == "SCA6"]
KS <- ks.test(loc_scores[group == "control"]$rep, loc_scores[sample == "SCA2-1"]$rep, 
  alternative = "greater", exact = NULL)

ks_tests <- data.table(rep_score_data_ks_tests(str_score))
setkey(ks_tests, p.value)
head(ks_tests)
# Lots of SCA3 locus results, but the samples are similar visually to the normals

loc_scores <- str_score$data[locus == "SCA3"]
KS <- ks.test(loc_scores[group == "control"]$rep, loc_scores[sample == "SCA6-1"]$rep, 
  alternative = "greater", exact = NULL)

plot(str_score, "SCA3", sample_col = c("SCA6-1" = "red"))



