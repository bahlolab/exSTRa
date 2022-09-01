library(exSTRa)
library(plotly)

wg2hd <- copy(exstra_wgs_pcr_2["HD"])

ggexstra_ecdf(wg2hd)
ggplotly()
