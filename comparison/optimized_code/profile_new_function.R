library(profvis)


setwd("/home/etienne/Desktop/repos/Github_Enterprise/BSc_project/comparison/optimized_code")



p <- profvis(expr = {
  source("function_reworked.R")
})

htmlwidgets::saveWidget(p, "profile_report.html", selfcontained = TRUE)