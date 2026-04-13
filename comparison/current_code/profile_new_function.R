library(profvis)

p <- profvis(expr = {
	source("function_current.R")

})

htmlwidgets::saveWidget(p, "profile_report.html", selfcontained = TRUE)
