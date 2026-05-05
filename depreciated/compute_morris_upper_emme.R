library(sensitivity)
library(tidyverse)

getwd()


sa <- readRDS("morris_screening/sa_object.rds")
lumped_statistics <- read.csv("compute_morris/lumped_stats_upper_emme.csv")



normalized_lumped_statistics <- lumped_statistics %>%
    group_by(location) %>%
    mutate(
      normalized_total_area = total_area/ mean(total_area),
      normalized_temperature = total_averaged_temperature / mean(total_averaged_temperature )
    ) %>% ungroup()

# Check the parameter matrix
head(sa$X)
nrow(sa$X)



# One fresh copy per location
sa_ilfis           <- sa
sa_zollbrueck      <- sa
sa_schuepbachkanal <- sa
sa_schuepbach      <- sa


 

# Tell each one its own Y
tell(sa_ilfis, normalized_lumped_statistics %>% filter(location == "ilfis") %>% arrange(morris_row_index) %>% pull(normalized_total_area))
tell(sa_zollbrueck,     normalized_lumped_statistics %>% filter(location == "zollbrueck")     %>% arrange(morris_row_index) %>% pull(normalized_total_area))
tell(sa_schuepbachkanal, normalized_lumped_statistics %>% filter(location == "schuepbachkanal") %>% arrange(morris_row_index) %>% pull(normalized_total_area))
tell(sa_schuepbach,     normalized_lumped_statistics %>% filter(location == "schuepbach")     %>% arrange(morris_row_index) %>% pull(normalized_total_area))


# Now each can be inspected independently
print("zollbrueck")
print(sa_zollbrueck)
plot(sa_zollbrueck)

# Now each can be inspected independently
print("ilfis")
print(sa_ilfis)
plot(sa_ilfis)

# Now each can be inspected independently
print("schuepbachkanal")
print(sa_schuepbachkanal)
plot(sa_schuepbachkanal)

# Now each can be inspected independently
print("schuepbach")
print(sa_schuepbach)
plot(sa_schuepbach)




