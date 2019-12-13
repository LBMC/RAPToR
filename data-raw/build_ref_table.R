ref_table <- read.table("data-raw/ref_table.csv", h=T, sep = ',', stringsAsFactors = F)

save("ref_table", file = "data/ref_table.RData")

rm(ref_table)