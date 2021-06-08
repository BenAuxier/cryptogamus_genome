library(polymapR)
setwd("Sabine_manuscript_Figure4/")

#superstrict map
finalfinalfinal_map <- read.csv("outmap_1417markers.csv")

par(mfrow=c(5,3))

# LG1
LG1 <- read.csv("LG1_man_2_inv.csv", header = TRUE)
compare_maps(maplist=list("LG1" = finalfinalfinal_map[finalfinalfinal_map$LG == 1,],
                          "LG1_man" = LG1))

# LG2
LG2 <- read.csv("LG2_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG2" = finalfinalfinal_map[finalfinalfinal_map$LG == 2,],
                          "LG2_man" = LG2))

## LG3
LG3b <- read.csv("LG3b_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG3" = finalfinalfinal_map[finalfinalfinal_map$LG == 3,],
                          "LG3b_man" = LG3b))

# LG4
LG4 <- read.csv("LG4_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG4" = finalfinalfinal_map[finalfinalfinal_map$LG == 4,],
                          "LG4_man" = LG4))

# LG5
LG5 <- read.csv("LG5_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG5" = finalfinalfinal_map[finalfinalfinal_map$LG == 5,],
                          "LG5_LG11_man" = LG5,
                          "LG11" = finalfinalfinal_map[finalfinalfinal_map$LG == 11,]))

# LG6
LG6 <- read.csv("LG6_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG6" = finalfinalfinal_map[finalfinalfinal_map$LG == 6,],
                          "LG6_man" = LG6))

# LG7
LG7 <- read.csv("LG7_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG7" = finalfinalfinal_map[finalfinalfinal_map$LG == 7,],
                          "LG7_LG8man" = LG7,
                          "LG8" = finalfinalfinal_map[finalfinalfinal_map$LG == 8,]))

# LG9
LG9 <- read.csv("LG8_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG9" = finalfinalfinal_map[finalfinalfinal_map$LG == 9,],
                          "LG9_man" = LG9))

# LG10
LG10 <- read.csv("LG9_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG10" = finalfinalfinal_map[finalfinalfinal_map$LG == 10,],
                          "LG10_TIG058_man" = LG10))

# LG11
LG11 <- read.csv("LG10_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG12" = finalfinalfinal_map[finalfinalfinal_map$LG == 12,],
                          "LG12_man" = LG11))

# LG12
LG12 <- read.csv("LG11_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG13" = finalfinalfinal_map[finalfinalfinal_map$LG == 13,],
                          "LG13_man" = LG12))

# LG14
LG14 <- read.csv("LG12_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG14_man" = LG14,
                          "LG14_man" = LG14))
# LG15
LG15 <- read.csv("LG13_man_2.csv", header = TRUE)
compare_maps(maplist=list("LG15_man" = LG15,
                          "LG15_man" = LG15))



#############################################################################

LG3plot
par(mfrow=c(2,6))

#orde polymapr vs manual
#compare_maps(maplist=list("LG3_final" = finalfinalfinal_map[finalfinalfinal_map$LG == 3,],
#                          "LG3_man" = LG3map_man))

#LG3_inversion <- read.csv("//WURNET.NL/Homes/vreeb004/My Documents/1PhD/Analysis/GBS/filter_04022020/PerLG/Manually_calculated_RF/final/LG3_man_2_inversion.csv")

#orde polymapr vs manual
#compare_maps(maplist=list( "LG3_man" = LG3map_man,
#                          "LG3_final" = finalfinalfinal_map[finalfinalfinal_map$LG == 3,],
#                         "LG3_man_inv" = LG3_inversion))
#LG3map_man <- read.csv("//WURNET.NL/Homes/vreeb004/My Documents/1PhD/Writing/1_lastchapter/data&scripts/Manually2239/LG3_man_2.csv", header = TRUE)

"LG1_final" = finalfinalfinal_map[finalfinalfinal_map$LG == 1,]
LG1_final
LG1

#############################################################################

finalfinalfinal_map
