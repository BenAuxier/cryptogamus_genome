# redo_fig4_sabine_May2020. Peter Bourke, Plant breeding WUR
######################################################################################

#Preliminary work: I commented out, included here FYI

# setwd("Sabine_manuscript_Figure4/") #where I have all the files you sent by email 15.05.2020
# source("compareLG_polymapr.R")

## First I want to put the map data together into a more convenient structure, like a list
# dim(finalfinalfinal_map) #1417 4
# finalfinalfinal_map[1:5,]
# unique(finalfinalfinal_map$LG)


## Rename this map1, and the _man forced order map is map2
# map1 <- setNames(lapply(unique(finalfinalfinal_map$LG), function(lg)
#   matrix(unlist(finalfinalfinal_map[finalfinalfinal_map$LG == lg,2:3]),
#          ncol=2,
#          dimnames = list(as.character(finalfinalfinal_map[finalfinalfinal_map$LG == lg,1]),
#                          colnames(finalfinalfinal_map)[2:3]))),paste0("LG",unique(finalfinalfinal_map$LG)))

#Ik heb de names van Table 2 gehaalt, alleen in plaats van LG5_LG11, gewoon LG5_11 om het koorter te maken. Kan hier makkelijk zelf veranderen

# map2 <- list("LG1" = matrix(unlist(LG1[,2:3]), 
#                             ncol=2,
#                             dimnames = list(as.character(LG1[,1]),colnames(LG1)[2:3])),
#              "LG2" = matrix(unlist(LG2[,2:3]), 
#                             ncol=2,
#                             dimnames = list(as.character(LG2[,1]),colnames(LG2)[2:3])),
#              "LG3b" = matrix(unlist(LG3b[,2:3]), 
#                              ncol=2,
#                              dimnames = list(as.character(LG3b[,1]),colnames(LG3b)[2:3])),
#              "LG4" = matrix(unlist(LG4[,2:3]), 
#                             ncol=2,
#                             dimnames = list(as.character(LG4[,1]),colnames(LG4)[2:3])),
#              "LG5_11" = matrix(unlist(LG5[,2:3]), 
#                                ncol=2,
#                                dimnames = list(as.character(LG5[,1]),colnames(LG5)[2:3])),
#              "LG6" = matrix(unlist(LG6[,2:3]), 
#                             ncol=2,
#                             dimnames = list(as.character(LG6[,1]),colnames(LG6)[2:3])),
#              "LG7_8" = matrix(unlist(LG7[,2:3]), 
#                               ncol=2,
#                               dimnames = list(as.character(LG7[,1]),colnames(LG7)[2:3])),
#              "LG9" = matrix(unlist(LG9[,2:3]), 
#                             ncol=2,
#                             dimnames = list(as.character(LG9[,1]),colnames(LG9)[2:3])),
#              "LG10_TIG058" = matrix(unlist(LG10[,2:3]), 
#                                     ncol=2,
#                                     dimnames = list(as.character(LG10[,1]),colnames(LG10)[2:3])),
#              "LG12" = matrix(unlist(LG11[,2:3]), 
#                              ncol=2,
#                              dimnames = list(as.character(LG11[,1]),colnames(LG11)[2:3])),#the numbering here is a bit confusing!
#              "LG13" = matrix(unlist(LG12[,2:3]), 
#                              ncol=2,
#                              dimnames = list(as.character(LG12[,1]),colnames(LG12)[2:3])),
#              "LG14" = matrix(unlist(LG14[,2:3]), 
#                              ncol=2,
#                              dimnames = list(as.character(LG14[,1]),colnames(LG14)[2:3])),
#              "LG15" = matrix(unlist(LG15[,2:3]), 
#                              ncol=2,
#                              dimnames = list(as.character(LG15[,1]),colnames(LG15)[2:3])))
# 
# length(map1)
# length(map2)


## Update versions of the compare_maps function to make a nice plot
## Switched to lists of matrices to make things easier later..

## Version 1 - ordinary comparison
# test_ls1 <- list("m1" = map1$LG1,
#                  "m2" = map2$LG1)

compare_maps_1 <- function(maplist,
                            chm.wd = 0.2,
                            bg.col = "white",
                            links.col = "grey42",
                            thin.links = NULL,
                           yposn = NULL,
                           displace = NULL,
                           label.cex = 1,
                            ...){
  if(!is.null(displace)){
    if(length(displace) != length(maplist)) stop("displace must be a vector of equal length to maplist")
  }
  
  ##As links are drawn between neighbouring maps only, generate a list of common markers between neighbouring maps
  common_marks <- lapply(1:(length(maplist) - 1), function(n) 
    intersect(rownames(maplist[[n]]),rownames(maplist[[n + 1]])))
  
  if(any(sapply(common_marks,length) == 0)) stop(paste0("Input had the following number of linking markers: ",
                                                        paste0(sapply(common_marks,length), collapse = ", "),
                                                        ".\nWithout linking markers, a comparison is impossible"))
  
  maxy <- max(sapply(maplist, function(x) max(x[,"position"])))
  
  ## Generalise bg.col into a vector to handle vector input on bg.col:
  if(length(bg.col) != 1){
    if(length(bg.col) != length(maplist)) stop("If supplying multiple background colours, please specify only the required number (same number of elements as maplist)!")
  } else{
    bg.col <- rep(bg.col, length(maplist))
  }
  
  ## Set up the plot area:
  plot(NULL,xlim = c(0,2*length(maplist)), 
       ylim = c(1.1*maxy,-0.1*maxy),
       xlab = "", ylab = "cM", axes = FALSE, ...)
  
  axis(2)
  
  if(!is.null(displace)){
    for(i in 1:length(displace)){
      maplist[[i]][,"position"] <- maplist[[i]][,"position"] + displace[i] 
    }
  }
  

  if(is.null(yposn)){
    yposn <- rep(par("usr")[4], length(maplist))
  } else{
    yposn <- rep(yposn, length(maplist))
  } 
  
  if(!is.null(displace)){
    yposn <- yposn + displace
  }
  
  ## Draw chromosome outlines:
  for(i in seq(length(maplist))){
    symbols(x = 2*i - 1, y = min(maplist[[i]][,"position"]), circles= chm.wd, bg=bg.col[i], add=TRUE, inches = FALSE)
    symbols(x = 2*i - 1, y = max(maplist[[i]][,"position"]), circles= chm.wd, bg=bg.col[i], add=TRUE, inches = FALSE)
    rect(2*i - 1 - chm.wd, min(maplist[[i]][,"position"]), 2*i - 1 + chm.wd, max(maplist[[i]][,"position"]), col = bg.col[i])
    
    ## Add the rest of the marker positions:
    for(j in 1:nrow(maplist[[i]])){
      segments(x0 = 2*i - 1 - chm.wd,y0 = maplist[[i]][,"position"][j],
               x1 = 2*i - 1 + chm.wd,y1 = maplist[[i]][,"position"][j])
    }
    

    
    if(!is.null(names(maplist)[i])) {
      par("xpd" = NA)
      text(x = 2*i - 1, y = yposn[i],
           labels = names(maplist)[i], font = 2,
           cex = label.cex)
      par("xpd" = FALSE)
    }

  }
  
  ## Add links between markers
  for(i in seq(length(common_marks))){
    common.mrks <- common_marks[[i]]
    
    for(mark in common.mrks){
      segments(x0 = 2*i - 1 + chm.wd,y0 = maplist[[i]][mark,"position"],
               x1 = 2*i + 1 - chm.wd,y1 = maplist[[i+1]][mark,"position"],
               col = links.col, lty = 3)
    }
  }
  
  
} 

#test:
# compare_maps_1(test_ls1,yposn = -12) #can specify the y position of the labels now


#######################################
## Next - single linkage (no comparison) for LG 14 and 15..

compare_maps_single <- function(maplist,
                           chm.wd = 0.2,
                           bg.col = "white",
                           links.col = "grey42",
                           thin.links = NULL,
                           yposn = NULL,
                           label.cex = 1,
                           ...){

  maxy <- max(maplist[[1]][,"position"])
  
  ## Set up the plot area:
  plot(NULL,xlim = c(0,4*length(maplist)), 
       ylim = c(1.1*maxy,-0.1*maxy),
       xlab = "", ylab = "cM", axes = FALSE, ...)
  
  axis(2)
  
  if(is.null(yposn)){
    yposn <- rep(par("usr")[4], length(maplist))
  } else{
    yposn <- rep(yposn, length(maplist))
  } 
  
  ## Draw chromosome outlines:
  for(i in seq(length(maplist))){
    symbols(x = 3*i - 1, y = min(maplist[[i]][,"position"]), circles= chm.wd, bg=bg.col[i], add=TRUE, inches = FALSE)
    symbols(x = 3*i - 1, y = max(maplist[[i]][,"position"]), circles= chm.wd, bg=bg.col[i], add=TRUE, inches = FALSE)
    rect(3*i - 1 - chm.wd, min(maplist[[i]][,"position"]), 3*i - 1 + chm.wd, max(maplist[[i]][,"position"]), col = bg.col[i])
    
    ## Add the rest of the marker positions:
    for(j in 1:nrow(maplist[[i]])){
      segments(x0 = 3*i - 1 - chm.wd,y0 = maplist[[i]][,"position"][j],
               x1 = 3*i - 1 + chm.wd,y1 = maplist[[i]][,"position"][j])
    }
    
    
    
    if(!is.null(names(maplist)[i])) {
      par("xpd" = NA)
      text(x = 3*i - 1, y = yposn[i],
           labels = names(maplist)[i], font = 2,
           cex = label.cex)
      par("xpd" = FALSE)
    }
    
  }
  
} 

#test
# compare_maps_single(list("LG14" = map2$LG14),
#                     bg.col = rgb(0.5,0.95,0.5,1))


## Version- stacked comparison, like for LG5 and LG11
# Note: this is not a fully general function, only suitable for your data!! Makes it faster for me to re-write..
compare_maps_stack <- function(maplist,
                           chm.wd = 0.2,
                           bg.col = "white",
                           links.col = "grey42",
                           thin.links = NULL,
                           lab_space = 15,
                           displace = c(0,20), #this can be used to tweak position of stacked maps
                           label.cex = 1,
                           maxy = NULL,
                           ...){
  
  #input has to have 2 elements:
  if(length(maplist) != 2) stop("Only 2 elements allowed")
  if(length(maplist[[1]]) != 2) stop("First element must also be a list")
  
  ##As links are drawn between neighbouring maps only, generate a list of common markers between neighbouring maps
  common_marks <- lapply(1:2, function(n) 
    intersect(rownames(maplist[[1]][[n]]),rownames(maplist[[2]])))
  
  if(any(sapply(common_marks,length) == 0)) stop(paste0("Input had the following number of linking markers: ",
                                                        paste0(sapply(common_marks,length), collapse = ", "),
                                                        ".\nWithout linking markers, a comparison is impossible"))
  
  if(is.null(maxy)) maxy <- max(maplist[[2]][,"position"])
  
  #May need to invert the smaller maps:
  if(cor(maplist[[1]][[1]][common_marks[[1]],"position"],
         maplist[[2]][common_marks[[1]],"position"]) < 0){
    print("inverting map!")
    maplist[[1]][[1]][,"position"] <- max(maplist[[1]][[1]][,"position"]) - maplist[[1]][[1]][,"position"]
  }
    
  if(cor(maplist[[1]][[2]][common_marks[[2]],"position"],
         maplist[[2]][common_marks[[2]],"position"]) < 0){
    print("inverting map!")
    maplist[[1]][[2]][,"position"] <- max(maplist[[1]][[2]][,"position"]) -  maplist[[1]][[2]][,"position"] 
  }
  

  ## Generalise bg.col into a vector to handle vector input on bg.col:
  if(length(bg.col) != 1){
    if(length(bg.col) != length(maplist)) stop("If supplying multiple background colours, please specify only the required number (same number of elements as maplist)!")
  } else{
    bg.col <- rep(bg.col, length(maplist))
  }
  
  ## Set up the plot area:
  plot(NULL,xlim = c(0,2*length(maplist)), 
       ylim = c(1.1*maxy,-0.1*maxy),
       xlab = "", ylab = "cM", axes = FALSE, ...)
  
  axis(2)
  
  ## Generate y positions:
  
  #At the very least, displace by max of sub-map1, use displace to add a suitable buffer
  maplist[[1]][[2]][,"position"] <- maplist[[1]][[2]][,"position"] + max(maplist[[1]][[1]][,"position"]) + displace[2]
  
  maplist[[1]][[1]][,"position"] <- maplist[[1]][[1]][,"position"] + displace[1]
  

  ## Draw chromosome outlines:
  symbols(x = 1, y = min(maplist[[1]][[1]][,"position"]), circles= chm.wd, bg=bg.col[1], add=TRUE, inches = FALSE)
  symbols(x = 1, y = max(maplist[[1]][[1]][,"position"]), circles= chm.wd, bg=bg.col[1], add=TRUE, inches = FALSE)
  
  symbols(x = 1, y = min(maplist[[1]][[2]][,"position"]), circles= chm.wd, bg=bg.col[1], add=TRUE, inches = FALSE)
  symbols(x = 1, y = max(maplist[[1]][[2]][,"position"]), circles= chm.wd, bg=bg.col[1], add=TRUE, inches = FALSE)
  
  symbols(x = 2*2 - 1, y = min(maplist[[2]][,"position"]), circles= chm.wd, bg=bg.col[2], add=TRUE, inches = FALSE)
  symbols(x = 2*2 - 1, y = max(maplist[[2]][,"position"]), circles= chm.wd, bg=bg.col[2], add=TRUE, inches = FALSE)
  
    
  rect(1 - chm.wd, min(maplist[[1]][[1]][,"position"]), 1 + chm.wd, max(maplist[[1]][[1]][,"position"]), col = bg.col[1])
  rect(1 - chm.wd, min(maplist[[1]][[2]][,"position"]), 1 + chm.wd, max(maplist[[1]][[2]][,"position"]), col = bg.col[1])
  rect(2*2 - 1 - chm.wd, min(maplist[[2]][,"position"]), 2*2 - 1 + chm.wd, max(maplist[[2]][,"position"]), col = bg.col[2])
  
    
  ## Add the rest of the marker positions:
  for(j in 1:nrow(maplist[[1]][[1]])){
      segments(x0 = 1 - chm.wd,y0 = maplist[[1]][[1]][,"position"][j],
               x1 = 1 + chm.wd,y1 = maplist[[1]][[1]][,"position"][j])
  }
  
  for(j in 1:nrow(maplist[[1]][[2]])){
    segments(x0 = 1 - chm.wd,y0 = maplist[[1]][[2]][,"position"][j],
             x1 = 1 + chm.wd,y1 = maplist[[1]][[2]][,"position"][j])
  }
  
  for(j in 1:nrow(maplist[[2]])){
    segments(x0 = 2*2 - 1 - chm.wd,y0 = maplist[[2]][,"position"][j],
             x1 = 2*2 - 1 + chm.wd,y1 = maplist[[2]][,"position"][j])
  }
    
  
  ## Add chromosome names:

  par("xpd" = NA)
  
  text(x = 1, y = min(maplist[[1]][[1]][,"position"]) + lab_space,
      labels = names(maplist[[1]][1]), font = 2,
      cex = label.cex)
  text(x = 1, y = min(maplist[[1]][[2]][,"position"]) + lab_space,
       labels = names(maplist[[1]][2]), font = 2,
       cex = label.cex)
  
  text(x = 3, y = min(maplist[[2]][,"position"]) + lab_space,
       labels = names(maplist[2]), font = 2,
       cex = label.cex)
      
  par("xpd" = FALSE)
    
    
  ## Add links between markers
  for(i in seq(length(common_marks))){
    common.mrks <- common_marks[[i]]
    
    for(mark in common.mrks){
      segments(x0 = 1 + chm.wd,y0 = maplist[[1]][[i]][mark,"position"],
               x1 = 3 - chm.wd,y1 = maplist[[2]][mark,"position"],
               col = links.col, lty = 3)
    }
  }
  
  
} 

#test:
# compare_maps_stack(maplist = list("stack1" = list("LG5" = map1$LG5,
#                                                   "LG11" = map1$LG11),
#                                   "LG5_11" = map2$LG5_11),
#                    lab_space = -12, displace = c(0,25),
#                    bg.col = c("white",rgb(0.5,0.95,0.5,1)),
#                    label.cex = 0.9
#                    )

##########################################################################################
## -- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- 
##########################################################################################
## Put it all together into a publication-quality plot (Figure 4)

# save(map1, map2, file = "Fig4data.Rda", compress = "xz")
load("Fig4data.Rda")
################
## PDF version
################
#set up background colours (can change)
background_cols <- c("white",rgb(0.5,0.95,0.5,1))
labsize = 0.8 #select size of labels, default 1

pdf("Figure4.pdf",width = 5,height = 6)
layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9,
                10:13),ncol=4,byrow=T),
       widths = c(1,1,0.5,0.5))

par(mar=c(0.5,4,3.5,0.5))
compare_maps_1(list("LG1" = map1$LG1,
                    "LG1" = map2$LG1),
               bg.col = background_cols,
               yposn = -12,
               label.cex = labsize)

compare_maps_1(list("LG2" = map1$LG2,
                    "LG2" = map2$LG2),
               bg.col = background_cols,
               yposn = -22,
               label.cex = labsize)

compare_maps_1(list("LG3" = map1$LG3,
                    "LG3b" = map2$LG3),
               bg.col = background_cols,
               yposn = -18,
               displace = c(20,0),
               label.cex = labsize) 

compare_maps_1(list("LG4" = map1$LG4,
                    "LG4" = map2$LG4),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)


compare_maps_stack(maplist = list("stack1" = list("LG5" = map1$LG5,
                                                  "LG11" = map1$LG11),
                                  "LG5_11" = map2$LG5_11),
                   lab_space = -16, displace = c(0,35),
                   bg.col = background_cols,
                   label.cex = labsize,
                   maxy = 115
)

compare_maps_1(list("LG6" = map1$LG6,
                    "LG6" = map2$LG6),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)

par(mar=c(0,4,3.5,0.5))
compare_maps_stack(maplist = list("stack1" = list("LG7" = map1$LG7,
                                                  "LG8" = map1$LG8),
                                  "LG7_8" = map2$LG7_8),
                   lab_space = -22, displace = c(0,50),
                   bg.col = background_cols,
                   label.cex = labsize,
                   maxy = 175)

par(mar=c(0.5,4,3.5,0.5))
compare_maps_1(list("LG9" = map1$LG9,
                    "LG9" = map2$LG9),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)

compare_maps_1(list("LG10" = map1$LG10,
                    "LG10_TIG058" = map2$LG10_TIG058),
               bg.col = background_cols,
               yposn = -15,
               displace = c(70,0),
               label.cex = labsize)

compare_maps_1(list("LG12" = map1$LG12,
                    "LG12" = map2$LG12),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)

compare_maps_1(list("LG13" = map1$LG13,
                    "LG13" = map2$LG13),
               bg.col = background_cols,
               yposn = -20,
               displace = c(15,0),
               label.cex = labsize)

# par(mar=c(0.5,3,3.5,0.5))
compare_maps_single(list("LG14" = map2$LG14),
                    bg.col = background_cols[2],chm.wd = 0.9,
                    label.cex = labsize)

compare_maps_single(list("LG15" = map2$LG15),
                    bg.col = background_cols[2],chm.wd = 0.9,
                    label.cex = labsize)

dev.off()

##########################################################################################
## -- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- -- --- 
##########################################################################################



################
## PNG version
################
#set up background colours (can change)
background_cols <- c("white",rgb(0.5,0.95,0.5,1))
labsize = 0.8 #select size of labels, default 1

png("Figure4.png",width = 2500,height = 3000, res = 500)
layout(matrix(c(1,2,3,3,4,5,6,6,7,8,9,9,
                10:13),ncol=4,byrow=T),
       widths = c(1,1,0.5,0.5))

par(mar=c(0.5,4,3.5,0.5))
compare_maps_1(list("LG1" = map1$LG1,
                    "LG1" = map2$LG1),
               bg.col = background_cols,
               yposn = -12,
               label.cex = labsize)

compare_maps_1(list("LG2" = map1$LG2,
                    "LG2" = map2$LG2),
               bg.col = background_cols,
               yposn = -22,
               label.cex = labsize)

compare_maps_1(list("LG3" = map1$LG3,
                    "LG3b" = map2$LG3),
               bg.col = background_cols,
               yposn = -18,
               displace = c(20,0),
               label.cex = labsize) 

compare_maps_1(list("LG4" = map1$LG4,
                    "LG4" = map2$LG4),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)


compare_maps_stack(maplist = list("stack1" = list("LG5" = map1$LG5,
                                                  "LG11" = map1$LG11),
                                  "LG5_11" = map2$LG5_11),
                   lab_space = -16, displace = c(0,35),
                   bg.col = background_cols,
                   label.cex = labsize,
                   maxy = 115
)

compare_maps_1(list("LG6" = map1$LG6,
                    "LG6" = map2$LG6),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)

par(mar=c(0,4,3.5,0.5))
compare_maps_stack(maplist = list("stack1" = list("LG7" = map1$LG7,
                                                  "LG8" = map1$LG8),
                                  "LG7_8" = map2$LG7_8),
                   lab_space = -22, displace = c(0,50),
                   bg.col = background_cols,
                   label.cex = labsize,
                   maxy = 175)

par(mar=c(0.5,4,3.5,0.5))
compare_maps_1(list("LG9" = map1$LG9,
                    "LG9" = map2$LG9),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)

compare_maps_1(list("LG10" = map1$LG10,
                    "LG10_TIG058" = map2$LG10_TIG058),
               bg.col = background_cols,
               yposn = -15,
               displace = c(70,0),
               label.cex = labsize)

compare_maps_1(list("LG12" = map1$LG12,
                    "LG12" = map2$LG12),
               bg.col = background_cols,
               yposn = -20,
               label.cex = labsize)

compare_maps_1(list("LG13" = map1$LG13,
                    "LG13" = map2$LG13),
               bg.col = background_cols,
               yposn = -20,
               displace = c(15,0),
               label.cex = labsize)

# par(mar=c(0.5,3,3.5,0.5))
compare_maps_single(list("LG14" = map2$LG14),
                    bg.col = background_cols[2],chm.wd = 0.9,
                    label.cex = labsize)

compare_maps_single(list("LG15" = map2$LG15),
                    bg.col = background_cols[2],chm.wd = 0.9,
                    label.cex = labsize)

dev.off()
