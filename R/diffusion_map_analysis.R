####################################################################  
#
# Diffusion map of all epithelial cells using the first 2 diffusion
# components.
#
# Dot represents single cell and arrows indicate 5 branches 
# starting from Stage NOR to the other stages.
#
#################################################################### 



library(dplyr)
library(ggplot2)


###### Load the data ######
dm <- readRDS("Epi_dm.rds")  # The cell embeddings of first two dimensions and the stage infomation needed 


###### Calculate the start point
df <- dm %>% group_by(stage) %>% summarise_all(mean) %>% as.data.frame()
stage <- as.character(unique(dm$stage))  # including NOR, INF, HYP, DYS and CIS

w0 <- subset(dm, stage == "NOR", select = c("DM1", "DM2"))
wot <- subset(df, stage != "NOR", select = c("DM1", "DM2"))
len <- apply(w0, 1, function(x) {
    a1 = x[1]
    b1 = x[2]
    tmp <- apply(wot, 1, function(y) {
        a2 = y[1]
        b2 = y[2]
        l <- sqrt((a1 - a2)^2 + (b1 - b2)^2)
        return(l)
    })
    return(tmp)
})
sum <- colSums(len)
ind <- dm[names(sum[sum == max(sum)]), ]


###### Visualization ######
df <- subset(df, stage != "NOR")
df <- df[5:1, ]
df$dm12 <- df$DM1 * 2 - ind[1]
df$dm22 <- df$DM2 * 2 - ind[2]

col <- c("purple", "#204426", "#1d3159", "#286aad", "#92bbe2", "#b02727")
thecol <- c("#AA22FF", "#387742", "#2E4D8C", "#3489E0", "#A5D3FF", 
    "#E33232")
names(col) <- stage
names(thecol) <- stage

ggplot(dm, aes(DM1, DM2)) + geom_point(aes(fill = stage), size = 0.1, 
    stroke = NA, shape = 21, alpha = 0.7) + geom_segment(data = df, 
    aes(x = ind[1], y = ind[2], xend = dm12, yend = dm22, color = stage), 
    arrow = arrow(length = unit(1, "cm")), size = 5) + scale_color_manual(values = col) + 
    scale_fill_manual(values = thecol)
