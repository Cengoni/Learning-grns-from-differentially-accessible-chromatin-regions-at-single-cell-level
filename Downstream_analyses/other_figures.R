
library(ggplot2)
library(ggrepel)

# Figure 3.5 C
x<-read.csv("/home/claudia/Desktop/UGent/Thesis/networks_second/stats_named.csv", header=F)
x
rownames(x) <- x$V3
x <- x[, -3]
x <- x[order(x$V1), ]
cor(x$V1,x$V2)
plot(x$V1, x$V2) 
cor.test(x$V1, x$V2)
 
barplot(x$V1,  # Modify x-axis labels
        las = 2,
        names.arg = rownames(x),
        cex.names = 0.7)

x_ <- ggplot(x, aes(x = V1, y = V2)) +
    geom_point() +
    geom_text_repel(aes(label = rownames(x)), size = 4) +
    theme_classic() +
    labs(x = "Number of DARs", y = "Number of cells") +  
    geom_smooth(method = "lm", se = TRUE, alpha=.15, fill='blue')
x_

#####################################################################################

# Figure 3.9 A
tricho <- c(3, 14, 31) 
cortex <- c(3, 16, 37)
endo <- c(4, 6, 21)
tot <- data.frame(tricho,cortex, endo)
barplot(as.matrix(tot),  col = c("#3274a0ff", "#3a923aff","#e1812cff"))
legend("topright",                                    # Add legend to barplot
       legend = c("Root", "Unknown", "Known"),
       fill = c("#3274a0ff", "#3a923aff","#e1812cff"))
