library(ggplot2)
library(ggrepel)
library(colorspace)

corr.intersect <- read.table("results/intersection_correlations/darmanis_significant_correlations_panther_protein_classes.txt",
                             sep = "\t",
                             header = FALSE,
                             row.names = 1)

corr.intersect <- corr.intersect %>%
  mutate(csum = rev(cumsum(rev(V3))), 
                      pos = V3/2 + lead(csum, 1),
                      pos = if_else(is.na(pos), V3/2, pos))

ggplot(corr.intersect,
       aes(y = V3,
           x = reorder(V2, -V3),
           fill = V2)) +
  geom_col(width = 1, 
           color = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1),
        legend.position = "none") +
  ylab("Number of proteins") +
  xlab("Protein class") +
  scale_fill_manual(values = qualitative_hcl(length(rownames(corr.intersect)), palette = "pastel1")) 










+
  coord_polar(theta = "y", 
              start = 0) +
  theme_void() +
  scale_fill_manual(values = qualitative_hcl(length(rownames(corr.intersect)), palette = "pastel1")) + 
  theme(legend.position="bottom") +
  geom_label_repel(data = corr.intersect,
                   aes(y = pos, label = V2),
                   size = 4.5, 
                   nudge_x = 1, 
                   show.legend = FALSE,
                   max.overlaps = 50) +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()
