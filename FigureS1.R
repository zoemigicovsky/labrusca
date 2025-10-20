library(tidyverse)
library(cowplot)
# Load in eigenvalue table
eigenval <- scan("filtered_chr1_19_ld_pca.eigenval")
# create pve object for pca
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)

# Load in supplemental Data to plot PCs

labrusca_info <-  read_delim("TableS1.csv")

geo_colors <- c("Connecticut" = "#1b9e77", 
                "Maine" = "#d95f02",
                "Massachusetts" = "#7570b3",
                "New Hampshire" = "#e7298a",
                "New Jersey" = "#66a61e",
                "New York" = "#e6ab02",
                "North Carolina" = "#a6761d",
                "Pennsylvania" = "#666666",
                "Vermont" = "#1f78b4",
                "Virginia" = "#b2df8a",
                "Unknown" = "black")

# PCA plot (top panel)
pca_plot <- labrusca_info %>%
  ggplot(aes(PC1, PC2, colour = GEO, shape = Source)) + 
  geom_point(alpha = 0.6, size = 3) +  
  theme_light() +
  scale_colour_manual(values = geo_colors) +
  labs(x = paste0("PC1 (", signif(pve$pve[1], 3), "%)"),
       y = paste0("PC2 (", signif(pve$pve[2], 3), "%)"),
       colour = "Location",
       shape = "Source") +
  theme(
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    legend.title = element_text(size = 8, color = "black", face = "bold"),
    legend.text = element_text(size = 7, color = "black")
  )

# Boxplot of PC1 (bottom panel)
pc1_box <- labrusca_info %>%
  ggplot(aes(x = GEO, y = PC1, fill = GEO)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.6, alpha = 0.6) +
  facet_wrap(~Source, scales = "free_x") +
  scale_fill_manual(values = geo_colors) +
  theme_bw() +
  labs(x = "State", y = "PC1 Value") +
  theme(
    strip.text = element_text(face = "bold", colour="black"),
    strip.background = element_rect(fill = "white", color = "black", size=1),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    text = element_text(color = "black")
  )

pc2_box <- labrusca_info %>%
  ggplot(aes(x = GEO, y = PC2, fill = GEO)) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.7) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.6, alpha = 0.6) +
  facet_wrap(~Source, scales = "free_x") +
  scale_fill_manual(values = geo_colors) +
  theme_bw() +
  labs(x = "State", y = "PC2 Value") +
  theme(
    strip.text = element_text(face = "bold", colour="black"),
    strip.background = element_rect(fill = "white", color = "black", size=1),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    text = element_text(color = "black")
  )
# Combine using cowplot
combined_plot <- plot_grid(pca_plot, pc1_box, pc2_box, ncol = 1, rel_heights = c(2, 1, 1), labels = c("A", "B", "C"))


# Save to file
ggsave("FigureS1.pdf", combined_plot, width = 8, height = 11)

