library(tidyverse)
library(grid)

# Load in supplemental Data to plot 
labrusca_info <-  read_delim("TableS1.csv")

# put Cluster data into long format for plotting
fast_long_data <- labrusca_info %>%
  pivot_longer(cols = starts_with("Cluster"), 
               names_to = "Cluster", 
               values_to = "AncestryProportion")

# Define color palettes
cluster_colors <- c("Cluster1" = "#E69F00",
                    "Cluster2" = "#56B4E9",
                    "Cluster3" = "#009E73")

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

all_colors <- c(cluster_colors, geo_colors)

# Custom legend function for GEO
geo_legend_custom <- function() {
  legend_data <- data.frame(State = factor(names(geo_colors), levels = sort(names(geo_colors))),
                            y = 1)
  
  ggplot(legend_data, aes(x = State, y = y, fill = State)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = geo_colors, name = "State") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 12),  # Smaller text size
      legend.key.size = unit(0.3, "cm"),     # Smaller colored squares
      legend.key.height = unit(0.5, "cm"),
      legend.box = "horizontal"
    ) +
    guides(fill = guide_legend(nrow = 1, byrow = TRUE))
}

# Cluster legend (to be placed on the right)
cluster_legend <- get_legend(
  ggplot(data.frame(Cluster = factor(c("Cluster1", "Cluster2", "Cluster3")), y = c(1, 1, 1)), aes(x = Cluster, y = y, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_fill_manual(values = cluster_colors, name = "Cluster") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.key.size = unit(0.5, "cm")
    )
)

# Plotting function for structure bars
make_faststructure_plot <- function(data, title, show_legend = TRUE) {
  ggplot(data, aes(x = fct_reorder(ind, GEO), y = AncestryProportion)) + 
    geom_bar(aes(fill = Cluster), stat = "identity", position = "stack", width = 1) +
    geom_tile(aes(x = fct_reorder(ind, GEO), y = -0.05, fill = GEO), 
              width = 1, height = 0.05, inherit.aes = FALSE) +
    scale_fill_manual(
      values = all_colors,
      breaks = names(cluster_colors),
      guide = if (show_legend) "legend" else "none"
    ) +
    theme_bw() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Black outline
      axis.text.x = element_blank(),  
      axis.text.y = element_text(size = 12),
      axis.ticks.x = element_blank(),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 12),
      panel.grid = element_blank()) +
    labs(x = "Individual (sorted by state)", y = "Ancestry proportion", title = title)
}

# Generate the custom GEO legend
geo_legend_grob <- geo_legend_custom()

# Create main plots
germplasm_plot <- fast_long_data %>%
  filter(Source == "Conserved") %>%
  make_faststructure_plot("A. Conserved", show_legend = FALSE)

wild_plot <- fast_long_data %>%
  filter(Source == "Wild") %>%
  make_faststructure_plot("B. Wild", show_legend = TRUE)

# Combine plots
combined_plots <- plot_grid(
  germplasm_plot,
  wild_plot,
  ncol = 2,
  rel_widths = c(1, 3),
  align = "h",
  axis = "tb"
)

# Create the final plot with GEO legend underneath
final_plot <- plot_grid(
  combined_plots,
  geo_legend_grob,  # Add the geo legend here
  ncol = 1,
  rel_heights = c(1, 0.1)  # Adjust relative heights to control space taken by the legend
)

# Save the final plot to a file
ggsave("Figure2.pdf",
       final_plot,
       width = 18, height = 9)
