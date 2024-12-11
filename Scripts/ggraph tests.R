
# Read the CSV file with semicolon as delimiter
data <- read.csv(here::here("Data/bigdiv_snpduo_results_full2.csv"), sep = ";")

# View the first few rows of the data
head(data)

# Color of nodes: country of origin ----

# Load required libraries
library(ggraph)
library(tidygraph)
library(igraph)
library(dplyr)

# Use KING coefficient for filtering
subset_data <- data[1:100, 108:110]
ids <- data$ID.x[1:100]
countries <- trimws(data$Country_origin.x[1:100])

# Create edges dataframe
edges <- data.frame(from = character(), to = character(), 
                    weight = numeric(), relationship = character())

for(i in 1:(length(ids)-1)) {
  for(j in (i+1):length(ids)) {
    king_value <- subset_data[i, "KING.snpduo"]
    
    if(!is.na(king_value)) {
      if(king_value >= 0.0884) {
        relationship <- if(king_value >= 0.177) "First degree" else "Second degree"
        
        edges <- rbind(edges, data.frame(
          from = ids[i],
          to = ids[j],
          weight = king_value,
          relationship = relationship
        ))
      }
    }
  }
}

# Create nodes dataframe
nodes <- data.frame(
  id = ids,
  country = countries
)

# Create tidygraph object
graph_tidy <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

# Create a rainbow color palette for countries
n_countries <- length(unique(countries))
country_colors <- setNames(
  rainbow(n_countries, alpha = 0.7),
  sort(unique(countries))
)

# Previous code remains the same until the ggraph_plot creation

# Create a modified version for plotly conversion
ggraph_plotly <- ggraph(graph_tidy, layout = "fr") +
  geom_edge_link(aes(
    color = relationship,
    linetype = relationship,
    width = relationship,
    text = paste("Relationship:", relationship)  # Add text for hover
  ), alpha = 0.6) +
  geom_node_point(aes(
    color = country,
    size = centrality_degree(),
    text = paste("ID:", name,  # Add text for hover
                 "\nCountry:", country,
                 "\nConnections:", centrality_degree())
  )) +
  scale_edge_color_manual(
    values = c("First degree" = "red", "Second degree" = "blue")
  ) +
  scale_edge_linetype_manual(
    values = c("First degree" = "solid", "Second degree" = "dashed")
  ) +
  scale_edge_width_manual(
    values = c("First degree" = 1, "Second degree" = 0.2)
  ) +
  scale_color_manual(
    values = country_colors
  ) +
  scale_size_continuous(
    range = c(2, 8)
  ) +
  theme_minimal() +  # Use minimal theme instead of theme_graph()
  theme(
    legend.position = "right",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Related Pairs Network",
    edge_color = "Relationship type",
    edge_linetype = "Relationship type",
    edge_width = "Relationship type",
    color = "Country of origin",
    size = "Number of connections"
  )

# Convert to plotly with specific configuration (Doesn't work) ----  
interactive_plot <- ggplotly(ggraph_plotly, tooltip = "text") %>%
  layout(
    showlegend = TRUE,
    hovermode = "closest"
  )

# Display the interactive plot
interactive_plot