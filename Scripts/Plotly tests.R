# Read the CSV file with semicolon as delimiter
data <- read.csv(here::here("Data/bigdiv_snpduo_results_full2.csv"), sep = ";")

# View the first few rows of the data
head(data)

# Color of nodes: country of origin ----
# Load required library
library(igraph)

# Use KING coefficient for filtering
subset_data <- data[1:100, 108:110]
ids <- data$ID.x[1:100]
# Trim whitespace from country names
countries <- trimws(data$Country_origin.x[1:100])

# Create an empty edge list
edges <- data.frame(from = character(), to = character(), 
                    weight = numeric(), relationship = character())

# Create edges between IDs based on the KING values
for(i in 1:(length(ids)-1)) {
  for(j in (i+1):length(ids)) {
    king_value <- subset_data[i, "KING.snpduo"]
    
    if(!is.na(king_value)) {
      # Only keep first and second degree relationships
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

# Create the graph object
g <- graph_from_data_frame(edges, directed = FALSE, vertices = ids)

# Calculate degree for each vertex
node_degrees <- degree(g)

# Normalize node sizes to be between 2 and 8
node_sizes <- 2 + (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 6

# Create color palette for unique countries
unique_countries <- sort(unique(countries))  # Sort countries alphabetically
country_colors <- rainbow(length(unique_countries), alpha = 0.7)
names(country_colors) <- unique_countries

# Assign colors to nodes based on country
node_colors <- country_colors[countries]

# Define edge styles based on relationship (changed to match previous chart)
edge_colors <- ifelse(edges$relationship == "First degree", 
                      adjustcolor("red", alpha.f = 0.6), 
                      adjustcolor("blue", alpha.f = 0.6))
edge_widths <- ifelse(edges$relationship == "First degree", 1, 0.2)
edge_types <- ifelse(edges$relationship == "First degree", 1, 2)

# Create custom layout with more spread
layout <- layout_with_fr(g)
layout <- layout * 1.5

# Set up the plotting area with extra space on the right
par(mar = c(5, 4, 4, 12), xpd = TRUE)  # Increase right margin and allow plotting outside

# Create the plot
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.label.cex = 0.4,
     vertex.label = NA,
     vertex.color = node_colors,
     vertex.frame.color = "gray30",
     edge.width = edge_widths,
     edge.color = edge_colors,
     edge.lty = edge_types,
     main = "Related Pairs Network\nNode size: number of connections, Color: country of origin")

# Add two legends
# Relationship legend (modified to match previous chart)
legend(1.2, 1,  # Adjust these values to position the legend 
       title = "Relationship",
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = c(1, 0.2),
       lty = c(1, 2),
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       bty = "n",
       title.adj = 0.15)

# Country legend
legend(1.2, 0.4,  # Adjust these values to position the legend
       title = "Country of Origin",
       legend = unique_countries,
       col = country_colors,
       pch = 19,
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       bty = "n",
       title.adj = 0.15)

# Create regular intervals for size legend
legend_degrees <- seq(0, 80, by = 20)
# Use exactly the same scaling as for the nodes
legend_sizes <- 2 + (legend_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 6 * 0.25

# Add node size legend
legend(1.2, -0.6,  # Adjusted position
       title = "Number of Connections",
       legend = paste(legend_degrees, "connections"),
       pt.cex = legend_sizes,
       pch = 21,
       pt.bg = adjustcolor("lightblue", alpha.f = 0.7),
       pt.lwd = 1,
       col = "gray30",
       cex = 0.6,
       bty = "n",
       title.adj = 0.15)

# Reset the plotting parameters to default (important for subsequent plots)
par(mar = c(5, 4, 4, 2), xpd = FALSE)

# Use plotly ----
library(plotly)

# Get the layout coordinates from your igraph object
layout_coords <- layout_with_fr(g) * 1.5

# Create a data frame for nodes
nodes_df <- data.frame(
  x = layout_coords[,1],
  y = layout_coords[,2],
  id = V(g)$name,
  country = countries[match(V(g)$name, ids)],
  degree = node_degrees
)

# Create edge data frame with correct mapping
edges_df <- data.frame(
  x = layout_coords[match(edges$from, V(g)$name), 1],
  y = layout_coords[match(edges$from, V(g)$name), 2],
  xend = layout_coords[match(edges$to, V(g)$name), 1],
  yend = layout_coords[match(edges$to, V(g)$name), 2],
  relationship = edges$relationship
)

# Create the plot
p <- plot_ly() %>%
  # Add first degree edges
  add_segments(
    data = subset(edges_df, relationship == "First degree"),
    x = ~x, y = ~y, xend = ~xend, yend = ~yend,
    line = list(
      color = "rgba(255,0,0,0.6)",
      width = 1
    ),
    name = "First degree",
    showlegend = TRUE,
    legendgroup = "relationships",
    legendgrouptitle = list(text = "Relationship Type")
  ) %>%
  # Add second degree edges
  add_segments(
    data = subset(edges_df, relationship == "Second degree"),
    x = ~x, y = ~y, xend = ~xend, yend = ~yend,
    line = list(
      color = "rgba(0,0,255,0.6)",
      width = 0.2,
      dash = "dash"
    ),
    name = "Second degree",
    showlegend = TRUE,
    legendgroup = "relationships"
  )

# Add size legend entries with evenly spaced values
size_values <- seq(0, 80, by = 20)
for(size_val in size_values) {
  p <- p %>% add_trace(
    x = c(Inf), y = c(Inf),  # Place point outside visible area
    type = "scatter",
    mode = "markers",
    marker = list(
      size = size_val/4,  # Linear scaling
      color = "rgb(128, 128, 128)",  # Solid gray fill
      line = list(color = "rgb(80, 80, 80)", width = 1)  # Darker gray border
    ),
    name = paste(size_val, "connections"),
    showlegend = TRUE,
    legendgroup = "sizes",
    legendgrouptitle = list(text = "Number of Connections"),
    inherit = FALSE
  )
}

# Add first country with legend group title
first_country <- unique_countries[1]
country_nodes <- nodes_df[nodes_df$country == first_country,]
p <- p %>% add_trace(
  data = country_nodes,
  x = ~x, y = ~y,
  type = "scatter",
  mode = "markers",
  marker = list(
    size = node_sizes[match(country_nodes$id, ids)] * 5,
    color = country_colors[first_country],
    line = list(color = "gray30", width = 1),
    symbol = "circle"
  ),
  name = first_country,
  hoverinfo = "text",
  hovertext = ~paste("ID:", id, "<br>Country:", country, 
                     "<br>Connections:", degree),
  showlegend = TRUE,
  legendgroup = "countries",
  legendgrouptitle = list(text = "Country of Origin")
)

# Add remaining countries
for(country in unique_countries[-1]) {
  country_nodes <- nodes_df[nodes_df$country == country,]
  p <- p %>% add_trace(
    data = country_nodes,
    x = ~x, y = ~y,
    type = "scatter",
    mode = "markers",
    marker = list(
      size = node_sizes[match(country_nodes$id, ids)] * 5,
      color = country_colors[country],
      line = list(color = "gray30", width = 1),
      symbol = "circle"
    ),
    name = country,
    hoverinfo = "text",
    hovertext = ~paste("ID:", id, "<br>Country:", country, 
                       "<br>Connections:", degree),
    showlegend = TRUE,
    legendgroup = "countries"
  )
}

# Update layout with specific legend settings
p <- p %>% layout(
  title = "Related Pairs Network",
  xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
  yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
  hovermode = "closest",
  legend = list(
    x = 1.1,
    y = 0.9,
    tracegroupgap = 30,
    grouptitlefont = list(size = 12),
    itemsizing = list(countries = "constant")  # Apply constant sizing only to country legend group
  )
)

p
