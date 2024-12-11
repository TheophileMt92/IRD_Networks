# Hello IRD 

# Read the CSV file with semicolon as delimiter
data <- read.csv(here::here("Data/bigdiv_snpduo_results_full2.csv"), sep = ";")

# View the first few rows of the data
head(data)

# First igraph test ----
# Load required library
library(igraph)

# Use only first 30 rows
subset_data <- data[1:30, 108:110]
ids <- data$ID.x[1:30]

# Create an empty edge list
edges <- data.frame(from = character(), to = character(), weight = numeric())

# Create edges between IDs based on the SNP duo values
for(i in 1:(length(ids)-1)) {
  for(j in (i+1):length(ids)) {
    # Calculate average of the three measures for each pair
    weight <- mean(c(
      subset_data[i, "R0.snpduo"],
      subset_data[i, "R1.snpduo"],
      subset_data[i, "KING.snpduo"]
    ))
    
    # Only add edge if weight is not NA
    if(!is.na(weight)) {
      edges <- rbind(edges, data.frame(
        from = ids[i],
        to = ids[j],
        weight = abs(weight)
      ))
    }
  }
}

# Normalize weights to be between 1 and 3 for better visualization
edges$weight <- 1 + (edges$weight - min(edges$weight)) / (max(edges$weight) - min(edges$weight)) * 2

# Create the graph object
g <- graph_from_data_frame(edges, directed = FALSE, vertices = ids)

# Create the plot with some basic styling
plot(g,
     layout = layout_with_fr(g),
     vertex.size = 8,
     vertex.label.cex = 0.6,
     vertex.color = "lightblue",
     vertex.frame.color = "lightblue",
     edge.width = E(g)$weight,
     edge.color = "darkgrey",
     main = "SNP Duo Network (First 30 Rows)",
     bg = "white")

# Filter out samples: keep only first / second-degree relationships ----
# Load required library
library(igraph)

# Use KING coefficient for filtering
# First degree: KING >= 0.177
# Second degree: KING >= 0.0884
subset_data <- data[1:100, 108:110]
ids <- data$ID.x[1:100]

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

# Define colors for relationship types
edge_colors <- ifelse(edges$relationship == "First degree", "red", "blue")
# Define edge widths based on relationship
edge_widths <- ifelse(edges$relationship == "First degree", 1, 0.2)

# Create the plot
plot(g,
     layout = layout_with_fr(g),
     vertex.size = 3,
     vertex.label.cex = 0.5,
     vertex.label = NA,  # Added this line to remove labels
     vertex.color = "lightblue",
     vertex.frame.color = "darkblue",
     edge.width = edge_widths,  # Now using different widths for different relationships
     edge.color = edge_colors,
     main = "Related Pairs Network (First/Second Degree Only; First 100 rows)")

# Add legend with modified parameters
legend("bottomright", 
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = c(1, 0.2),
       cex = 0.7,        # Reduces the text size
       bty = "n",        # Removes the box around the legend
       seg.len = 2,      # Reduces the length of the lines in the legend
       x.intersp = 0.5,  # Reduces horizontal spacing
       y.intersp = 0.8)  # Reduces vertical spacing
# Size of the nodes: adjust based on number of clonal accessions in the node ----
# Load required library
library(igraph)

# Use KING coefficient for filtering
subset_data <- data[1:100, 108:110]
ids <- data$ID.x[1:100]

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

# Node size scaling function
scale_size <- function(x) {
  2 + (x - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 6
}

# Normalize node sizes to be between 2 and 8
node_sizes <- scale_size(node_degrees)

# Define colors for relationship types with transparency
edge_colors <- ifelse(edges$relationship == "First degree", 
                      adjustcolor("red", alpha.f = 0.6), 
                      adjustcolor("blue", alpha.f = 0.6))

# Define edge widths based on relationship
edge_widths <- ifelse(edges$relationship == "First degree", 1, 0.2)

# Create custom layout with more spread
layout <- layout_with_fr(g)
layout <- layout * 1.5

# Create the plot with extra space on the right for legends
par(mar = c(5, 4, 4, 8))

# Create the plot
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.label = NA,
     vertex.color = adjustcolor("lightblue", alpha.f = 0.7),
     vertex.frame.color = "darkblue",
     edge.width = edge_widths,
     edge.color = edge_colors,
     main = "Related Pairs Network")

# Add relationship legend
legend("bottomright", 
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = c(1, 0.2),
       cex = 0.7,
       bty = "n",
       seg.len = 2,
       x.intersp = 0.5,
       y.intersp = 0.8,
       title = "Relationship type",
       title.adj = 0.15)

# Create regular intervals for legend (0 to 80 by 20)
legend_degrees <- seq(0, 80, by = 20)

# Calculate legend point sizes with adjustment factor
legend_sizes <- scale_size(legend_degrees) * 0.5

# Add node size legend with regular intervals
legend("bottomleft",
       legend = paste(legend_degrees, "connections"),
       pt.cex = legend_sizes,
       pch = 21,
       pt.bg = adjustcolor("lightblue", alpha.f = 0.7),
       pt.lwd = 1,
       col = "darkblue",
       cex = 0.7,
       bty = "n",
       x.intersp = 0.5,
       y.intersp = 0.8,
       title = "Number of connections",
       title.adj = 0.15)

# Branches thickness and type based on relationships ----
# Load required library
library(igraph)

# Use KING coefficient for filtering
subset_data <- data[1:100, 108:110]
ids <- data$ID.x[1:100]

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

# Define edge styles based on relationship
edge_colors <- ifelse(edges$relationship == "First degree", 
                      adjustcolor("red", alpha.f = 0.6), 
                      adjustcolor("blue", alpha.f = 0.6))
edge_widths <- ifelse(edges$relationship == "First degree", 1, 0.2)
edge_types <- ifelse(edges$relationship == "First degree", 1, 2)  # 1 for solid, 2 for dashed

# Create custom layout with more spread
layout <- layout_with_fr(g)
layout <- layout * 1.5

# Create the plot
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.label.cex = 0.4,
     vertex.label = NA,
     vertex.color = adjustcolor("lightblue", alpha.f = 0.7),
     vertex.frame.color = "darkblue",
     edge.width = edge_widths,
     edge.color = edge_colors,
     edge.lty = edge_types,  # Add line type
     main = "Related Pairs Network\nNode size indicates number of connections")

# Add relationship legend
legend("bottom", 
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = c(1, 0.2),
       lty = c(1, 2),
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       bty = "n")

# Create regular intervals for legend
legend_degrees <- seq(0, 80, by = 20)
# Use exactly the same scaling as for the nodes
legend_sizes <- 2 + (legend_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 6 * 0.25  # Added 0.25 multiplier to match visual size

# Add node size legend
legend("bottomleft",
       legend = paste(legend_degrees, "connections"),
       pt.cex = legend_sizes,
       pch = 21,
       pt.bg = adjustcolor("lightblue", alpha.f = 0.7),
       pt.lwd = 1,
       col = "darkblue",
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       bty = "n")

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

# Color the nodes based on species ----

# Load required library
library(igraph)

# Use KING coefficient for filtering
subset_data <- data[1:1000, 108:110]
ids <- data$ID.x[1:1000]
# Trim whitespace from species names
species <- trimws(data$Species_corrected.x[1:1000])

# Create an empty edge list
edges <- data.frame(from = character(), to = character(), 
                    weight = numeric(), relationship = character())

# Create edges between IDs based on the KING values
for(i in 1:(length(ids)-1)) {
  for(j in (i+1):length(ids)) {
    king_value <- subset_data[j, "KING.snpduo"]
    
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

# Remove any rows with NA or duplicate edges
edges <- na.omit(edges)
edges <- unique(edges)

# Get unique vertices that are actually used in edges
used_vertices <- unique(c(edges$from, edges$to))

# Create a mapping between IDs and species for used vertices only
species_map <- species[match(used_vertices, ids)]

# Create the graph object
g <- graph_from_data_frame(edges, directed = FALSE, vertices = used_vertices)

# Calculate degree for each vertex
node_degrees <- degree(g)

# Normalize node sizes to be between 2 and 8
node_sizes <- 2 + (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 6

# Create color palette for unique species - using 4 distinct colors
unique_species <- sort(unique(species_map))
species_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")  # Red, Blue, Green, Purple
names(species_colors) <- unique_species

# Assign colors to nodes based on species
node_colors <- species_colors[species_map]

# Define edge styles based on relationship
edge_colors <- ifelse(edges$relationship == "First degree", 
                      adjustcolor("red", alpha.f = 0.6), 
                      adjustcolor("blue", alpha.f = 0.6))
edge_widths <- ifelse(edges$relationship == "First degree", 0.05, 0.02)
edge_types <- ifelse(edges$relationship == "First degree", 1, 2)

# Create custom layout with more spread
layout <- layout_with_fr(g)
layout <- layout * 1.5

# Set up the plotting area with extra space on the right
par(mar = c(5, 4, 4, 12), xpd = TRUE)

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
     main = "Related Pairs Network\nNode size: number of connections, Color: species")

# Add legends
# Relationship legend
legend(1.2, 1,  # Position on the right
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

# Species legend
legend(1.2, 0.4,  # Position on the right
       title = "Species",
       legend = unique_species,
       col = species_colors,
       pch = 19,
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       bty = "n",
       title.adj = 0.15)

# Create regular intervals for size legend
legend_degrees <- seq(0, 1000, by = 250)  # This will create: 0, 250, 500, 750, 1000
legend_sizes <- 2 + (legend_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 6 * 0.25

# Add node size legend
legend(1.2, -0.6,  # Position on the right
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

# Reset the plotting parameters to default
par(mar = c(5, 4, 4, 2), xpd = FALSE)

# Measure the number of first degree relationships in the dataset ----

# Load required library
library(igraph)

# Create an empty edge list
edges <- data.frame(from = character(), to = character(), 
                    weight = numeric(), relationship = character())

# Create edges between IDs based on the KING values
for(i in 1:(nrow(data)-1)) {
  for(j in (i+1):nrow(data)) {
    king_value <- data$KING.snpduo[j]
    
    if(!is.na(king_value)) {
      # Only keep first and second degree relationships
      if(king_value >= 0.0884) {
        relationship <- if(king_value >= 0.177) "First degree" else "Second degree"
        
        edges <- rbind(edges, data.frame(
          from = data$ID.x[i],
          to = data$ID.x[j],
          weight = king_value,
          relationship = relationship
        ))
      }
    }
  }
}

# Create summary of relationships
relationship_counts <- table(edges$relationship)

# Calculate percentages
relationship_percentages <- prop.table(relationship_counts) * 100

# Create a data frame with the results
relationship_summary <- data.frame(
  Relationship = names(relationship_counts),
  Count = as.numeric(relationship_counts),
  Percentage = round(as.numeric(relationship_percentages), 2)
)

# Print the results
print(relationship_summary)

# Calculate total number of individuals involved
n_individuals <- length(unique(c(edges$from, edges$to)))
print(paste("Total number of individuals with relationships:", n_individuals))

# Calculate number of individuals by species
individuals_by_species <- table(data$Species_corrected.x[match(unique(c(edges$from, edges$to)), data$ID.x)])
print("Number of individuals by species:")
print(individuals_by_species)