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
     vertex.frame.color = "darkblue",
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
subset_data <- data[1:30, 108:110]
ids <- data$ID.x[1:30]

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

# Create the plot
plot(g,
     layout = layout_with_fr(g),
     vertex.size = 8,
     vertex.label.cex = 0.6,
     vertex.color = "lightblue",
     vertex.frame.color = "darkblue",
     edge.width = 2,
     edge.color = edge_colors,
     main = "Related Pairs Network (First/Second Degree Only)")

# Add legend
legend("bottomright", 
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = 2)


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

# Create the plot
plot(g,
     layout = layout_with_fr(g),
     vertex.size = 8,
     vertex.label.cex = 0.6,
     vertex.color = "lightblue",
     vertex.frame.color = "darkblue",
     edge.width = 2,
     edge.color = edge_colors,
     main = "Related Pairs Network (First/Second Degree Only)")

# Add legend
legend("bottomright", 
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = 2)

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

# Normalize node sizes to be between 2 and 8 (reduced from 3-15)
node_sizes <- 2 + (node_degrees - min(node_degrees)) / (max(node_degrees) - min(node_degrees)) * 6

# Define colors for relationship types with transparency
edge_colors <- ifelse(edges$relationship == "First degree", 
                      adjustcolor("red", alpha.f = 0.6), 
                      adjustcolor("blue", alpha.f = 0.6))

# Create custom layout with more spread
layout <- layout_with_fr(g)
# Multiply layout coordinates to spread nodes further apart
layout <- layout * 1.5

# Create the plot
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.label.cex = 0.4,
     vertex.color = adjustcolor("lightblue", alpha.f = 0.7),
     vertex.frame.color = "darkblue",
     edge.width = 1.5,
     edge.color = edge_colors,
     main = "Related Pairs Network\nNode size indicates number of connections")

# Add legend
legend("bottomright", 
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = 2,
       bg = "white")
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
edge_widths <- ifelse(edges$relationship == "First degree", 2, 1)
edge_types <- ifelse(edges$relationship == "First degree", 1, 2)  # 1 for solid, 2 for dashed

# Create custom layout with more spread
layout <- layout_with_fr(g)
layout <- layout * 1.5

# Create the plot
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.label.cex = 0.4,
     vertex.color = adjustcolor("lightblue", alpha.f = 0.7),
     vertex.frame.color = "darkblue",
     edge.width = edge_widths,
     edge.color = edge_colors,
     edge.lty = edge_types,  # Add line type
     main = "Related Pairs Network\nNode size indicates number of connections")

# Add smaller legend in the bottom right corner
legend("bottomright", 
       legend = c("First degree", "Second degree"),
       col = c("red", "blue"), 
       lwd = c(2, 1),
       lty = c(1, 2),
       cex = 0.6,  # Smaller text size
       inset = c(0.02, 0.02),  # Adjust position
       bg = adjustcolor("white", alpha.f = 0.8),  # Semi-transparent background
       box.lwd = 0.5)  # Thinner box border

# Color of nodes: contry of origin ----
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

# Define edge styles based on relationship
edge_colors <- ifelse(edges$relationship == "First degree", 
                      adjustcolor("gray30", alpha.f = 0.6), 
                      adjustcolor("gray60", alpha.f = 0.6))
edge_widths <- ifelse(edges$relationship == "First degree", 2, 1)
edge_types <- ifelse(edges$relationship == "First degree", 1, 2)

# Create custom layout with more spread
layout <- layout_with_fr(g)
layout <- layout * 1.5

# Create the plot
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.label.cex = 0.4,
     vertex.color = node_colors,
     vertex.frame.color = "gray30",
     edge.width = edge_widths,
     edge.color = edge_colors,
     edge.lty = edge_types,
     main = "Related Pairs Network\nNode size: number of connections, Color: country of origin")

# Add two legends
# Relationship legend
legend("bottomright", 
       title = "Relationship",
       legend = c("First degree", "Second degree"),
       col = c("gray30", "gray60"), 
       lwd = c(2, 1),
       lty = c(1, 2),
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       box.lwd = 0.5)

# Country legend
legend("topright", 
       title = "Country of Origin",
       legend = unique_countries,
       col = country_colors,
       pch = 19,
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       box.lwd = 0.5)

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
                      adjustcolor("gray30", alpha.f = 0.6), 
                      adjustcolor("gray60", alpha.f = 0.6))
edge_widths <- ifelse(edges$relationship == "First degree", 2, 1)
edge_types <- ifelse(edges$relationship == "First degree", 1, 2)

# Create custom layout with more spread
layout <- layout_with_fr(g)
layout <- layout * 1.5

# Create the plot
plot(g,
     layout = layout,
     vertex.size = node_sizes,
     vertex.label.cex = 0.4,
     vertex.color = node_colors,
     vertex.frame.color = "gray30",
     edge.width = edge_widths,
     edge.color = edge_colors,
     edge.lty = edge_types,
     main = "Related Pairs Network\nNode size: number of connections, Color: species")

# Add two legends
# Relationship legend
legend("bottomright", 
       title = "Relationship",
       legend = c("First degree", "Second degree"),
       col = c("gray30", "gray60"), 
       lwd = c(2, 1),
       lty = c(1, 2),
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       box.lwd = 0.5)

# Species legend
legend("bottomleft", 
       title = "Species",
       legend = unique_species,
       col = species_colors,
       pch = 19,
       cex = 0.6,
       inset = c(0.02, 0.02),
       bg = adjustcolor("white", alpha.f = 0.8),
       box.lwd = 0.5)
