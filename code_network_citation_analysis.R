
# PRELIMINARY FUNCTIONS ########################################################

# Load the packages
# install.packages("sensobol") # Run this line if sensobol is not installed.

sensobol::load_packages(c("sensobol", "data.table", "tidyverse", "janitor", 
                          "igraph", "ggraph", "tidygraph", "cowplot", "viridis", 
                          "wesanderson", "parallel", "doParallel", "tm", "scales"))

# Custom theme for plots
theme_AP <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(fill = "transparent",
                                           color = NA),
          legend.key = element_rect(fill = "transparent",
                                    color = NA), 
          strip.background = element_rect(fill = "white"), 
          legend.margin = margin(0.5, 0.1, 0.1, 0.1),
          legend.box.margin = margin(0.2,-4,-7,-7), 
          plot.margin = margin(3, 4, 0, 4), 
          legend.text = element_text(size = 8), 
          axis.title = element_text(size = 10),
          legend.key.width = unit(0.4, "cm"), 
          legend.key.height = unit(0.4, "cm"), 
          legend.title = element_text(size = 9)) 
}

# Function to identify rows with NA
identify_rows_with_na <- function(data_table) {
  
  # Check for NA values in each row
  rows_with_na <- apply(data_table, 1, function(row) any(is.na(row)))
  
  # Return row indices with NA values
  return(which(rows_with_na))
}

# Define color palette
selected_wesanderson <- "Chevalier1"

################################################################################

# DIMENSIONS DATA ##############################################################

water.models <- c("WaterGAP", "PCR-GLOBWB", "LPJmL", "CLM4.5", "DBHM", 
                  "TOPMODEL", "H08", "JULES-W1", "MPI-HM", "VIC", "SWAT", 
                  "GR4J", "HYPE", "HBV", "MATSIRO", "SACRAMENTO", "MHM", 
                  "CWatM", "ORCHIDEE")

dt <- list()

for (i in 1:length(water.models)) {
  
  dt[[i]] <- fread(paste(water.models[[i]], ".csv", sep = ""), skip = 1) %>%
    clean_names() %>%
    data.table()
  
}

names(dt) <- water.models
dt.water <- rbindlist(dt, idcol = "Model")

wos.dt <- fread("final.dt.csv")
wos.titles <- wos.dt[Model %in% water.models]

# REMOVE DUPLICATED REFERENCES #################################################

# Number of papers in more than one model
n_occur <- data.frame(table(dt.water$publication_id))
papers_repeated <- data.table(n_occur[n_occur$Freq > 1,])
length(papers_repeated$Var1) # number of repeated papers

# Fraction of repeated papers over the total
length(papers_repeated$Var1) / nrow(dt.water)

# How many papers are repeated twice, three times, etc...
papers_repeated[, .(N.repeated.papers = .N), Freq]

# Extract which papers are repeated for which model
dt.sample.repeated <- dt.water[publication_id %in% papers_repeated$Var1] %>%
  .[, .(publication_id, Model, title, source_title_anthology_title)] %>%
  .[order(publication_id)]

dt.sample.repeated

# Randomly retrieve only one of the repeated studies per model
set.seed(6)
dt.no.repeated <- dt.sample.repeated[,.SD[sample(.N, min(1,.N))], publication_id]

# Setkey to filter and retrieve
res <- setkey(dt.water, publication_id, Model) %>%
  .[J(dt.no.repeated$publication_id, dt.no.repeated$Model)]

# Make the dataset without repeated papers across models
final.dt <- rbind(res, dt.water[!publication_id %in% papers_repeated$Var1])

# Check which papers do not have cited bibliography metadata and exclude them
final.dt <- final.dt[, empty_cited_references:= grepl("^\\s*$", cited_references)] %>%
  .[empty_cited_references == FALSE] %>%
  # Filter dataset to ensure all titles use a water model
  .[tolower(.$title) %in% wos.titles$title.large] %>%
  setnames(., "authors", "from.authors")

# Check the WOS and the Dimensions dataset
wos.dimensions <- merge(wos.dt[Model %in% water.models] %>%
  .[, .(WOS = .N), Model], 
  final.dt[, .(Dimensions = .N), Model], 
  by = "Model")

wos.dimensions[order(-Dimensions)]

# PLOT DIFFERENCES BETWEEN WOS AND DIMENSIONS ##################################

plot.models <- wos.dimensions %>%
  melt(., measure.vars = c("WOS", "Dimensions")) %>%
  ggplot(., aes(reorder(Model, value), value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  coord_flip() +
  scale_fill_manual(name = "Database", 
                    values = wes_palette(name = selected_wesanderson, 2)) +
  labs(y = "Count", x = "") +
  theme_AP() +
  theme(legend.position = c(0.8, 0.3))

wos.dimensions.dt <- wos.dimensions %>%
  melt(., measure.vars = c("WOS", "Dimensions"), variable.name = "dataset") %>%
  .[, .(total = sum(value)), dataset]

wos.dimensions.dt

plot.databases <- wos.dimensions.dt %>%
  ggplot(., aes(dataset, total, fill = dataset)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(name = "Database", 
                    values = wes_palette(name = selected_wesanderson, 2)) +
  labs(x = "", y = "Count") +
  theme_AP() + 
  theme(legend.position = "none")

plot_grid(plot.models, plot.databases, ncol = 2, rel_widths = c(0.65, 0.35), 
          labels = "auto")

# EXTRACT REFERENCES ###########################################################

setnames(final.dt, c("Model", "publication_id"), c("citing_model", "citing_id"))

column_names <- c("authors",
                  "author_id",
                  "source",
                  "year",
                  "volume",
                  "issue",
                  "pagination",
                  "to.doi",
                  "publication_id",
                  "times_cited")

direct.citation <- final.dt %>%
  .[, .(citing_id, cited_references, citing_model, title, doi, pub_year, from.authors)] %>%
  separate_rows(cited_references, sep = ";(?=\\[)") %>%
  separate(col = cited_references, into = column_names, sep = "\\|") %>%
  data.table() %>%
  setnames(., "authors", "to.authors")

# ARRANGE DATA FOR NETWORK CITATION ANALYSIS ###################################

# Create a directed graph from the dataset
edges <- data.table(from = direct.citation$citing_id, 
                    to = direct.citation$publication_id, 
                    from.model = direct.citation$citing_model, 
                    year = direct.citation$year, 
                    n.citations = direct.citation$times_cited, 
                    to.doi = direct.citation$to.doi, 
                    to.authors = direct.citation$to.authors)

# Merge data from citing papers with data from cited papers
network.dt <- merge(edges, final.dt[, .(citing_id, doi, pub_year, times_cited, from.authors)], 
                    by.x = "from", by.y = "citing_id")

colnames(network.dt)

# Change column names to clarify
new_colnames <- c("from", "to", "from.model", "to.year", "to.n.citations",
                  "to.doi", "to.authors", "from.doi", "from.year", "from.n.citations", "from.authors")
setnames(network.dt, colnames(network.dt), new_colnames)

# Reorder columns
setcolorder(network.dt, c("from", "to", "from.year", "to.year", "from.authors", 
                          "to.authors", "from.n.citations", "to.n.citations", 
                          "from.doi", "to.doi", "from.model"))

# Remove square brackets from the to.authors column
network.dt[, to.authors:= gsub("\\[|\\]", "", to.authors)]

# Identify the model of the cited paper
tmp <- network.dt[, .(from.model, from)] %>%
  unique()
setkey(tmp, "from")

# Define parallel computing
cl <- makeCluster(floor(detectCores() * 0.75))
registerDoParallel(cl)

# Search in parallel
to.model <- foreach(i=1:nrow(network.dt),
             .packages = "data.table", 
             .combine = "c") %dopar%
  {
    tmp[network.dt[[i, "to"]]]$from.model
  }

# Stop parallel cluster
stopCluster(cl)

# Add vector of model names
network.dt[, to.model:= to.model]

# Turn some columns into numeric
columns_to_modify <- grep("citation", names(network.dt), value = TRUE)
network.dt[, (columns_to_modify):= lapply(.SD, as.numeric), .SDcols = columns_to_modify]

# Export dataset
fwrite(network.dt, "network.dt.csv")

# PLOT DESCRIPTIVE FIGURES #####################################################

plot.n.citing <- network.dt[, unique(from.n.citations), from] %>%
  ggplot(., aes(V1)) +
  geom_histogram(color = "black", fill = "grey") + 
  scale_x_log10() +
  theme_AP() + 
  labs(x = "Citations", y = "Counts")

plot.n.cited <- network.dt[, unique(to.n.citations), to] %>%
  ggplot(., aes(V1)) +
  geom_histogram() + 
  scale_x_log10() +
  geom_histogram(color = "black", fill = "grey") + 
  theme_AP() + 
  labs(x = "Citations", y = "")

da <- list(plot.n.cited, plot.n.citing) 

plot_grid(plotlist = da, ncols = 2, labels = "auto")

# CHECK WHETHER WATER MODEL COMMUNITIES CITE THEIR OWN MODEL THE MOST ##########

da <- network.dt[, .N, .(from.model, to.model)] %>%
  .[!is.na(to.model)] 

keycol <-c("from.model","N")
setorderv(da, keycol, order = -1)
dcast(da, from.model~to.model, value.var = "N")

# Identify which model is the most cited
largest.citation <- da[da[, .I[N == max(N)], from.model]$V1]
sum(largest.citation$from.model == largest.citation$to.model) / nrow(largest.citation) 

# Tile plot
ggplot(da, aes(x = from.model, y = to.model, fill = N))+
  geom_tile(col = "black") +
  scale_fill_viridis() + 
  theme_AP() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# EXTRACT THE LARGEST AND THE SECOND LARGEST CITATION ##########################

merge(da[, .(maxima = max(N)), from.model], 
      da[, .(second.largest = sort(N, partial = length(N) - 1)[length(N)-1]), from.model], by = "from.model") %>%
  data.table() %>%
  .[, diff:= maxima / second.largest] %>%
  print()

# NETWORK CO-CITATION ANALYSIS #################################################

edges <- biblionetwork::biblio_cocitation(network.dt, 
                                          "from", 
                                          "to",
                                          weight_threshold = 3)
edges

references_filtered <- references_filtered %>% 
  relocate(publication_id, .before = authors) # first column has to be the identifier

graph <- tbl_main_component(nodes = direct.citation, 
                            edges = edges, 
                            directed = FALSE)
graph




















da[, sort(N), from.model]


fwrite(network.dt, "network.dt.csv")







# Calculate network metrics
citation_graph <- graph_from_data_frame(d = edges, directed = TRUE)

network_metrics <- data.table(node = V(citation_graph)$name,
                              degree = degree(citation_graph, mode = "in"),
                              betweenness = betweenness(citation_graph),
                              closeness = closeness(citation_graph),
                              pagerank = page_rank(citation_graph)$vector
)

network_metrics[order(-degree)]
network_metrics[order(-betweenness)]

# MORE NETWORK METRICS ---------------------------------------------------------

da <- tbl_graph(merge(edges, edges.dt, by = c("from", "to")), directed = TRUE)



da <- da[!identify_rows_with_na(da)]


da %>% 
  activate(nodes) %>% 
  mutate(assort = graph_assortativity(attr = from.model)) %>% 
  pull(assort) %>% 
  head(1)




# Merge metrics with the original dataset
dt.metrics <- merge(direct.citation, network_metrics, by.x = "citing_id", by.y = "node")
colnames(dt.metrics)

dt.metrics[, .(citing_model, citing_id, publication_id, doi, times_cited, year, 
               degree, betweenness, closeness, pagerank)] %>%
  .[order(citing_model)] 

# CO-CITATION NETWORK ##########################################################















edges <- biblionetwork::biblio_cocitation(edges, 
                                          source = "from", 
                                          ref = "to", 
                                          weight_threshold = 3)

references_filtered <- direct.citation %>% 
  distinct(publication_id, .keep_all = TRUE) %>% 
  select(-citing_id) %>%
  relocate(publication_id, .before = citing_model)

graph <- tbl_main_component(nodes = references_filtered, 
                            edges = edges, 
                            directed = FALSE)
graph




set.seed(1234)
graph <- add_clusters(graph, clustering_method = "leiden", 
                      objective_function = "CPM") # identifying clusters of nodes 

nb_communities <- graph %>% 
  activate(nodes) %>% 
  as_tibble %>% 
  select(Com_ID) %>% 
  unique %>% 
  nrow

palette <- scico::scico(n = nb_communities, palette = "hawaii") %>% # creating a color palette
  sample()

graph <- community_colors(graph, palette, community_column = "Com_ID")


# PLOT NETWORK #################################################################




graph <- as_tbl_graph(edges.without.na)

ggraph(graph, layout = "igraph", algorithm = "nicely") + 
  geom_edge_link(aes(colour = model)) + 
  geom_node_point(size = 0.5)













      

citation_graph <- graph_from_data_frame(d = edges.without.na, directed = TRUE)

# Calculate network metrics
network_metrics <- data.table(node = V(citation_graph)$name,
                              degree = degree(citation_graph, mode = "in"),
                              betweenness = betweenness(citation_graph),
                              closeness = closeness(citation_graph),
                              pagerank = page_rank(citation_graph)$vector
)

# Merge metrics with the original dataset
df_with_metrics <- merge(citations_dt, network_metrics, by.x = "citing_id", by.y = "node")


merge(citations_dt, network_metrics
# NETWORK CITATION ANALYSIS ####################################################










graph <- as_tbl_graph(edges.without.na)

ggraph(graph, layout = "sparse_stress") + 
  geom_edge_link(aes(colour = model)) + 
  geom_node_point(size = 0.5)













# Display the resulting dataset with network metrics


# Visualize the network
set.seed(123)  # For reproducibility

ggraph(citation_graph, layout = "fr") +
  geom_edge_link() +
  geom_node_point(aes(colour = citing_model)) +
  theme_void()























dimensions_references <- dimensions_direct_citation %>% 
  distinct(publication_id, .keep_all = TRUE) %>% 
  select(-citing_id)

prove.dt <- data.table(dimensions_direct_citation)[, .(citing_id, publication_id)]
prove.dt[, score:= 1]
spread(prove.dt, "citing_id", "score", fill = 0)

















# AFFILIATIONS -----------------------------------------------------------------

affiliations <- dimensions_dt[, .(publication_id,
                                  authors_affiliations_name_of_research_organization,
                                  authors_affiliations_country_of_research_organization)] %>%
  separate_rows(authors_affiliations_name_of_research_organization, sep = "; ")

# CLEAN REFERENCES -------------------------------------------------------------

setnames(dimensions_dt, "publication_id", "citing_id")

references_extract <- dimensions_dt %>%
  .[, .(citing_id, cited_references)] %>%
  separate_rows(cited_references, sep = ";(?=\\[)") %>%
  data.table()

column_names <- c("authors",
                  "author_id",
                  "source",
                  "year",
                  "volume",
                  "issue",
                  "pagination",
                  "doi",
                  "publication_id",
                  "times_cited")

dimensions_direct_citation <- references_extract %>% 
  separate(col = cited_references, into = column_names, sep = "\\|")


dimensions_references <- dimensions_direct_citation %>% 
  distinct(publication_id, .keep_all = TRUE) %>% 
  select(-citing_id)

prove.dt <- data.table(dimensions_direct_citation)[, .(citing_id, publication_id)]
prove.dt[, score:= 1]
spread(prove.dt, "citing_id", "score", fill = 0)


spread()


pg <- igraph::graph_from_data_frame(prove.dt)
plot(pg)


# CO-CITATION NETWORK ##########################################################

citations <- dimensions_direct_citation %>% 
  add_count(publication_id) %>% 
  select(publication_id, n) %>% 
  unique

references_filtered <- dimensions_references %>% 
  left_join(citations) %>% 
  filter(n >= 5)

edges <- biblionetwork::biblio_cocitation(filter(dimensions_direct_citation, publication_id %in% references_filtered$publication_id), 
                                          "citing_id", 
                                          "publication_id",
                                          weight_threshold = 3)
edges

references_filtered <- references_filtered %>% 
  relocate(publication_id, .before = authors) # first column has to be the identifier

graph <- tbl_main_component(nodes = references_filtered, 
                            edges = edges, 
                            directed = FALSE)
graph

library(tidygraph)
library(networkflow)

set.seed(1234)
graph <- add_clusters(graph, clustering_method = "louvain") # identifying clusters of nodes 

nb_communities <- graph %>% 
  activate(nodes) %>% 
  as_tibble %>% 
  select(ID) %>% 
  unique %>% 
  nrow
palette <- scico::scico(n = nb_communities, palette = "hawaii") %>% # creating a color palette
  sample()

graph <- community_colors(graph, palette, community_column = "Com_ID")

graph <- graph %>% 
  activate(nodes) %>%
  mutate(size = n,# will be used for size of nodes
         first_author = str_extract(authors, "(?<=\\[)(.+?)(?=,)"),
         label = paste0(first_author, "-", year)) 

graph <- community_names(graph, 
                         ordering_column = "size", 
                         naming = "label", 
                         community_column = "Com_ID")

graph <- vite::complete_forceatlas2(graph, 
                                    first.iter = 10000)


top_nodes  <- top_nodes(graph, 
                        ordering_column = "size", 
                        top_n = 15, 
                        top_n_per_com = 2,
                        biggest_community = TRUE,
                        community_threshold = 0.02)
community_labels <- community_labels(graph, 
                                     community_name_column = "Community_name",
                                     community_size_column = "Size_com",
                                     biggest_community = TRUE,
                                     community_threshold = 0.02)

# AUTHORS, AFFILIATIONS AND COUNTRY ############################################

dt.affiliations <- dt[, .(id, Title, Authors, Affiliations, `Authors with affiliations`)]
setnames(dt.affiliations, colnames(dt.affiliations), 
         c("id", "title", "Authors", "Affiliations", "Authors.Affiliations"))

dt.affiliations <- dt.affiliations %>%
  separate_longer_delim(., cols = "Authors.Affiliations", delim = ";") %>%
  data.table() %>%
  .[, `:=`(authors = sub(",.*", "", Authors.Affiliations), 
           title = tolower(title),
           department = sub("^[^,]*,\\s(.*?)(?:,\\s|$).*", "\\1", Authors.Affiliations),
           university = sub("^(?:[^,]*,\\s){2}(.*?)(?:,\\s|$).*", "\\1", Authors.Affiliations),
           country = sub(".*,\\s*([^,]+)$", "\\1", Authors.Affiliations))] %>%
  .[, .(id, authors, title, department, university, country)] 

# REFERENCES ###################################################################

dt.references <- dt[, .(id, References)] %>%
  separate_longer_delim(., cols = "References", delim = ";") %>%
  data.table()

da <- dt.references %>%
  mutate(authors = str_extract_all(References, extract_authors), 
         references = tolower(map(References, extract_longest_part)), 
         year = str_extract(References, "(?<=\\()([^()]*?)(?=\\)[^()]*$)"))


str(da)












da[4]




extract_longest_part <- function(input_string) {
  # Split the string by commas
  parts <- strsplit(input_string, ",")[[1]]
  
  # Filter out parts with commas and find the index of the longest remaining part
  longest_index <- which.max(nchar(parts[!grepl(",", parts)]))
  
  # Extract the longest remaining part
  longest_part <- parts[!grepl(",", parts)][longest_index]
  
  # Remove words within square brackets, including the brackets (case-insensitive)
  cleaned_part <- gsub("\\[[^\\]]*\\]", "", longest_part, ignore.case = TRUE)
  
  return(trimws(cleaned_part))
}


out <- list()
for(i in 1:nrow(da)) {
  
  out[[i]] <- extract_longest_part(da[[i, "References"]])
}

do.call(rbind, out)


da[714, ]
out[714]

# Example usage
input_string <- dt.references$References
longest_part <- extract_longest_part(input_string)

print(longest_part)







extract_title_rowwise <- function(references) {
  extract_title <- function(reference) {
    # Extract the part of the string after the last author's name
    last_author <- sub(".*\\b([^,]+),\\s*([^,]+)$", "\\1 \\2", reference)
    title_part <- sub(paste0(".*", last_author, ","), "", reference)
    
    # Use the previous logic to find the longest consecutive sequence without commas
    pattern <- "([^,]+)"
    matches <- regmatches(title_part, gregexpr(pattern, title_part))[[1]]
    
    if (length(matches) > 0) {
      # Filter out matches that contain digits (assumption: title does not contain digits)
      non_digit_matches <- matches[!grepl("\\d", matches)]
      if (length(non_digit_matches) > 0) {
        longest_sequence <- non_digit_matches[which.max(nchar(non_digit_matches))]
        return(trimws(longest_sequence))
      }
    }
    
    return(NULL)
  }
  
  titles <- apply(references, 1, extract_title)
  return(titles)
}


titles <- extract_title_rowwise(da$References)
print(titles)



extract_authors <- "\\b[A-Z][a-z]+\\s[A-Z]\\."
extract_year_brackets <- "(?<=\\()\\d{4}(?=\\))"
extract_pages <- "(?<= (p)?p\\. )([A-Z])?\\d+(-([A-Z])?\\d+)?"
extract_volume_and_number <- "(?<=( |^)?)\\d+ \\(\\d+(-\\d+)?\\)"




extract_title <- function(reference) {
  pattern <- "([^,]+)"
  matches <- regmatches(reference, gregexpr(pattern, reference))[[1]]
  
  if (length(matches) > 0) {
    longest_sequence <- matches[which.max(nchar(matches))]
    return(trimws(longest_sequence))
  } else {
    return(NULL)
  }
}


split_delim <- "[[:punct:]]"

split_fun <- function(x, split_delim) {
  
  st <- NA
  x1 <- as.character(x)
  x2 <- collapse_chars(x = x1, sep = sep)
  regex_fd <- split_delim
  x3 <- unlist(strsplit(x2, split = regex_fd, perl = TRUE))
  x4 <- unlist(strsplit(x3, split = "^( ){1,}"))
  x5 <- unlist(strsplit(x4, split = "$( ){1,}"))
  st <- x5[x5 != ""]
  return(st)

}
da$authors

da[1, ] %>%
  mutate(no_authors = str_remove(References, authors))

remove_authors_from_references(authors_string <- da[1, "authors"], 
                               references_string = da[1, "References"])

remove_authors_from_references <- function(authors_string, references_string) {
  # Split the authors and references strings into vectors
  authors <- unlist(strsplit(gsub(",", "", authors_string), " "))
  references <- strsplit(references_string, ",")[[1]]
  
  # Remove common authors from references
  cleaned_references <- gsub(paste0("\\b", authors, "\\b", collapse = "|"), "", references)
  
  # Remove leading and trailing commas and spaces
  cleaned_references <- trimws(gsub("^,|,$", "", cleaned_references))
  
  return(cleaned_references)
}





di <- split_fun(da[14, ], split_delim = ",")
di[which.max(nchar(di))]









str(da)
max(da)

words <- unlist(strsplit(gsub("[,.]", "", input_string), " "))

# Example usage:
input_string <- "Drouen L., Charpentier J.F., Semail E., Clenet S., Study of an innovative electrical machine fitted to marine current turbines, OCEANS 2007-Europe, pp. 1-6, (2007)"
result <- extract_largest_string_without_commas(input_string)
print(result)








sub(paste0(".*", da$authors, "(.*?),.*"), "\\1", da$References)

extract_text_after_pattern_until_comma(input_string = da$References[[1]], pattern = da$authors[[1]])


extract_text_after_pattern_until_comma <- function(input_string, pattern) {
  # Find the position of the first occurrence of the pattern
  pattern_position <- regexpr(pattern, input_string)[[1]]
  
  if (pattern_position > 0) {
    # Find the position of the first comma after the pattern
    first_comma_position <- regexpr(",", input_string[pattern_position:nchar(input_string)])[1]
    
    if (first_comma_position > 0) {
      # Extract text after the pattern until the first comma
      result <- substr(input_string, pattern_position + attr(first_comma_position, "match.length"), pattern_position + first_comma_position - 2)
      
      # Remove leading and trailing whitespaces
      return(trimws(result))
    }
  }
  
  return(NULL)
}

# Example usage:
text_string <- "John D. Smith, Spencer M. Johnson, Alice K. Brown, James P. White"
pattern <- "Spencer M\\."

result <- extract_text_after_pattern_until_comma(text_string, pattern)
print(result)


extract_text_after_pattern_until_comma <- function(input_string, pattern) {
  # Find the match details of the first occurrence of the pattern
  pattern_match <- regexec(pattern, input_string)
  
  if (pattern_match[[1]][1] > 0) {
    # Find the position of the first comma after the pattern
    first_comma_position <- regexpr(",", input_string, start = pattern_match[[1]][2] + attr(pattern_match, "match.length"))[1]
    
    if (first_comma_position > 0) {
      # Extract text after the pattern until the first comma
      result <- substr(input_string, pattern_match[[1]][2] + attr(pattern_match, "match.length"), pattern_match[[1]][2] + first_comma_position - 2)
      
      # Remove leading and trailing whitespaces
      return(trimws(result))
    }
  }
  
  return(NULL)
}

# Example usage:
text_string <- "John D. Smith, Spencer M. Johnson, Alice K. Brown, James P. White"
pattern <- "Spencer M\\."

result <- extract_text_after_pattern_until_comma(text_string, pattern)
print(result)







