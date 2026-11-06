suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(httr))
suppressMessages(library(jsonlite))
suppressMessages(library(xml2))


get_ncbi_lineage <- function(species_name) {
  search_url <- paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=', 
                       URLencode(species_name), '&retmode=json')
  
  search_req <- GET(search_url)
  search_res <- content(search_req, as = 'parsed')
  
  if (length(search_res$esearchresult$idlist) > 0) {
    taxid <- search_res$esearchresult$idlist[[1]]
    
    fetch_url <- paste0('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=', 
                        taxid, '&retmode=xml')
    
    fetch_req <- GET(fetch_url)
    fetch_res <- content(fetch_req, as = 'text', encoding = 'UTF-8') %>% read_xml()
    
    lineage_nodes <- xml_find_all(fetch_res, '//LineageEx/Taxon')
    
    lineage_df <- data.frame(list(
      rank = xml_text(xml_find_first(lineage_nodes, './/Rank')),
      taxonomy = xml_text(xml_find_first(lineage_nodes, './/ScientificName'))
    ))
    
    major_ranks <- c('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus')
    major_lineage <- lineage_df %>%
      filter(rank %in% major_ranks) %>% 
      pivot_wider(names_from = rank, values_from = taxonomy)

    return(major_lineage)
  }
}


refgenie_url <- 'https://api.refgenie.org/v4/species/summary'
refgenie_req <- GET(refgenie_url)
refgenie_res <- content(refgenie_req, as = 'text', encoding = 'UTF-8') %>% fromJSON()

genomes <- names(refgenie_res)

taxonomies <- data.frame(matrix(nrow = 0, ncol = 8)) %>%
  setnames(c('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'))

for (genome in genomes) {
  new_row <- get_ncbi_lineage(genome)
  if (!is.null(new_row)) {
    new_row$species <- genome
    if (ncol(new_row) == ncol(taxonomies)) {
      taxonomies <- rbind(taxonomies, new_row)
    }
  }
  Sys.sleep(1)
}

json_output <- toJSON(taxonomies, pretty = TRUE, auto_unbox = TRUE)
write(json_output, '../taxa.json')
