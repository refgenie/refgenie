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

build_taxonomy_tree <- function(taxonomies) {
  
  build_subtree <- function(data, level_col) {
    if(nrow(data) == 0) return(NULL)
    
    unique_taxa <- unique(data[[level_col]])
    unique_taxa <- unique_taxa[!is.na(unique_taxa)]
    
    if(length(unique_taxa) == 0) return(NULL)
    
    children <- map(unique_taxa, function(taxon) {
      subset_data <- data[data[[level_col]] == taxon & !is.na(data[[level_col]]), ]
      
      next_level <- case_when(
        level_col == 'domain' ~ 'kingdom',
        level_col == 'kingdom' ~ 'phylum',
        level_col == 'phylum' ~ 'class',
        level_col == 'class' ~ 'order',
        level_col == 'order' ~ 'family', 
        level_col == 'family' ~ 'genus',
        level_col == 'genus' ~ 'species',
        TRUE ~ NA_character_
      )
      
      # Create the node structure
      node <- list(
        name = taxon,
        taxonomicLevel = level_col
      )
      
      # Add genomeAlias for species (use NA instead of NULL)
      if(level_col == 'species') {
        node$genomeAlias <- NA_character_
      }
      
      # Add children if not at species level AND there are actually children
      if(!is.na(next_level) && level_col != 'species') {
        children_nodes <- build_subtree(subset_data, next_level)
        if(!is.null(children_nodes) && length(children_nodes) > 0) {
          node$children <- children_nodes
        }
        # If no children found, don't add the children property at all
      }
      
      return(node)
    })
    
    children <- children[!sapply(children, is.null)]
    return(children)
  }
  
  # Build tree starting from domain
  root_children <- build_subtree(taxonomies, 'domain')
  
  # Create root
  tree <- list(
    name = 'Root',
    children = root_children
  )
  
  return(tree)
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
write(json_output, '/Users/sam/Documents/Work/refgenie-ui/src/assets/taxa.json')

taxonomy_tree <- build_taxonomy_tree(taxonomies)
json_output <- toJSON(taxonomy_tree, pretty = TRUE, auto_unbox = TRUE)

write(json_output, '/Users/sam/Documents/Work/refgenie-ui/src/assets/tree.json')


# sample_genomes <- c(
#   "Homo sapiens",
#   "Mus musculus", 
#   "Rattus norvegicus",
#   "Danio rerio",
#   "Drosophila melanogaster",
#   "Caenorhabditis elegans",
#   "Saccharomyces cerevisiae",
#   "Arabidopsis thaliana",
#   "Escherichia coli",
#   "Gallus gallus",
#   "Sus scrofa",
#   "Bos taurus",
#   "Canis lupus familiaris",
#   "Felis catus",
#   "Macaca mulatta",
#   "Pan troglodytes",
#   "Gorilla gorilla",
#   "Pongo abelii",
#   "Equus caballus",
#   "Ovis aries",
#   "Capra hircus",
#   "Rattus rattus",
#   "Cricetulus griseus",
#   "Mesocricetus auratus",
#   "Cavia porcellus",
#   "Oryctolagus cuniculus",
#   "Xenopus tropicalis",
#   "Xenopus laevis",
#   "Takifugu rubripes",
#   "Tetraodon nigroviridis",
#   "Oryzias latipes",
#   "Gasterosteus aculeatus",
#   "Anopheles gambiae",
#   "Aedes aegypti",
#   "Culex quinquefasciatus",
#   "Bombyx mori",
#   "Tribolium castaneum",
#   "Apis mellifera",
#   "Nasonia vitripennis",
#   "Drosophila simulans",
#   "Drosophila yakuba",
#   "Drosophila erecta",
#   "Drosophila ananassae",
#   "Drosophila pseudoobscura",
#   "Drosophila virilis",
#   "Anopheles darlingi",
#   "Anopheles funestus",
#   "Strongylocentrotus purpuratus",
#   "Ciona intestinalis",
#   "Ciona savignyi",
#   "Branchiostoma floridae",
#   "Petromyzon marinus",
#   "Callorhinchus milii",
#   "Latimeria chalumnae",
#   "Lepisosteus oculatus",
#   "Oryza sativa",
#   "Zea mays",
#   "Sorghum bicolor",
#   "Triticum aestivum",
#   "Hordeum vulgare",
#   "Glycine max",
#   "Medicago truncatula",
#   "Lotus japonicus",
#   "Populus trichocarpa",
#   "Vitis vinifera",
#   "Solanum lycopersicum",
#   "Solanum tuberosum",
#   "Nicotiana tabacum",
#   "Brachypodium distachyon",
#   "Setaria italica",
#   "Physcomitrella patens",
#   "Chlamydomonas reinhardtii",
#   "Volvox carteri",
#   "Ostreococcus lucimarinus",
#   "Thalassiosira pseudonana",
#   "Phaeodactylum tricornutum",
#   "Emiliania huxleyi",
#   "Schizosaccharomyces pombe",
#   "Candida albicans",
#   "Neurospora crassa",
#   "Aspergillus nidulans",
#   "Aspergillus fumigatus",
#   "Magnaporthe oryzae",
#   "Ustilago maydis",
#   "Cryptococcus neoformans",
#   "Pneumocystis jirovecii",
#   "Plasmodium falciparum",
#   "Plasmodium vivax",
#   "Toxoplasma gondii",
#   "Leishmania major",
#   "Trypanosoma brucei",
#   "Trypanosoma cruzi",
#   "Giardia lamblia",
#   "Trichomonas vaginalis",
#   "Entamoeba histolytica",
#   "Dictyostelium discoideum",
#   "Nematostella vectensis",
#   "Hydra vulgaris",
#   "Trichoplax adhaerens",
#   "Amphimedon queenslandica",
#   "Monosiga brevicollis",
#   "Capsaspora owczarzaki"
# )
# 
# for (genome in sample_genomes) {
#   new_row <- get_ncbi_lineage(genome)
#   new_row$species <- genome
#   if (ncol(new_row) == ncol(taxonomies)) {
#     taxonomies <- rbind(taxonomies, new_row)
#   }
#   Sys.sleep(1)
# }
# 
# taxonomy_tree <- build_taxonomy_tree(taxonomies)
# json_output <- toJSON(taxonomy_tree, pretty = TRUE, auto_unbox = TRUE)
# 
# write(json_output, '/Users/sam/Documents/Work/refgenie-ui/src/assets/sample_tree.json')
# 
