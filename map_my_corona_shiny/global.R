library(dplyr)
library(leaflet)
#------------------------------------------------------------------------------#
#           read data - Country Polygons GeoJSON as sp                         #
#------------------------------------------------------------------------------#
# Source
# https://datahub.io/core/geo-countries
# https://datahub.io/core/geo-countries/r/countries.geojson
countries <- readRDS("data/countries.RDS")


color_area_IDs <- c(col_collect = "Collection date of top hit in each country",
                    col_release = "Release date of top hit in each country",
                    none        = "none")



validate_fasta <- function(path_fa, type){
  x <- readLines(path_fa)
  #print(x)
  f <- substr(x, 1, 1)
  # Check if it has header
  if (f[1] != ">") {
    stop("Sequence is not a valid fasta sequence")
  }
  
  # Check if nuc or prot
  
  # Remove header
  x <- x[!substr(x, 1, 1) == ">"]
  # collapse and get unique letters
  x <- paste(x, collapse = "")
  x <- unique(strsplit(x, "")[[1]]) 
  
  if (type == "nucleotide") {
    nucs <- c("A", "C", "G", "T", "U", "R", "Y", "K", "M", 
              "S", "W", "B", "D", "H", "V", "N", "-")
  } else if (type == "protein") {
    nucs <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", 
              "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "Y",
              "Z", "X", "*", "-")
  }
  
  
  if (all(x %in% nucs)) {
    print("valid sequence")
    TRUE
  } else {
    print("Not a valid sequence")
    FALSE
  }
  
}


