rm(list=ls())

library(Biostrings)

#task 1
Score <- function(start_ind, DNA_string, motif_length) {
  n <- length(DNA_string) #number of sequences
  
  # alignment matrix
  motifs <- character(n)
  for (i in seq_len(n)) {
    motifs[i] <- as.character(subseq(DNA_string[[i]], start = start_ind[i], width = motif_length))
  }
  
  #rows = sequences, columns = motif positions
  motif_matrix <- do.call(rbind, strsplit(motifs, "")) # do.call = (funkce, na co ji volame)-- treba (sum, values))
  
  # inicializace profile matrix: rows = A,C,G,T ; columns = motif positions
  bases <- c("A", "C", "G", "T")
  profile <- matrix(0, nrow = 4, ncol = motif_length, dimnames = list(bases, NULL))
  
  for (j in seq_len(motif_length)) {
    columns <- motif_matrix[,j]
    counts <- table(columns)
    for (b in bases){
      #ifelse returns a value with the same shape as test which is 
      #filled with elements selected from either yes or no depending on whether the element of test is TRUE or FALSE.
      profile[b, j] <- ifelse(b %in% names(counts), counts[[b]], 0) 
    }
    #profile ma radky A C G T a sloupce jako "pozice v motivu" -- cisla v matici znamenaji vyskyt base na tom miste
    
  }
  consensus <- character(motif_length)
  score <- 0
  for (j in seq_len(motif_length)) {
    col_values <- profile[, j]
    max_base <- names(which.max(col_values))  # písmeno s max počtem
    consensus[j] <- max_base
    score <- score + max(col_values)          # přičte počet výskytů
  }
  consensus_string <- paste0(consensus, collapse = "")
  
  return(consensus_string)
}

