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
  
  return(list(consensus_string, score))
}

seqs <- readDNAStringSet("seq_score.fasta")
start_ind <- c(15,5,10,2,4) # pocet start ind musi odpovidat poctu sekvenci!!!!
Score(start_ind, seqs, 6)

## task 2
NextLeaf <- function(s, t, k ){ 
  for (i in t:1) {          # loop od posledni dna seq (t je pocet sekvenci)
    if (s[i] < k) {         # muzem dat uz inkrement?? 
      s[i] <- s[i] + 1
      return(s)             
    }
    s[i] <- 1
  }
  return(s)
}
k<- 70 - 5
NextLeaf(start_ind, 5, k)

## task 3
BFMotifSearch <- function(DNA, t, n, l){
  s <- rep(1,t) #1 se opakuje t-krat
  bestScore <- Score(s, DNA, l)$score # vytahne primo score, mam tam i ten consensus string 
  bestMotif <- s #?? nevim
  repeat {
    
    s <- NextLeaf(s, t, n - l + 1)
    
    # Compute score for current s
    currentScore <- Score(s, DNA, l)$score
    
    # update skore pokud je to potreba
    if (currentScore > bestScore) {
      bestScore <- currentScore
      bestMotif <- s
    }
    
    # stop pokud jsou to zas jen jednicky
    if (all(s == 1)) {
      break
    }
  }
  return(bestMotif)
}

## task 4
NextVertex <- function(s, i, t, k) {
  #i = level of vertex, t= num of sequences
  if (i < t) {
    s[i+1] <- 1
    return(list(s = s, i = i + 1))
  } else { #inkrement ze zadu
    for (j in t:1) {
      if (s[j] < k) {
        s[j] <- s[j] + 1
        return(list(s = s, i = j))
        }
      }
      return(list(s = s, i = 0))
  }
}

## task 5
ByPass <- function(s, i, t, k) {
  for (j in i:1) {
    if (s[j] < k){
      s[j] <- s[j] + 1
      return(list(s =s , i=j))
    }
  }
  return(list(s = s, i = 0))
}

## task 6

