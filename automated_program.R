library(tidyverse)
library(meanShiftR)

predict_SV = function(data) {
  
  #Clean the data
  data_tidy <<- data %>%
    rename("Reads" = 1, "FLAG" = 2, "POS" = 4, "CIGAR" = 6, "PNEXT" = 8, "TLEN" = 9) %>% #Name columns
    select(Reads, FLAG, POS, CIGAR, PNEXT, TLEN) %>%
    filter(!is.na(POS) & !is.na(PNEXT)) %>% #Clean
    #Include only reads that have one mapping
    group_by(Reads) %>%
    mutate(count = n()) %>%
    filter(count == 2) %>%
    filter(TLEN > 0)
  
  print("Clean")
  
  data_marked <<- data_tidy %>%
    #Mark the data
    mutate(direction = ifelse(FLAG %in% c(161, 97, 163, 99) & TLEN > 0, "I", NA)) %>%
    mutate(direction = ifelse(FLAG %in% c(161, 97, 163, 99) & TLEN < 0, "O", direction)) %>%
    mutate(direction = ifelse(FLAG %in% c(81, 145, 83, 147) & TLEN < 0, "I", direction)) %>%
    mutate(direction = ifelse(FLAG %in% c(81, 145, 83, 147) & TLEN > 0, "O", direction)) %>%
    mutate(direction = ifelse(FLAG %in% c(65,129,113,177), "S", direction))

  print("Marked")
  
  match_length <- max(as.numeric(str_extract(data_marked$CIGAR, "[:digit:]+")))
  
  print(as.character(match_length))
  
  
  data_lengths <<- data_marked %>%
    mutate(PNEXT = ifelse(direction == "I" & TLEN > 0, PNEXT + match_length, PNEXT)) %>%
    mutate(POS = ifelse(direction == "I" & TLEN < 0, POS + match_length, POS)) %>% 
    mutate(PNEXT = ifelse(direction == "O" & TLEN < 0, PNEXT + match_length, PNEXT)) %>%
    mutate(POS = ifelse(direction == "O" & TLEN > 0, POS + match_length, POS)) %>%
    mutate(Length = abs(PNEXT - POS))
  
  print("Length Corrected")
  
  max_Length <- median(data_lengths$Length) + 3.22*mad(data_lengths$Length)

  
  print(as.character(max_Length))
  
  data_selected <- data_lengths %>%
    filter(Length > match_length*2) %>% #2*Read length
    filter(Length < max_Length) #Select empirical max distance
    
  coverage <- simulate_coverage(data_selected)
  
  min_coverage <- min(coverage)
  
  print("Coveraged Simulated")
    
  data_off_diagonal <- data_lengths %>%
    filter(Length > max_Length)
  
  points <- cbind(data_off_diagonal$POS, data_off_diagonal$PNEXT)
  
  candidates <- meanShift(points)
  
  candidates_points <- as.data.frame(candidates$value)
  
  candidates_points <- candidates_points %>%
    rename("POS" = 1, "PNEXT" = 2)
  
  data_off_diagonal$POS_center <- candidates_points$POS
  
  data_off_diagonal$PNEXT_center <- candidates_points$PNEXT
  
  print("Mean Shift Clustered")
  
  counts <<- data_off_diagonal %>%
    mutate(near = ifelse(abs(POS - POS_center) < max_Length & abs(PNEXT - PNEXT_center) < max_Length, TRUE, FALSE)) %>%
    filter(near) %>%
    group_by(POS_center, PNEXT_center) %>%
    mutate(count = sum(near)) %>%
    ungroup() %>%
    ungroup()
  
  median <-  data_off_diagonal %>%
    mutate(near = ifelse(abs(POS - POS_center) < max_Length & abs(PNEXT - PNEXT_center) < max_Length, TRUE, FALSE)) %>%
    filter(near) %>%
    group_by(POS_center, PNEXT_center) %>%
    summarize(count = n(), POS_center, PNEXT_center) %>%
    ungroup() %>%
    ungroup() %>%
    summarize(median = median(count)) %>%
    pull
  
  mad <- data_off_diagonal %>%
    mutate(near = ifelse(abs(POS - POS_center) < max_Length & abs(PNEXT - PNEXT_center) < max_Length, TRUE, FALSE)) %>%
    filter(near) %>%
    group_by(POS_center, PNEXT_center) %>%
    summarize(count = n(), POS_center, PNEXT_center) %>%
    ungroup() %>%
    ungroup() %>%
    summarize(mad = mad(count)) %>%
    pull
  
  sv_predictions <<- counts %>%
    mutate(median = median, mad = mad, relative_signal = count / mean(coverage)) %>%
    mutate(SV = ifelse(direction == "I", "Deletion", "None")) %>%
    mutate(SV = ifelse(direction == "O", "Duplication", SV )) %>%
    mutate(SV = ifelse(direction == "S", "Inversion", SV))


  print("SV Predicted")
  
  return(sv_predictions)
  
}

simulate_coverage = function(data) {
  
  data <- data %>%
    ungroup()
  
  n <- data %>%
    summarize(n()) %>%
    pull
  
  #Max d
  max_d <- data %>%
    summarise(max(Length)) %>%
    pull
  
  #Length of genome studied
  L <- data %>%
    summarise(max(PNEXT) - min(POS)) %>%
    pull
  
  
  #Run as a block to find empirical distribution
  
  simulationSums <- {} #Create an empty vector
  
  for (i in 1:10000) { #Run 1000 times to get a distribution
    vector <- {} #Create a holding vector
    for (d in 1:max_d) { #Iterate from 1 to max d
      poson <- rpois(1, (n / L)) #Poisson distribution with average n/L
      vector <- append(vector, sample(data$Length, size = poson, replace = TRUE) > d) #Sample a random amount of lengths and see if they are greater than current d
    }
    simulationSums[i] <- sum(vector) #Sum gives count of reads, append to sums
  }    
  
  return(simulationSums)
}





  
