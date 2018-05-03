devtools::load_all()

#Loading inputs:
target_input <- system.file('extdata', 'GFP sequence.fasta', package = "lambdaPrimeR")
vector_input <- system.file('extdata', 'pUC19.fasta', package = "lambdaPrimeR")
inputs <- read_inputs(target_input, vector_input)

#Creating primers
inputs <- create_primers(inputs)
primer_seqs <- get_sequences(inputs@primers)

tmd_values <- primer_seqs %>%
  group_by(orientation, primer_region) %>%
  mutate(tm = melt_temperature(toupper(seq)))
  do(tm = melt_temperature(.$seq)) 

vector <- inputs@vector
vector_seq <- vector$seq

  
#Thermodynamical calculations:
evaluate_primers <- function(run_object) {
  #Extracting primer object:
  primers <- run_object@primers
  #Get primer sequences:
  primer_sequences <- get_sequences(primers)
  
  #Extracting region of template sequence around the insertion position:
  vector <- run_object@vector
  insert_position <- vector$vector_anneal_left_end
  vector_seq <- substr(vector$seq, (insert_position - 49), (insert_position + 50))
  
  
  #Calculating melting temperatures:
  tmd_values <- primer_sequences %>%
    mutate(seq = toupper(seq)) %>%
    group_by(orientation, primer_region) %>%
    mutate(tm = melt_temperature(seq))

  #Annealing temparature
  tmd_values <- tmd_values %>%
    group_by(primer_region) %>%
    mutate(ta = min(tm) + 3)

  ta_vector <- unique(tmd_values$ta[tmd_values$primer_region=='vector'])
  
  #NUPACK calculations of complete sequences:
  seqs <- primer_sequences %>%
    filter(primer_region == 'complete') %>%
    pull(seq)
  names(seqs) <- c('fw', 'rv')
  
  #Self-dimer formation free energy:
  dimers <- tribble(
    ~primer, ~complex, ~complex_type, ~ta, ~dG,
    'forward', 'foward', 'sdimer', ta_vector, nupack_pfunc(c(seqs['fw'], seqs['fw']), ta_vector),
    'reverse', 'reverse', 'sdimer', ta_vector, nupack_pfunc(c(seqs['rv'], seqs['rv']), ta_vector),
    'forward', 'reverse', 'hdimer', ta_vector, nupack_pfunc(c(seqs['fw'], seqs['rv']), ta_vector),
    'reverse', 'forward', 'hdimer', ta_vector, nupack_pfunc(c(seqs['rv'], seqs['fw']), ta_vector)
  ) %>%
    unnest() %>%
    select(-partition_func)

  #MFE of individual primers
  mfe <- tribble(
    ~primer, ~complex, ~complex_type, ~ta, ~mfe,
    'forward', 'forward', 'hp', ta_vector, nupack_mfe(seqs['fw'], ta_vector),
    'reverse', 'reverse', 'hp', ta_vector, nupack_mfe(seqs['rv'], ta_vector)
  ) %>%
    unnest() %>%
    select(-structure, -base_pairs)

  #Primer-template complex free energy:
  primer_template <- tribble(
    ~primer, ~complex, ~complex_type, ~ta, ~dG,
    'forward', 'temp', 'ptemp', ta_vector, nupack_pfunc(c(seqs['fw'], vector_seq), ta_vector),
    'reverse', 'temp', 'ptemp', ta_vector, nupack_pfunc(c(seqs['rv'], vector_seq), ta_vector)
  ) %>%
    unnest() %>%
    select(-partition_func)
  
  #Creating output;
  results <- bind_rows(dimers, mfe, primer_template) %>%
    select(-complex) %>%
    spread(complex_type, dG)
  
  return(results)

  # return(tmd_values)
}

test <- evaluate_primers(inputs)


#Calculate primer score:
primer_score <- results %>%
  select(-complex, -partition_func) %>%
  spread(complex_type, dG) %>%
  mutate_at(vars(hdimer, hp, ptemp, sdimer), funs(as.numeric(.))) %>%
  mutate(z = ptemp - (sdimer + hdimer + hp))
