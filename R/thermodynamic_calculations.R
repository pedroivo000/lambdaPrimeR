###########################
#Thermodynamic calculations
###########################


#' Calculate primer binding free energies.
#' 
#' This function is used to calculate the binding free energy of each primer in
#' primer pair for different complexes: hairping formation, homodimer (primer
#' binding to itself), heterodimer (forward and reverse primer binding to each
#' other) and the primer-vector complex.
#' 
#'
#' @param run_object An object of the Run class, containing the target and vector
#' inputs. 
#'
#' @return A dataframe containing the calculated free energy for each primer's
#' hairpin, homodimer, heterodimer and primer-vector complex. 
#' 
#' @export
#'
#' @examples
evaluate_primers <- function(run_object) {
  #Extracting primer object:
  primers <- run_object@primers
  #Get primer sequences:
  # primer_sequences <- get_sequences(primers)
  
  #Extracting region of template sequence around the insertion position:
  vector <- asS3(run_object@vector) %>%
    ungroup() %>%
    distinct(id, name, seq)
  
  vector_seqs <- primers %>%
    select(id, vector_region_origin) %>%
    group_by(id) %>%
    mutate(
      vector_seq = vector$seq[vector$name == vector_region_origin],
      vector_seq = substr(vector_seq, 
                          primers$forward_vector_anneal_beg - 29, 
                          primers$reverse_vector_anneal_end + 30
    ))
  
  primers <- left_join(primers, vector_seqs)
  
  #Evaluating each primer pair:
  primer_eval <- primers %>%
    group_by(id) %>%
    do(result = evaluate_primer_pair(.)) %>%
    unnest()
  
  # #Calculating melting temperatures:
  # tmd_values <- primers %>%
  #   select(id, contains('seq')) %>%
  #   gather(seqtype, seq, -id) %>%
  #   mutate(seq = toupper(seq)) %>%
  #   group_by(id, seqtype) %>%
  #   mutate(tm = melt_temperature(sequence = seq))
  # 
  # #Annealing temperature:
  # tmd_values <- tmd_values %>%
  #   separate(seqtype, into = c('orientation', 'primer_region'), 
  #            sep = '_', extra ='merge') %>%
  #   group_by(id, orientation, primer_region) %>%
  #   mutate(ta = min(tm) + 3)
  # 
  # # tmd_values <- primer_sequences %>%
  # #   mutate(seq = toupper(seq)) %>%
  # #   group_by(orientation, primer_region) %>%
  # #   mutate(tm = melt_temperature(seq))
  # # 
  # # #Annealing temparature
  # # tmd_values <- tmd_values %>%
  # #   group_by(primer_region) %>%
  # #   mutate(ta = min(tm) + 3)
  # 
  # ta_vector <- unique(tmd_values$ta[tmd_values$primer_region=='vector_anneal_seq'])
  # 
  # #NUPACK calculations of complete sequences:
  # seqs <- primer_sequences %>%
  #   filter(primer_region == 'complete') %>%
  #   pull(seq)
  # names(seqs) <- c('fw', 'rv')
  # 
  # #Self-dimer formation free energy:
  # dimers <- tribble(
  #   ~primer, ~complex, ~complex_type, ~ta, ~dG,
  #   'forward', 'foward', 'sdimer', ta_vector, nupack_pfunc(c(seqs['fw'], seqs['fw']), ta_vector),
  #   'reverse', 'reverse', 'sdimer', ta_vector, nupack_pfunc(c(seqs['rv'], seqs['rv']), ta_vector),
  #   'forward', 'reverse', 'hdimer', ta_vector, nupack_pfunc(c(seqs['fw'], seqs['rv']), ta_vector),
  #   'reverse', 'forward', 'hdimer', ta_vector, nupack_pfunc(c(seqs['rv'], seqs['fw']), ta_vector)
  # ) %>%
  #   unnest() %>%
  #   select(-partition_func)
  # 
  # #MFE of individual primers
  # mfe <- tribble(
  #   ~primer, ~complex, ~complex_type, ~ta, ~mfe,
  #   'forward', 'forward', 'hp', ta_vector, nupack_mfe(seqs['fw'], ta_vector),
  #   'reverse', 'reverse', 'hp', ta_vector, nupack_mfe(seqs['rv'], ta_vector)
  # ) %>%
  #   unnest() %>%
  #   select(-structure, -base_pairs)
  # 
  # #Primer-template complex free energy:
  # primer_template <- tribble(
  #   ~primer, ~complex, ~complex_type, ~ta, ~dG,
  #   'forward', 'temp', 'ptemp', ta_vector, nupack_pfunc(c(seqs['fw'], vector_seq), ta_vector),
  #   'reverse', 'temp', 'ptemp', ta_vector, nupack_pfunc(c(seqs['rv'], vector_seq), ta_vector)
  # ) %>%
  #   unnest() %>%
  #   select(-partition_func)
  # 
  # #Creating output;
  # results <- bind_rows(dimers, mfe, primer_template) %>%
  #   select(-complex) %>%
  #   spread(complex_type, dG) %>%
  #   mutate_at(vars(hdimer, hp, ptemp, sdimer), funs(as.numeric(.))) %>%
  #   mutate(z = ptemp - (sdimer + hdimer + hp))
  
  return(primer_eval)
}

evaluate_primer_pair <- function(primer_df) {

  vector_seq <- primer_df$vector_seq
  
  #Calculating melting temperatures:
  tmd_values <- primer_df %>%
    select(contains('seq'), -vector_seq) %>%
    gather(seqtype, seq) %>%
    mutate(seq = toupper(seq)) %>%
    filter(!grepl('primer_seq', seqtype)) %>%
    group_by(seqtype) %>%
    mutate(tm = melt_temperature(sequence = seq))
  
  #Annealing temperature:
  tmd_values <- tmd_values %>%
    separate(seqtype, into = c('orientation', 'primer_region'), 
             sep = '_', extra ='merge') %>%
    group_by(primer_region) %>%
    mutate(ta = min(tm) + 3)
  
  ta_vector <- unique(tmd_values$ta[tmd_values$primer_region=='vector_anneal_seq'])
  
  #NUPACK calculations of complete sequences:
  fw <- primer_df$forward_primer_seq
  rv <- primer_df$reverse_primer_seq
  
  #Self-dimer formation free energy:
  dimers <- tribble(
    ~primer, ~complex, ~complex_type, ~ta, ~dG,
    'forward', 'foward', 'sdimer', ta_vector, nupack_pfunc(c(fw, fw), ta_vector),
    'reverse', 'reverse', 'sdimer', ta_vector, nupack_pfunc(c(rv, rv), ta_vector),
    'forward', 'reverse', 'hdimer', ta_vector, nupack_pfunc(c(fw, rv), ta_vector),
    'reverse', 'forward', 'hdimer', ta_vector, nupack_pfunc(c(rv, fw), ta_vector)
  ) %>%
    unnest() %>%
    select(-partition_func)
  # return(dimers)
  
  #MFE of individual primers
  mfe <- tribble(
    ~primer, ~complex, ~complex_type, ~ta, ~mfe,
    'forward', 'forward', 'hp', ta_vector, nupack_mfe(fw, ta_vector),
    'reverse', 'reverse', 'hp', ta_vector, nupack_mfe(rv, ta_vector)
  ) %>%
    unnest() %>%
    select(-structure, -base_pairs)
  
  #Primer-template complex free energy:
  primer_template <- tribble(
    ~primer, ~complex, ~complex_type, ~ta, ~dG,
    'forward', 'temp', 'ptemp', ta_vector, nupack_pfunc(c(fw, vector_seq), ta_vector),
    'reverse', 'temp', 'ptemp', ta_vector, nupack_pfunc(c(rv, vector_seq), ta_vector)
  ) %>%
    unnest() %>%
    select(-partition_func)
  
  #Creating output;
  results <- bind_rows(dimers, mfe, primer_template) %>%
    select(-complex) %>%
    spread(complex_type, dG) %>%
    mutate_at(vars(hdimer, hp, ptemp, sdimer), funs(as.numeric(.))) %>%
    mutate(z = ptemp - (sdimer + hdimer + hp))
  
  return(results)
}

evaluate_primers_noclass <- function(primers_df) {
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
    spread(complex_type, dG) %>%
    mutate_at(vars(hdimer, hp, ptemp, sdimer), funs(as.numeric(.))) %>%
    mutate(z = ptemp - (sdimer + hdimer + hp))
  
  return(results)
}

#' Calculate primer melting temperature. 
#'
#' @param sequence A string corresponding to the primer sequence.
#' @param template_conc A numeric value corresponding to the concentration of 
#' the DNA molecule in excess. During PCR experiments, this value correspond to 
#' the primer concentrations. Defauts to $10^-6$ M. 
#' @param ion_conc A string corresponding to the concentration of ion (see 
#' MELTING documentation for details). Defaults to \code{Mg=2e-03}.
#'
#' @return A dataframe containing the primer sequence and the calculated melting
#' temperature \code{tm}.
#' @export
#'
#' @examples
melt_temperature <- function(sequence, template_conc = 1e-06, ion_conc = 'Mg=2e-03') {
  melting_path <- "/Users/pedro_work/MELTING5.1.1/executable/melting"
  sequence <- toupper(sequence)
  melting_command <- paste(melting_path, '-S', sequence, '-H dnadna', '-P',
                           template_conc, '-E', ion_conc)
  #Running MELTING:
  out <- system(melting_command, intern = T)
  #Extracting melting temperature value from command output:
  melt_temp <- as.numeric(stringr::str_extract(out[5], '\\d+\\.*\\d*'))
  
  return(melt_temp)
}

#' Gather NUPACK's tools parameters.
#'
#' @param tool Name of the NUPACK tool.
#' @param input_sequence A string corresponding to the input sequence.
#' @param melting_temperature A numeric value indicating the melting temperature
#' to be used during the NUPACK's calculations.
#' @param multi Multiple input sequences? (Default is \code{FALSE}).
#' @param dangles A string with the way in which dangle energies are 
#' incorporated. Possible values: 
#' 
#' \itemize{
#'   \item \code{none}: No dangle energy are incorporated
#'   \item \code{some} (default): A dangle energy is incorporated for each 
#'   unpaired base flanking a duplex (a base flanking two duplexes contributes 
#'   only the minimum of the two possible dangle energies).
#'   \item \code{all}: A dangle energy is incorporated for each base flanking a 
#'   duplex regardless of whether it is paired.
#' }
#'
#' @return A data frame with base parameters for a NUPACK function
nupack_args <- function(tool, input_sequence, melting_temperature,
                        dangles=c('some', 'none', 'all')){
  
  #Creating temporary file for NUPACK:
  temp <- tempfile(fileext = '.in')
  
  #Checking number of input sequences:
  if(length(input_sequence) == 1) {
    #"-multi" flag:
    multi <- FALSE
    multi_opt <- ''
    
    #writing input to temporary file:
    write_file(input_sequence, temp)
  } else {
   
    #"-multi" flag:
    multi <- TRUE
    multi_opt <- "-multi"
    
    #Writing sequences on temporary files:
    n_seqs <- length(input_sequence)
    distinct_seqs <- paste0(1:n_seqs, collapse = ' ')
    file_content <- c(n_seqs, input_sequence, distinct_seqs)
    write_lines(file_content, temp)
  }
  
  #File prefix:
  file_prefix <- tools::file_path_sans_ext(temp)

  #dangles:
  dangles_opt <- match.arg(dangles)
  
  #Base command: 
  base <- paste0(c(tool, '-T', melting_temperature, '-material dna', '-dangles', 
                   dangles_opt, file_prefix, multi_opt), collapse = ' ')
  
  #Storing args in a data frame:
  args <- tibble(
    multi_bool = multi,
    base_command = base,
    prefix = file_prefix
  )
}

#' Calculate base-pairing probability using NUPACK's \code{pairs} tool.
#'
#' @param cutoff A numeric value (default 0.001). Only probabilities above 
#' \code{cutoff} are written to output. 
#' @inheritParams nupack_args
#'
#' @return A data frame with three columns: \code{p} corresponds to the 
#' probability of \code{ibase} being paired to \code{jbase}. If \code{paired} =
#' \code{FALSE}, column \code{p} contains the probability of 
#' \code{ibase} being unpaired. 
#' 
#' @export
#'
#' @examples
nupack_pairs <- function(input_sequence, melting_temperature, paired=TRUE, 
                         dangles=c('some', 'none', 'all'), cutoff = NULL) {
  
  #Parsing function variables to get NUPACK's base arguments:
  args <- nupack_args(tool = 'pairs', input_sequence, melting_temperature, dangles)
  
  #Probaility cutoff (unique parameter of pais function):
  if(is.null(cutoff)) {
    cutoff <- 0.001
  }
  
  #Asssembling command string:
  command <- paste(args$base_command, '-cutoff', cutoff, sep = ' ')
  
  #Executing NUPACK command:
  system(command, intern = T)
  
  #Import output file:
  #Depending on "-multi" flag, a different input file will be generated - suffix 
  #".ppairs" if multi=FALSE and ".epairs" if multi=TRUE.
  output_file <- ifelse(args$multi_bool, 
                        paste(args$prefix, ".epairs", sep = ''),
                        paste(args$prefix, ".ppairs", sep = ''))
  nupack_output <- read_lines(output_file)
  
  #Deleting temp files
  temp_files <- paste(args$prefix, ".*", sep = '')
  unlink(temp_files)
  
  #Parse output:
  #Get total number of bases input sequence(s):
  i <- grep('^\\d+$', nupack_output, perl = T) #get index
  n_bases <- as.numeric(nupack_output[i])

  #Base pairing probabilities:
  pair_probabilities <- tibble(lines = nupack_output) %>%
    filter(grepl('^\\d+\\t', lines)) %>%
    separate(lines, into = c('ibase', 'jbase', 'p'), sep = '\t') %>%
    mutate_all(as.numeric)
  
  #Renaming 'p' column to 'exp_num' if multi=TRUE
  if(args$multi_bool) {
    pair_probabilities <- pair_probabilities %>%
      rename(exp_num = p)
  }
  
  #Extracting base pair probabilities by default (otherwise retrieve probability
  #of bases being unpaired):
  if(paired) {
    pair_probabilities %>%
      filter(jbase < n_bases+1)
  } else {
    pair_probabilities %>%
      filter(jbase == n_bases+1)
  }
}

#' Minimum free energy (MFE) structure determination using NUPACK's \code{mfe} 
#' tool.
#'
#'
#' @param degenerate Compute all structures with MFE? Default is \code{FALSE}.
#' @inheritParams nupack_args
#'
#' @return A data frame with the MFE value and the dot-parens-plus drawing of 
#' MFE structure. 
#' @export
#'
#' @examples
nupack_mfe <- function(input_sequence, melting_temperature, 
                       dangles=c('some', 'none', 'all'), degenerate=FALSE) {
  
  #Parsing function variables to get NUPACK's base arguments:
  args <- nupack_args(tool = 'mfe', input_sequence, melting_temperature, dangles)
  degenerate_opt <- ifelse(degenerate, ' -degenerate', '')
  
  #Asssembling command string:
  command <- paste(args$base_command, degenerate_opt, sep = " ")
  
  #Executing NUPACK command:
  system(command)
  
  #Import output file 
  output_file <- paste(args$prefix, '.mfe', sep = '')
  nupack_output <- read_lines(output_file)
  
  #Deleting temp files
  temp_files <- paste(args$prefix, ".*", sep = '')
  unlink(temp_files)
  
  #Extracting base pairs:
  base_pairs <- tibble(lines = nupack_output) %>%
    filter(grepl('^\\d+\\t', lines)) %>%
    separate(lines, into = c('ibase', 'jbase')) 
  
  #Extract MFE value:
  mfe <- tibble(
    # seq = input_sequence,
    dG = nupack_output[15],
    structure = nupack_output[16],
    base_pairs = list(base_pairs)
  )
  
 return(mfe)
}

#' Calculate the partition function of a DNA complex using NUPACK's \code{pfunc}
#' tool.
#'
#' @inheritParams nupack_args
#'
#' @return A data frame containing the free energy of the complex and the
#' partition function. 
#' @export
#'
#' @examples
nupack_pfunc <- function(input_sequence, melting_temperature, 
                         dangles=c('some', 'none', 'all')) {
  
  #Parsing function variables to get NUPACK's base arguments:
  args <- nupack_args(tool = 'pfunc', input_sequence, melting_temperature, dangles)

  #Executing NUPACK command:
  nupack_output <- system(args$base_command, intern = T)
  
  #Extract dG and partion function:
  pfunc <- tibble(
    dG = nupack_output[14],
    partition_func = nupack_output[15]
  )
  
}