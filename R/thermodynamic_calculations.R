###########################
#Thermodynamic calculations
###########################

#' Calculate primer melting temperature. 
#'
#' @param sequence A string corresponding to the primer sequence.
#' @param template_conc A numeric value corresponding to the DNA template 
#' concentration. Defauts to $10^-4$ M. 
#' @param ion_conc A string corresponding to the concentration of ion. See 
#' MELTING documentation for details. 
#'
#' @return A dataframe containing the primer sequence and the calculated melting
#' temperature \code{tm}.
#' @export
#'
#' @examples
melt_temperature <- function(sequence, template_conc, ion_conc) {
  melting_path <- "/Users/pedro_work/MELTING5.1.1/executable/melting"
  melting_command <- paste(melting_path, '-S', sequence, '-H dnadna', '-P',
                           template_conc, '-E', ion_conc)
  out <- system(melting_command, intern = T)
  
  #Creting melting temperature dataframe:
  melt_temp <- tibble(raw = out[5]) %>%
    mutate(
      primer = sequence,
      tm = as.numeric(str_extract(raw, '\\d+\\.\\d+'))
    ) %>%
    select(-raw)
  
  melt_temp
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
  file_prefix <- gsub('.in$', '', temp)
  
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
  system(command, intern = T)
  
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