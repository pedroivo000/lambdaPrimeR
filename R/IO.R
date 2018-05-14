#############
#Input/Output
#############

#' @include Template_classes.R Primers.R
NULL

#'Read sequences from FASTA file.
#'
#'\code{read_sequences} reads a FASTA file into a dataframe, one sequence per
#'row. This function also adds a \code{type} tag to each sequence, indicating if
#'they are either a \code{target} gene or a \code{vector} sequence. Depending on
#'the value of \code{type}, an instance of either the \code{Target} or 
#'\code{Vector} class will be created. User can optionally select what part of 
#'the sequence header will be used as the sequence \code{name}.
#'
#'@param file A path to a FASTA file, a vector containing the file path or a
#'  string containing the sequence with header.
#'@param input_type A string used to indicate if the sequence is a 'target' or a
#'  'template'.
#'@param header_structure A character vector corresponding to the structure of 
#'sequence header line. Each vector element should correspond to a field in the 
#'header line, delimitated by \code{sep}.
#'@param name_field A string inticating what field of the header string should 
#'be used as the sequence \code{name}.
#'@param sep A single character or regex used to separate fields in the sequence
#'  header according to the elements of \code{header_structure}. Defaults to '\\|'.
#'
#'@return A dataframe.
#'@export
#'
#' @examples
#' header <- c("name", "id", "extra_info")
#' read_sequences("data/test.fasta", input_type = 'target')
#' read_sequences("data/test.fasta", input_type = 'vector', 
#'                header_structure = header, separator = '-')
read_sequences <- function(file, input_type = c('target', 'vector'),
                           header_structure = NULL, sep = NULL, 
                           name_field = NULL) {
  
  #Checking input type:
  if(is.null(input_type)) {
    stop("Input type is required.")
  }
  type <- match.arg(input_type)
  
  #Read input with differnent functions depending on type of file (fasta or csv):
  if(grepl('*.fasta', file)) {
    df <- read_sequences_fasta(file,input_type, header_structure, sep, name_field)
  } 
  else if (grepl('*.csv', file)) {
    df <- read_sequences_csv(file)
  }
  
  
  #Checking what type of input was provided and creating appropriate class obj.:
  object <- switch (input_type,
                    'target' = Target(df),
                    'vector' = Vector(df)
  )
  
  return(object)
}

read_sequences_fasta <- function(file, input_type = c('target', 'vector'),
                                 header_structure = NULL, sep = NULL, 
                                 name_field = NULL) {
  
  #Checking input type:
  if(is.null(input_type)) {
    stop("Input type is required.")
  }
  type <- match.arg(input_type)
  
  #Reading fasta file into a column:
  raw <- tibble(text = read_file(file))
  
  #Converting the raw text column into the correct format:
  df <- raw %>%
    separate_rows(text, sep = '\n>') %>%
    separate(text, into = c('header', 'seq'), sep = '\n', extra = 'merge') %>%
    mutate(
      path = file, 
      type = input_type,
      seq = stringr::str_replace_all(seq, '\n', ''),
      header = stringr::str_replace(header, '>', ''),
      length = nchar(seq), 
      id = row_number()
    )
  
  #Checking if optional parameters to parse the header string were passed:
  if(is.null(header_structure)) {
    header_structure <- 'name'
  }
  if(is.null(sep)) {
    sep <- '\\|'
  }
  if(is.null(name_field)) {
    name_field <- 'name'
  }
  
  #Spliting header string:
  df <- df %>%
    separate(header, into = header_structure, remove = F, extra = 'drop', sep = sep) 
  
  #Renaming name_field column to "name": 
  colnames(df)[colnames(df)==name_field] <- "name"
  
  #Removing extra columns in df:
  df <- df %>%
    select(id, name, path, type, header, seq, length)
  
  return(df)
}

read_sequences_csv <- function(file, input_type = c('target', 'vector')) {
  #Saving file path:
  filepath <- file
  
  #Checking input type:
  if(is.null(input_type)) {
    stop("Input type is required.")
  }
  type <- match.arg(input_type)
  
  # #Checking if sequence csv file is valid:
  # #Required columns: name, header, seq (type col. is optional)
  # errors <- check_sequence_csv(file)
  # if()
  
  #Reading csv file:
  df <- read_csv(file)
  
  #Checking if header column if present, otherwise add it:
  if(!'header' %in% colnames(df)) {
    df <- df %>%
      mutate(header = name)
  }
  #Checking if input type column was present in the original csv file:
  if(!'type' %in% colnames(df)) {
    df <- df %>%
      mutate(type = type)
  }
  #Calulate sequence length and assigne unique id to each sequence:
  df <- df %>%
    mutate(
      length = nchar(seq), 
      id = row_number(),
      path = filepath
    ) 
  
  return(df)
}

read_primers_csv <- function(file) {
  
  #Reading csv file:
  df <- read_csv(file)
  
  #Reorganizing columns to match required columns of Primer objects: 
  df <- df %>%
    gather(seq_type, seq, -target, -vector, -orientation) %>%
    unite(name, orientation, seq_type, sep = '_') %>%
    spread(name, seq) %>%
    rename(
      target_region_origin = target, 
      vector_region_origin = vector
    ) %>%
    mutate(
      id = row_number(),
      forward_length = nchar(forward_primer_seq),
      reverse_length = nchar(reverse_primer_seq)
    )
  
  #Creating Primer object with df:
  primer_object <- Primers(df)
  return(primer_object)
}


#'Read sequences from file.
#'
#'\code{read_sequences_noclass} reads a FASTA file into a dataframe, one sequence per
#'row. This function also adds a \code{type} tag to each sequence, indicating if
#'they are either a \code{target} gene or a \code{template} sequence. User can
#'optionally select what part of the sequence header will be used as internal
#'id.
#'
#'@param file A path to a FASTA file, a vector containing the file path or a
#'  string containing the sequence with header.
#'@param input_type A string used to indicate if the sequence is a 'target' or a
#'  'template'.
#'@param input_df A dataframe containing an already imported sequence input. New
#'  inputs will be appended to this object.
#'@param id_field A string with the name of the id field. Defaults to 'id'.
#'@param separator A single character used to separate fields in the sequence
#'  headers. Defaults to '\\|'.
#'
#'@return A dataframe.
#'@export
#'
#' @examples
#' read_sequences_noclass("data/test.fasta", input_type = 'target')
#' read_sequences_noclass("data/test.fasta", id_field = 'name', separator = '-')
read_sequences_noclass <- function(sequence_input, input_type, input_df=NULL,
                           id_field = NULL, separator = NULL) {
  raw <- tibble(text = read_file(sequence_input))
  
  df <- raw %>%
    separate_rows(text, sep = '\n>') %>%
    separate(text, into = c('header', 'seq'), sep = '\n', extra = 'merge') %>%
    mutate(
      type = input_type,
      seq = stringr::str_replace_all(seq, '\n', ''),
      header = stringr::str_replace(header, '>', ''),
      length = nchar(seq)
    )
  
  #Checking if optional parameters were passed
  if(is.null(id_field)) {
    id_field <- 'id'
  }
  if(is.null(separator)) {
    separator <- '\\|'
  }
  
  #Creating sequence id from seq header
  df <- df %>%
    separate(header, into = id_field, remove = F, extra = 'drop', sep = separator)
  
  #Appending new input to input dataframe if input_df object was provided:
  if(!is.null(input_df)) {
    input_df <- bind_rows(input_df, df)
  } else {
    df
  }
}


#' Append sequence input to global input dataframe
#' 
#' Internal use in the Shiny UI only.
#'
#' @param input 
#' @param global_value 
#'
#' @return
#' 
#'
#' @examples
collect_input_shiny <- function(input, global_value) {
  new_input <- input
  global_value <- rbind(global_value, new_input)
}



#03/03/18, Wife's comment here:
#love you