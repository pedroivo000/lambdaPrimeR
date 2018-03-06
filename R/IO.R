#############
#Input/Output
#############


#'Read sequences from file.
#'
#'\code{read_sequences} reads a FASTA file into a dataframe, one sequence per
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
#'
#' @examples
#' read_sequences("data/test.fasta", input_type = 'target')
#' read_sequences("data/test.fasta", id_field = 'name', separator = '-')
read_sequences <- function(sequence_input, input_type, input_df=NULL,
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