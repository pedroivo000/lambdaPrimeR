#############
#Input/Output
#############


#' Read sequences from file.
#'
#'\code{read_sequences} reads a FASTA file into a dataframe, one sequence per row.
#'This function also adds a \code{type} tag to each sequence, indicating if they
#'are either a \code{target} gene or a \code{template} sequence. User can optionally 
#'select what part of the sequence header will be used as internal id. 
#'
#' @param file A path to a FASTA file or a vector containing the file path.
#' @param input_type A string used to indicate if the sequence is a 'target' 
#'   or a 'template'.
#' @param id_field A string with the name of the id field. Defaults to 'id'.
#' @param separator A single character used to separate fields in the sequence
#'  headers. Defaults to '\\|'.
#'
#' @return A dataframe. 
#'
#' @examples
#' read_sequences("data/test.fasta", input_type = 'target')
#' read_sequences("data/test.fasta", id_field = 'name', separator = '-')
read_sequences <- function(file, input_type, id_field = NULL, separator = NULL) {
  raw <- as_data_frame(read_file(file))
  
  df <- raw %>%
    separate_rows(value, sep = '\n>') %>%
    separate(value, into = c('header', 'seq'), sep = '\n', extra = 'merge') %>%
    mutate(
      type = input_type,
      seq = str_replace_all(seq, '\n', ''),
      header = str_replace(header, '>', ''),
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
}


#' Append sequence input to global input dataframe
#' 
#' Internal use in the Shiny UI only.
#'
#' @param input 
#' @param global_value 
#'
#' @return
#' @export
#'
#' @examples
collect_input_shiny <- function(input, global_value) {
  new_input <- input
  global_value <- rbind(global_value, new_input)
}

#03/03/18, Wife's comment here:
#love you

get_overlaps <- function(seq_object, position=NULL, length=NULL) {
  #Initialize empty overlap object
  overlaps <- tibble(
    id = character(),
    type = character(),
    origin = character(),
    seq = character()
  )
  
  #Checking what type of object was passed to function:
  if(seq_object$type == 'template') {
    #Using default overlap length of 20 if not provided:
    overlap_length <- ifelse(is.null(length), 20, length)
    #Extracting overlaps
    overlap_left <- substr(seq_object$seq, (position-overlap_length), position)
    overlap_right <- substr(seq_object$seq, (position+1), (overlap_length+position+1))
    #Populating overlap object:
    overlaps <- overlaps %>% 
      add_row(id = 'left', type = 'overlap', origin = 'template', seq = overlap_left) %>%
      add_row(id = 'right', type = 'overlap', origin = 'template', seq = overlap_right)
  } else {
    #Using default overlap length of 20 if not provided:
    overlap_length <- ifelse(is.null(length), 15, length)
    #Extracting overlaps
    overlap_left <- substr(seq_object$seq, 1, overlap_length)
    overlap_right <- substr(seq_object$seq, seq_object$length-overlap_length, seq_object$length)
    #Populating overlap object:
    overlaps <- overlaps %>% 
      add_row(id = 'left', type = 'overlap', origin = 'target', seq = overlap_left) %>%
      add_row(id = 'right', type = 'overlap', origin = 'target', seq = overlap_right)
  }
}
