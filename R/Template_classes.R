#' Validade a \code{Template} object.
#' 
#' Checks if a \code{Template} object or an object of any subclass is valid 
#' or not. 
#'
#' @param object The input data frame (returned by \code{read_sequences}) to be 
#' checked against a template data frame. 
#'
#' @return \code{TRUE} if object is a valid data frame, otherwise return a list
#' of detected errors during input validation. 
#' @keywords interal
#'
#' @examples
validate_template <- function(object) {
  #Vector that will contain all the error messages:
  errors <- NULL
  
  #First, checking if object is data frame:
  if(!is.data.frame(object)) {
    errors <- c(errors, "Input was not a data frame")
    return(FALSE)
  }
  
  #Checking if required columns are present with the correct format:
  required_columns <- tribble(
    ~cols, ~req_class,
    'path','character', 
    'header','character',
    'name','character', 
    'seq','character', 
    'type','character', 
    'length','integer'
  )
  
  #Creating data frame with object column types:
  column_types <- as.data.frame(lapply(object, class)) %>%
    mutate_all(funs(as.character(.))) %>%
    gather(cols, obj_class)
  
  #Checking column types against required columns:
  col_class_test <- suppressMessages(left_join(required_columns, column_types)) %>%
    mutate(test = ifelse(req_class == obj_class, TRUE, FALSE))
  
  #Error messages:
  discrepancies <- col_class_test %>%
    do(
      msg = case_when(
        is.na(.$test) ~ paste("Missing required column: ", .$cols, ";", sep = ''),
        !(.$test) ~ paste("Column ", .$cols, " type is incorrect: ",
                          "required type is ", .$req_class, ", current type is ",
                          .$obj_class, ";", sep = ''),
        (.$test) ~ NA_character_
      )
    ) %>%
    unnest() %>%
    filter(!is.na(msg))
  
  errors <- c(errors, discrepancies$msg)
  
  if(length(errors) != 0) {
    return(errors)
  } else {
    return(TRUE)  
  }
}

#' S4 class representing a DNA sequence template.
#' 
#' This is the master class for the \code{lambdaPrimeR} package. It represents
#' a input DNA sequence that can either be the Target or the Vector sequence.
#' Both the target and vector sequences will be used as templates for the design 
#' of the lambda PCR primers.
#' 
#' This master-class is internal, as the two subclasses \code{Target} and
#' \code{Vector} will be the ones used to design the primers in this package.
#'
#' @return The \code{Template} constructor returns an object of class \code{Template},
#' an instance of a data frame.
#' @export
#'
setClass("Template", contains = c("data.frame"))

#'Read sequences from file.
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
#' read_sequences("data/test.fasta", input_type = 'vector', header_structure = header, separator = '-')
read_sequences <- function(file, input_type = c('target', 'vector'),
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
      length = nchar(seq)
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
    select(name, path, type, header, seq, length)
  
  #Checking what type of input was provided and creating appropriate class obj.:
  object <- switch (input_type,
    'target' = Target(df),
    'vector' = Vector(df)
  )

  return(object)
}

setClass("Target", contains = c("Template"), validity = validate_template)
Target <- function(df, ...) new("Target", df, ...)

setClass("Vector", contains = c("Template"), validity = validate_template)
Vector <- function(df, ...) new("Vector", df, ...)