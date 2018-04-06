##############
# Primer class
##############

validate_primers <- function(object) {
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
    'forward_seq','character', 
    'reverse_seq','character',
    'forward_length','integer', 
    'reverse_length','integer', 
    'target_region_origin','character', 
    'vector_region_origin','character'
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


#' S4 class representing a group of lambda and theta primers.
#' 
#' The \code{Primers} class contains the information about the primers designed 
#' to amplify and insert the target sequence into the vector (the lambda primers),
#' and the small primer that will be used to amplify the other vector strand
#' (the theta primer). 
#'
#' @return The \code{Primers} constructor return an object of class \code{Primers},
#' an instance of data frame.
#'  
#' @export 
setClass("Primers", contains = "data.frame", validity = validate_primers)

Primers <- function(df, ...) new("Primers", ...)

# setGeneric("create_primers", function(template_list, ...) {
#   standardGeneric("create_primers")
# })

create_primers <- function(template_objects) {
  
}
#' Design primer pair from overlaps.
#'
#' This function uses the overlap information in the \code{overlaps} object
#' to design a primer pair by concatenating the overlap sequences in the
#' correct order (see Details).
#' 
#' In order to design the primers for the first PCR step of LambdaPCR, where we
#' add template-overlaping flanking regions to our target gene sequence, 
#' we have to use the following rules:
#'
#' \itemize{
#'   \item \strong{Forward primer:} left overlap + target left overlap, 
#'   5'-3' orientation
#'   \item \strong{Reverse primer:} complement of target right overlap + 
#'   complement of template right overlap,  5'-3' orientation
#' }
#'
#' @param overlaps A dataframe containing the overlap information for both input
#' types. The \code{overlaps} object is created by \code{get_overlaps}.
#'
#' 
#' @return A dataframe containing the forward and reverse primer sequences.
#' @export
create_primers <- function(template_list) {
  primers <- overlaps %>%
    #Add columns with overlap seqs:
    spread(origin, seq) %>%
    select(-type) %>%
    mutate(
      #We need to get the reverse complement of the right-oriented overlaps
      #and merge the overlaps sequence to get the final primer seqs:
      seq = case_when(
        id == 'left' ~ paste(tolower(template), 
                             target, 
                             sep = ''),
        id == 'right' ~ paste(tolower(reverse_complement(template)), 
                              reverse_complement(target), 
                              sep = ''
        )
      )
    ) %>%
    #changing id to primer orientation
    mutate(id = case_when(
      id == 'left' ~ gsub('left', 'forward', id),
      id == 'right' ~ gsub('right', 'reverse', id)
    )) %>%
    mutate(type = 'primer') %>%
    #extracting template and target annealing sequences
    mutate(temp = str_replace(seq, '([[:upper:]])', ',\\1')) %>%
    separate(temp, into = c('template_annealing_seq', 'target_annealing_seq'),
             sep = ',') %>%
    select(-target, -template)
  
}

