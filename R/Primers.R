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
    'id', 'integer',
    'forward_primer_seq','character', 
    'reverse_primer_seq','character',
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

Primers <- function(df, ...) new("Primers", df, ...)

setGeneric("get_sequences", function(object, ...) standardGeneric("get_sequences"))

setMethod("get_sequences", 
  signature = "Primers", 
  function(object) {
    # primer_seqs <- asS3(object)
    # primer_seqs <- primer_seqs %>%
    #   sperate()
      
    
    df <- tribble(
      ~orientation, ~primer_region, ~seq,
      'forward', 'complete', object$forward_primer_seq, 
      'forward', 'target', object$forward_target_anneal_seq,
      'forward', 'vector', object$forward_vector_anneal_seq, 
      'reverse', 'complete', object$reverse_primer_seq, 
      'reverse', 'target', object$reverse_target_anneal_seq,
      'reverse', 'vector', object$reverse_vector_anneal_seq
    )
  
    return(df)  
})