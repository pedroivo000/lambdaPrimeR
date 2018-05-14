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
    'id', 'integer',
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

#' S4 class representing a target gene sequence. 
#' 
#' The class \code{Target} is a subclass of \code{Template}. It represents the 
#' sequence of a target gene to be inserted in a vector using Lambda-PCR. 
#' 
#' @return The \code{Target} constructor returns an object of class \code{Target},
#' a subclass of \code{Template}.
#' @export
setClass("Target", contains = c("Template"), validity = validate_template)
Target <- function(df, ...) new("Target", df, ...)

#' S4 class representing a vector sequence. 
#' 
#' The class \code{Vector} is a subclass of \code{Template}. It represents the 
#' sequence of a destination vector that the target gene will be inserted into
#' using Lambda-PCR. 
#' 
#' @return The \code{Vector} constructor returns an object of class \code{Vector},
#' a subclass of \code{Template}.
#' @export
setClass("Vector", contains = c("Template"), validity = validate_template)
Vector <- function(df, ...) new("Vector", df, ...)