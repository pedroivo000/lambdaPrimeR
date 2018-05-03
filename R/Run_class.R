#############
# Run class #
#############

#' @include Template_classes.R Primers.R
NULL

validate_run <- function(object) {
  errors <- NULL
  
  id <- object@id
  if (!is.character(id)) {
    msg <- paste( "Run id needs to be a character string. Type provided was ", 
                  class(id), ".", sep='')
    errors <- c(errors, msg)
  }
  
  #if no errors, create instance of class:
  if(length(errors) != 0) {
    return(errors)
  } else {
    return(TRUE)
  }
}

setClass("Run", 
  slots = c(
    target = "Target",
    vector = "Vector", 
    primers = "Primers", 
    # options = "Options",
    id = "character"
  ), 
  validity = validate_run
)

Run <- function(target_object, vector_object, primer_object = NULL, id=NULL, ...) {
  if(is.null(primer_object)) {
    #The Primer slot will only be populated when a group of primers is designed 
    #for the input templates. When we create a new Run object, we need to create
    #an empty Primer object to act as a place holder:
    empty_primer_df <- tibble(
      forward_primer_seq = character(),
      reverse_primer_seq = character(),
      forward_length = integer(),  
      reverse_length = integer(),  
      target_region_origin = character(),  
      vector_region_origin = character()
    )
    #creating empty Primer object:
    primer_object <- Primers(empty_primer_df)
  } 
  
  #id has to be a unique string identifying the current primer design run:
  if(is.null(id)) {
    #concatenate target and vector:
    target_name <- target_object$name
    vector_name <- vector_object$name
    id <- paste(target_name, vector_name, sep = '-')
    #remove whitespaces:
    id <- gsub('\\s+', '_', id)
  } 
  
  new("Run", target = target_object, vector = vector_object, 
      primers = primer_object, id = id)
}

read_inputs <-  function(target, vector, id = NULL, header_structure = NULL, 
                         name_field = NULL, sep = NULL) {
  
  #Creating Template objects from file inputs:
  Target <- read_sequences(target, input_type = 'target', header_structure, 
                           name_field, sep)
  Vector <- read_sequences(vector, input_type = 'vector', header_structure,
                           name_field, sep)
  
  #Creating a Run object from the Template inputs:
  Run <- Run(Target, Vector, id = id)
}

setGeneric("primers", function(x) standardGeneric("primers"))

setMethod("primers", signature = "Run", function(x) x@primers)

setGeneric("get_inputs", function(x) standardGeneric("get_inputs"))

setMethod("get_inputs", signature = "Run", function(x) {
  target <- x@target
  vector <- x@vector
  inputs <- list(target = target, vector = vector)
}) 
