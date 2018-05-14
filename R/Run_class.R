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
  #id has to be a unique string identifying the current primer design run:
  if(is.null(id)) {
    #concatenate target and vector:
    target_name <- target_object$name
    vector_name <- vector_object$name
    id <- paste(target_name, vector_name, sep = '-')
    #remove whitespaces:
    id <- gsub('\\s+', '_', id)
  }
  
  #Check if primer object was passed:
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
      vector_region_origin = character(), 
      id = integer()
    )
    #creating empty Primer object:
    primer_object <- Primers(empty_primer_df)
  } 

  #Create new Run object:
  new("Run", target = target_object, vector = vector_object, 
      primers = primer_object, id = id)
  
}

read_inputs <-  function(target, vector, primers = NULL, 
                         target_vector_pairs = NULL, id = NULL, 
                         header_structure = NULL, 
                         name_field = NULL, sep = NULL) {
  
  #Creating Template objects from file inputs:
  Target <- read_sequences(target, input_type = 'target', header_structure, 
                           name_field, sep)
  Vector <- read_sequences(vector, input_type = 'vector', header_structure,
                           name_field, sep)
  
  #Checking if number of targets/vectors in the Template objects is not 1. 
  #We need to return a list of Run objects corresponding to each target-vector
  #pair in the input sequences:
  if(length(Target$id) != 1 || length(Vector$id) != 1) {
    #We need to check if target_vector_pair table was provided:
    if(is.null(target_vector_pairs)) {
      stop("Target-vector pair table not provided. A table with each Target-Vector 
           pair is necessary if multiple inputs are provided\n")
    } else {
      #Coercing Target and Vector objects into S3 data frames and splitting them 
      #into each row:
      targets <- asS3(inputs_imen@target) %>%
        group_by(name) %>%
        do(target_object = Target(.)) %>%
        rename(target = name)
      vectors <- asS3(inputs_imen@vector) %>%
        group_by(name) %>%
        do(vector_object = Vector(.)) %>%
        rename(vector = name)
      
      #Merging with pair table:
      run_list <- target_vector_pairs %>%
        left_join(targets) %>%
        left_join(vectors) %>%
        rowwise() %>%
        do(run_objects = Run(.$target_object, .$vector_object)) %>%
        pull(run_objects)
      
      names(run_list) <- lapply(run_list, function(x)x@id)
      return(run_list)
        
    }
    
  }
  
  #Check if primer input file was provided:
  if(!is.null(primers)) {
    #Creating Primers object from primer sequence table:
    Primers <- read_primers_csv(primers)
    
    #Creating a Run object with all inputs:
    Run <- Run(Target, Vector, Primers, id = id)
    
    #Get target and vector annealing regions from primer sequences:
    Run <- get_annealing_regions(Run)
  } else {
    #Creating a Run object with just Template inputs:
    Run <- Run(Target, Vector, id = id)  
  }
  
  
}

setGeneric("primers", function(x) standardGeneric("primers"))

setMethod("primers", signature = "Run", function(x) x@primers)

setGeneric("get_inputs", function(x) standardGeneric("get_inputs"))

setMethod("get_inputs", signature = "Run", function(x) {
  target <- x@target
  vector <- x@vector
  inputs <- list(target = target, vector = vector)
}) 
