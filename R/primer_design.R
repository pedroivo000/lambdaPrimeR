##############
#Primer design
##############

#' Extract target or vector-annealing regions from Template sequences. 
#'
#' @param template 
#' @param position 
#' @param min_length 
#' @param max_length 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
setGeneric("get_annealing_regions", 
 function(template, position = NULL, min_length = 15, max_length = 30, 
          multiple = TRUE, ...) {
   standardGeneric("get_annealing_regions")
})

setMethod("get_annealing_regions", 
  signature = "Target",
  function(template, min_length = 15, max_length = 30, mutiple = TRUE) {
    
    #Coercing Target object to S3 data frame:
    targets <- asS3(template)
    
    #Expanding target data frame with annealing region lengths:
    targets <- targets %>%
      mutate(target_annealing_length = list(min_length:max_length)) %>%
      unnest()
    
    #Saving sequence slot from template class in a variable:
    seq <- template$seq
    seq_length <- template$length
    
    #Extracting the annealing regions for each length value:
    #When extracting the annealing regions from the target sequence, we need to 
    #get the 5' and 3' ends of the gene sequence:
    #Extracting overlaps
    targets <- targets %>%
      group_by(target_annealing_length) %>%
      mutate(
        target_anneal_left_seq  = substr(seq, 1, target_annealing_length),
        target_anneal_right_seq = substr(seq, length-target_annealing_length, length),
        target_anneal_left_beg = 1,
        target_anneal_left_end = target_annealing_length,
        target_anneal_right_beg = length - target_annealing_length,
        target_anneal_right_end = length 
      )
    
    #Creating target object:
    target_object <- Target(targets)
    return(target_object)
  }
)

setMethod("get_annealing_regions", 
          signature = "Vector",
          function(template, position = NULL, min_length = 15, max_length = 30,
                   multi = TRUE) {
            #Coercing Target object to S3 data frame:
            vector <- asS3(template)
            
            #Expanding vector data frame with annealing region lengths:
            vector <- vector %>%
              mutate(vector_annealing_length = list(min_length:max_length)) %>%
              unnest()

            #Extracting the annealing regions for each length value:
            #When extracting the annealing regions from the vector sequence, we need to 
            #get the 5' and 3' ends of the gene sequence:
            #Extracting overlaps
            vector <- vector %>%
              group_by(vector_annealing_length) %>%
              mutate(
                vector_anneal_left_seq = substr(seq, position-vector_annealing_length, position),
                vector_anneal_right_seq = substr(seq, (position+1), (position+vector_annealing_length+1)),
                vector_anneal_left_beg = position - vector_annealing_length,
                vector_anneal_left_end = position,
                vector_anneal_right_beg = position + 1,
                vector_anneal_right_end = position + vector_annealing_length + 1
              )
            
            #Creating target object:
            vector_object <- Vector(vector)
            return(vector_object)
          }
)

setMethod("get_annealing_regions", 
  signature = "Run", 
  function(template, position = NULL, min_length = NULL, max_length = NULL) {
    #Extracting slots:
    Target <- template@target
    Vector <- template@vector
    Primers <- template@primers
    
    #If a primer object was created from a .csv file with the primer sequences, 
    #we need to update the Target and Vector objects with the coordinantes of 
    #primer annealing regions. Second, if the primer object constains multiple
    #primers, we need to match each primer with its respective target and vector
    #names:
    
    #Coercing Target, Vector and Primer objects to S3 data frame
    targets <- asS3(Target) %>%
      rename(target_region_origin = name) %>%
      select(-id)
    vector <- asS3(Vector) %>%
      select(name, seq) %>%
      rename(
        vector_region_origin = name, 
        vector_seq = seq
      ) 
    primers <- asS3(Primers)
    
    #Locating primer annealing region in the Target and Vector sequences:
    annealing_regions <- left_join(primers, vector) %>%
      left_join(targets) %>%
      group_by(id) %>%
      do(
        forward_target_anneal = as.tibble(str_locate(.$seq, toupper(.$forward_target_anneal_seq))),
        reverse_target_anneal = as.tibble(str_locate(.$seq, reverse_complement(.$reverse_target_anneal_seq))),
        forward_vector_anneal = as.tibble(str_locate(.$vector_seq, toupper(.$forward_vector_anneal_seq))),
        reverse_vector_anneal = as.tibble(str_locate(.$vector_seq, reverse_complement(toupper(.$reverse_vector_anneal_seq))))
      ) %>%
      unnest(.sep = '_')
    names(annealing_regions) <- gsub('start', 'beg', names(annealing_regions))
    
    #Merging coordinates to primer table:
    primers <- primers %>%
      left_join(annealing_regions)
    
    #Upadating Primer object:
    Primers <- Primers(primers)
    # return(Primers)
    # #We also need to add the coordinates to the Target and Vector objects:
    # Target$target_anneal_min_length <- nchar(primers$forward_target_anneal_seq)
    # Target$target_anneal_max_length <- nchar(primers$forward_target_anneal_seq)
    # Target$target_anneal_left_beg <- primers$forward_target_anneal_beg
    # Target$target_anneal_left_end <- primers$forward_target_anneal_end
    # Target$target_anneal_right_beg <- primers$reverse_target_anneal_beg
    # Target$target_anneal_right_end <- primers$reverse_target_anneal_end
    # Target$target_anneal_left_seq <- primers$forward_target_anneal_seq
    # Target$target_anneal_right_seq <- reverse_complement(primers$reverse_target_anneal_seq)
    # 
    # #Adding annealing region info to vector object:
    # Vector$vector_anneal_min_length <- nchar(primers$forward_vector_anneal_seq)
    # Vector$vector_anneal_max_length <- nchar(primers$forward_vector_anneal_seq)
    # Vector$vector_anneal_left_beg <- primers$forward_vector_anneal_beg
    # Vector$vector_anneal_left_end <- primers$forward_vector_anneal_end
    # Vector$vector_anneal_right_beg <- primers$reverse_vector_anneal_beg
    # Vector$vector_anneal_right_end <- primers$reverse_vector_anneal_end
    # Vector$vector_anneal_left_seq <- primers$forward_vector_anneal_seq
    # Vector$vector_anneal_right_seq <- reverse_complement(primers$reverse_vector_anneal_seq)
    # 
    #Update Run object with annealing regions:
    run_object <- Run(Target, Vector, Primers)
    return(run_object)
  }
  
)

#' Create primers from target and vector DNA sequences.
#'
#' @param run_object 
#'
#' @return
#' @export
#'
#' @examples
create_primers <- function(run_object) {
  #Extracting template slots from Run object:
  Target <- run_object@target
  Vector <- run_object@vector
  # Primers <- run_object@primers
  
  #Extracting annealing region sequences from Target and Vector objects:
  Target <- get_annealing_regions(Target, min_length = 15, max_length = 30)
  Vector <- get_annealing_regions(Vector, position = 1500, 
                                    min_length = 15, max_length = 30)
  
  target_left <- Target$target_anneal_left_seq
  target_right <- Target$target_anneal_right_seq
  vector_left <- Vector$vector_anneal_left_seq
  vector_right <- Vector$vector_anneal_right_seq
  
  #Assembling primer sequences:
  primers <- tibble(
    forward_primer_seq = paste(tolower(vector_left), target_left, sep = ''),
    reverse_primer_seq = paste(tolower(reverse_complement(vector_right)), 
                        reverse_complement(target_right), sep = ''),
    forward_target_anneal_beg = Target$target_anneal_left_beg,
    forward_target_anneal_end = Target$target_anneal_left_end,
    forward_vector_anneal_beg = Vector$vector_anneal_left_beg,
    forward_vector_anneal_end = Vector$vector_anneal_left_end,
    reverse_target_anneal_beg = Target$target_anneal_right_beg,
    reverse_target_anneal_end = Target$target_anneal_right_end,
    reverse_vector_anneal_beg = Vector$vector_anneal_right_beg,
    reverse_vector_anneal_end = Vector$vector_anneal_right_end,
    forward_target_anneal_seq = target_left,
    forward_vector_anneal_seq = vector_left,
    reverse_target_anneal_seq = reverse_complement(target_right),
    reverse_vector_anneal_seq = reverse_complement(vector_right),
    forward_length = nchar(forward_primer_seq),  
    reverse_length = nchar(reverse_primer_seq),  
    target_region_origin = Target$name,  
    vector_region_origin = Vector$name
  ) 
  
  primers <- primers %>%
    mutate(id = row_number())
  
  #Creating Primers object:
  Primers <- Primers(primers)
  
  #Updating Run object:
  run_object@target <- Target
  run_object@vector <- Vector
  run_object@primers <- Primers
  
  return(run_object)
}

# add_overlaps <- function(overlap_df) {
#   cols <- colnames(overlap_df)
# }

#' Extract overlaps from input sequences.
#'
#' @param sequence_inputs A dataframe containing the target and template input
#' sequences, generated by \code{read_sequences}.
#' @param position Numeric value indicating the insertion coordinate
#' @param length Length of the overlaping region to be extract from the input
#' sequences. Defaults to different values depending on the value of \code{type}:
#' 
#' \itemize{
#'   \item \code{type} = 'template' => \code{length} = 20
#'   \item \code{type} = 'target' => \code{length} = 15
#' }
#'
#' @return A dataframe containing the overlap sequence and orientation for each
#' input type. 
#' 
#' @export
#' @keywords internal
#'
get_overlaps <- function(sequence_inputs, position=NULL, length=NULL) {
  #Initialize empty overlap object
  overlaps <- tibble(
    id = character(),
    type = character(),
    origin = character(),
    seq = character()
  )
  
  #Function to find overlaps on each input type:
  find_overlaps <- function(data) {
    #Checking what type of input was passed:
    if(data$type == 'template') {
      #Using default overlap length of 20 if not provided:
      overlap_length <- ifelse(is.null(length), 20, length)
      #Extracting overlaps
      overlap_left <- substr(data$seq, (position-overlap_length), position)
      overlap_right <- substr(data$seq, (position+1), (overlap_length+position+1))
      #Populating overlap object:
      overlaps <- overlaps %>%
        add_row(id = 'left', type = 'overlap', origin = 'template', seq = overlap_left) %>%
        add_row(id = 'right', type = 'overlap', origin = 'template', seq = overlap_right)
    } else {
      #Using default overlap length of 15 if not provided:
      overlap_length <- ifelse(is.null(length), 15, length)
      #Extracting overlaps
      overlap_left <- substr(data$seq, 1, overlap_length)
      overlap_right <- substr(data$seq, data$length-overlap_length, data$length)
      #Populating overlap object:
      overlaps <- overlaps %>%
        add_row(id = 'left', type = 'overlap', origin = 'target', seq = overlap_left) %>%
        add_row(id = 'right', type = 'overlap', origin = 'target', seq = overlap_right)
    }
  }
  
  #Extracting overlaps
  sequence_inputs %>%
    group_by(type) %>%
    do(find_overlaps(.)) %>%
    ungroup()
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
create_primer_pair <- function(overlaps) {
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