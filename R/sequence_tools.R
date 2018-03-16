################
#Sequence tools
################

complement <- function(sequence) {
  seq <- toupper(sequence)
  complement <- chartr('ATCG', 'TAGC', seq)
}

reverse_complement <- function(sequence) {
  complement <- complement(sequence)
  revcomp <- stringi::stri_reverse(complement)
}

simulate_pcr <- function(sequence_object, primers) {
  #Get original sequence from sequence input:
  orginal_seq <- sequence_object
  forward_primer <- filter(primers, id == 'forward')
  reverse_primer <- filter(primers, id == 'reverse')
  left_flank <- forward_primer$template_annealing_seq
  right_flank <- reverse_complement(reverse_primer$template_annealing_seq)
  
  extended_seq <- paste(left_flank, orginal_seq, tolower(right_flank), sep = '')
}