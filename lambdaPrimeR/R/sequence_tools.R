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
