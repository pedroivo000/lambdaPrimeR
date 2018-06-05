devtools::load_all()

#File inputs paths:
target_input <- system.file('extdata', 'GFP sequence.fasta', package = "lambdaPrimeR")
vector_input <- system.file('extdata', 'pUC19.fasta', package = "lambdaPrimeR")
imen_target_input <- system.file('extdata', 'imen_primer_targets.csv', package = "lambdaPrimeR")
imen_vector_input <- system.file('extdata', 'pUC19_imen.fasta', package = "lambdaPrimeR")
imen_primer_csv <- system.file('extdata', 'imen_primers_toeval.csv', package = "lambdaPrimeR")

imen_primers <- read_primers_csv(imen_primer_csv)
#Creating input objects:
# targets <- read_sequences_csv(target_csv, input_type = 'target')
inputs_imen <- read_inputs(imen_target_input, imen_vector_input, primer_csv)
primers <- get_annealing_regions(inputs_imen)
inputs_single <- read_inputs(target_input, vector_input)
# inputs_csv <- read_inputs(target_csv, vector_input)


#Creating primers
inputs_single <- create_primers(inputs_single)
primer_seqs <- get_sequences(inputs@primers)
  
#Thermodynamical calculations:
primer_scores_imen <- evaluate_primers(inputs_imen)
primer_scores <- evaluate_primers(inputs_single)
