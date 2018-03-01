##############
#Primer design
##############

create_primers <- function(overlaps_target, overlaps_template) {
  overlaps <- bind_rows(overlaps_target, overlaps_template)
  
  primers <- overlaps %>%
    spread(origin, seq) %>%
    select(-type) %>%
    mutate(
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
    mutate(id = case_when(
      id == 'left' ~ gsub('left', 'forward', id),
      id == 'right' ~ gsub('right', 'reverse', id)
    )) %>%
    mutate(type = 'primer') %>%
    select(-target, -template)
}