template_file = '../multidomain_design/domains_2023_09_07_0.5_0.166666666666667_0.0333333333333333/load_balance_default.txt'
template = readLines(template_file)

new_template = template

# Give 96 coarray images
new_template[1] = paste0(seq(1, 96), collapse=" ")

last_image_allocated = length(new_template) - 1

# From experience, these domains should be split to enable a load balance
spi = c(146, 233, 286, 250, 164, 273)
for(i in 1:length(spi)){
    ind = spi[i]
    new_template[ind] = paste0(new_template[ind], ' ', last_image_allocated + i)
}

output_file = './load_balance_files/load_balance_partition_96MPI_4NL_CONVERGENCE_DEFAULT.txt'

writeLines(new_template, output_file)
