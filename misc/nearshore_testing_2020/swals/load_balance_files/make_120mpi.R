# Scan the NSW default load-balance file with 8 nodes/32 mpi, and make a new one with 30nodes/120 MPI 

x = readLines('load_balance_default_NSW_8nodes_32mpi.txt')
nproc = sapply(strsplit(x, ' '), length)

nproc_120 = round(nproc*120/32)

# Manually convert to integers that have smaller factors, so there is greater
# flexibility in splitting the domain
nproc_120[nproc_120 == 15] = 16
nproc_120[nproc_120 == 45] = 48
nproc_120[nproc_120 == 75] = 72

inds_start_end = c(0, cumsum(nproc_120))

inds_start = inds_start_end[1:(length(inds_start_end)-1)] + 1
inds_end   = inds_start_end[2:(length(inds_start_end))]

l120 = rep(1:120, 100)

all_inds = mapply(f<-function(i1, i2) l120[i1:i2], inds_start, inds_end)

text_for_file = unlist(lapply(all_inds, f<-function(x) paste(x, collapse=" ")))

cat(text_for_file, file='load_balance_default_NSW_30nodes_120mpi.txt', sep="\n")
