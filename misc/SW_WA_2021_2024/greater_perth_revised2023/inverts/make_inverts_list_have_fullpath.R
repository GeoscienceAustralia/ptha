x = readLines('inverts_file_list_relative_path.txt')
writeLines(normalizePath(x), con='inverts_file_list.txt')
