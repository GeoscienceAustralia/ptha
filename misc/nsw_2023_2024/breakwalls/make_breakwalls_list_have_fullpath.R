x = readLines('breakwalls_file_list_relative_path.txt')
writeLines(normalizePath(x), con='breakwalls_file_list.txt')
