pandoc -f commonmark -V papersize=a4 -V geometry:margin=1in -V geometry:landscape $1  -o $1.pdf 
