# Make pdfs
STARTDIR=$(pwd)
for i in ./README.md ./*/README.md ./nthmp/*/README.md; do echo $i; cd $(dirname $i); pandoc -V papersize=a4 -V geometry:margin=1in README.md -o README.pdf -f markdown-implicit_figures ; cd $STARTDIR; done


# Combine to a single pdf
gs -dNOPAUSE -sDEVICE=pdfwrite -SOUTPUTFILE=all_test_problems.pdf -dBATCH ./nthmp/README.pdf ./nthmp/*/README.pdf README.pdf ./[a-m]*/README.pdf ./nesting_plane_wave/README.pdf ./[o-z]*/README.pdf
