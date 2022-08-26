# Make pdfs
STARTDIR=$(pwd)
for i in */README.md nthmp/*/README.md; do echo $i; cd $(dirname $i); pandoc -V papersize=a4 -V geometry:margin=1in README.md -o README.pdf; cd $STARTDIR; done


# Combine to a single pdf
gs -dNOPAUSE -sDEVICE=pdfwrite -SOUTPUTFILE=all_test_problems.pdf -dBATCH */README.pdf nthmp/*/README.pdf
