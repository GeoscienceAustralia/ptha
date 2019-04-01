
# Report the fortran test. Note this is not an analytical result,
# we are just checking for a change from the standard value.
outlog = readLines('outfile.log')
if(tail(outlog)[4] == ' PASS'){
    print('PASS')
}else{
    print('FAIL')
}
