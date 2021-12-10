# Report the fortran test. Note this is not an analytical result,
# just a regression test. In future we should get a high-resolution
# reference model result to compare against.
outlog = readLines('outfile.log')
if(tail(outlog)[3] == ' PASS'){
    print('PASS')
}else{
    print('FAIL')
}

