# Check the tests embedded in the fortran log worked.
outlog = tail( readLines('outfile.log'), n=10)

PASS_LINES = c(2, 6, 10)
for(i in 1:length(PASS_LINES)){
    if(outlog[PASS_LINES[i]] == ' PASS'){
        print('PASS')
    }else{
        print('FAIL')
    }
}

