
context('test_write_table_to_ncdf4')

test_that('test_write_table_to_ncdf4', {

    test_table = data.frame(x = c(1,2,3), y=c('a', 'b', 'csdf'), z=c(1.1, 1.2, 1.3),
        stringsAsFactors=FALSE)

    # Rely on internal logical checks to catch errors
    errflag = try({
        write_table_to_netcdf(test_table, 
            file='test.nc', 
            units=c('m', '', 'kg'), 
            long_names=c('head count', 'mychar', 'asdfasdfa'), 
            var_prec=c('double', 'char', 'double'))
        })

    expect_true(class(errflag) != 'try-error')

    test_table_new = read_table_from_netcdf('test.nc')
    expect_true(isTRUE(all.equal(test_table, test_table_new)))

    # Clean up
    if(class(errflag) != 'try-error') unlink('test.nc')
})
