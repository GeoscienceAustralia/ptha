# YOU MUST CHANGE THE NAME OF THIS VARIABLE TO BE THE DIRECTORY CONTAINING THE 'PREDICT_TIDE' PROGRAM
.OTPS_directory = '/home/gareth/Code_Experiments/TIDAL_PREDICTION/TPXO9/OTPS'
.TIDAL_MODEL_CONTROL_FILE = 'DATA/Model_tpxo9.v5a'

if(!file.exists(.OTPS_directory)){
    msg = paste0('Cannot file OTPS_directory: ',
        .OTPS_directory, '\n Please edit the OTPS_directory_name_TPXO9.R script \n',
        'to ensure that the ".OTPS_directory" variable points to the OTPS folder in your tpxo9 install')
    stop(msg)
}
