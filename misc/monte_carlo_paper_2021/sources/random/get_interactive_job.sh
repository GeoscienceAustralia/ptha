# On NCI we run this to start up an interactive job on project w85, with access to key disks required to create the scenarios
qsub -I -X -P w85 -q express -l ncpus=48 -l mem=192gb -l wd -l jobfs=1gb -l walltime=01:00:00 -lstorage=scratch/w85+gdata/w85+gdata/fj6+scratch/n74+gdata/n74
