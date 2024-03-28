The "hazard scenario" runs for this project have been copied to the MDSS (tape) on NCI for project w85, using the following script
```
#!/bin/bash
#PBS -P w85
#PBS -q copyq
#PBS -l walltime=10:00:00
#PBS -lmem=4GB
#PBS -lncpus=1
#PBS -l wd
#PBS -l storage=scratch/w85+gdata/w85


# Copy the OUTPUTS folder for all scenarios to mdss
mdss mkdir tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
mdss put -r ptha18-NSW2023b-ID710.5-sealevel110cm tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
```

You should be able to get the files using an appropriate `mdss get` command.
