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
#mdss put -r ptha18-NSW2023b-ID710.5-sealevel110cm tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r ptha18-NSW2023-ID1315.5-sealevel110cm tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r ptha18-NSW2023-ID4186.3-sealevel110cm tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/

#mdss put -r run_chile1960_FujiSatake2013_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_chile1960_HoEtAl_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_Chile2010_Lorito11_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_Chile2015_Williamson_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_extreme_source_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_Kermadec2021_Romano_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_kt43731_12h_final_NNL4_1arcminoffshore-full-ambient_sea_level_1.1 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_kt43731_12h_final_NNL4_CONVERGENCE-full-ambient_sea_level_1.1 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_NewHebrides2021_GusmanWithKajiura_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_Puysegur2009_PTHA18HS1567_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_Puysegur2009_PTHA18HS1788_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_small_source_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_solomon2007_11244_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_solomon2007_1_19_22m_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_tohoku2011_Yamazaki_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r runs_before_Nov28_2023 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
#mdss put -r run_Chile2014_PTHA18_HS80427_1arcminoffshore-full-ambient_sea_level_0.0 tsunami/MODELS/inundation/NSW_tsunami_inundation_with_NSWSES/stage1/full_coast/swals/OUTPUTS/
