#!/bin/bash

cd ..

qsub test-tides/tides_lady_elliot_kermadectonga2.pbs
qsub test-tides/tides_lady_elliot_solomon2.pbs
qsub test-tides/tides_rosslyn_bay_kermadectonga2.pbs
qsub test-tides/tides_rosslyn_bay_solomon2.pbs
qsub test-tides/tides_south_trees_kermadectonga2.pbs
qsub test-tides/tides_south_trees_solomon2.pbs
qsub test-tides/tides_lady_elliot_large_kermadec.pbs
qsub test-tides/tides_lady_elliot_southamerica.pbs
qsub test-tides/tides_rosslyn_bay_large_kermadec.pbs
qsub test-tides/tides_rosslyn_bay_southamerica.pbs
qsub test-tides/tides_south_trees_large_kermadec.pbs
qsub test-tides/tides_south_trees_southamerica.pbs

qsub test-tides/tides_vary_kermadectonga2.pbs
qsub test-tides/tides_vary_solomon2.pbs
qsub test-tides/tides_vary_southamerica.pbs
qsub test-tides/tides_vary_large_kermadec.pbs
