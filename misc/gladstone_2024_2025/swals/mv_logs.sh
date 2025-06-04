#!/bin/bash

mv tohoku*.pbs.e* validation/log/
mv tohoku*.pbs.o* validation/log/
mv tohoku*.log* validation/log/

mv solomon*.pbs.e* validation/log/
mv solomon*.pbs.o* validation/log/
mv solomon*.log* validation/log/

mv model*.pbs.o* test/log/
mv model*.pbs.e* test/log/
mv model*.log* test/log/

mv extreme*.log* test-full/log/
mv extreme*.pbs.o* test-full/log/
mv extreme*.pbs.e* test-full/log/

mv small*.log* test-full/log/
mv small*.pbs.o* test-full/log/
mv small*.pbs.e* test-full/log/

mv tide_varying*.log* test-full/log/
mv tide_varying*.pbs.o* test-full/log/
mv tide_varying*.pbs.e* test-full/log/

mv gladstone*.log* run_ptha/log/
mv gladstone*.pbs.o* run_ptha/log/
mv gladstone*.pbs.e* run_ptha/log/
