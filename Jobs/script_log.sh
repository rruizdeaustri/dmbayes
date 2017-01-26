cd ..
clusterlauncher -N "test" -n 96 --fast -q eternity Jobs/sub_job.sh
#cd ..
#clusterlauncher -N "test" -n 95 --fast -q eternity ./bin/dmbayes_chord step1_chord_mps_clust.ini \> chains/log_chord.out
#clusterlauncher -v -N "mn_1" -n 64 -q long ./superbayes SampleIniFile_zg.ini \> ../chains/m1p_mup.out
