#!/bin/sh
ulimit -c unlimited
ulimit -s unlimited

./bin/dmbayes_chord step1_chord_mps3.ini > chains/log_chord_ph200_ident_subc.out


