#!/usr/bin/bash

scriptdir=$(dirname $(realpath $0))

function clean_tmp_files()
{
    (
        cd $(dirname $1)
        rm -f CELL CTRL ELECTRONS KPOINTS POS SYSTEM SPECIES
    )
}

clean_tmp_files $1

cat $1 | sed -e 's/ATOMIC_SPECIES/\&ATOMIC_SPECIES/g' \
             -e 's/ATOMIC_POSITIONS/\&ATOMIC_POSITIONS\n#/g' \
             -e 's/K_POINTS/\&K_POINTS\n#/g' \
             -e 's/CELL_PARAMETERS/\&CELL_PARAMETERS\n#/g' | sed 's/\&/\&\n/' | sed '/^\/$/d' | tee out | ${scriptdir}/split.awk

python ${scriptdir}/parse.py $(dirname $1)
# clean_tmp_files $1
