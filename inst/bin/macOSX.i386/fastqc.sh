#!/bin/bash

#PBS -l mem=400mb
#PBS -l walltime=0:35:00
#PBS -l nodes=1:ppn=1
#PBS -m ae
#PBS -j oe

cd "${fastqclocation}"
perl fastqc -o "${reportdir}" --noextract -q "${inputfile}"

fname="${inputfile##*/}"
fname="${fname%%.*}"
if [ "${fname}" != "${sample}" ]; then
	cd "${reportdir}"
	unzip "${fname}_fastqc.zip"
	mv "${fname}_fastqc" "${sample}_fastqc"
	zip -9rq "${sample}_fastqc" "${sample}_fastqc"
fi
