#!/bin/bash

#PBS -l mem=3g
#PBS -l walltime=5:55:00
#PBS -l nodes=1:ppn=1
#PBS -m ae
#PBS -j oe

cd "${bowtielocation}"
if [[ "${fastqfile}" = *.gz ]]; then
	gunzip -c "${fastqfile}" | ./bowtie --sam ${genome} - "${samfile}"
else
	cat "${fastqfile}" | ./bowtie --sam ${genome} - "${samfile}"
fi

"${samtools}/samtools" view -bS -o "${bamfile}" "${samfile}"
if [[ ${keepsam} = 'no' ]]; then
	rm -Rf "${samfile}"
fi
