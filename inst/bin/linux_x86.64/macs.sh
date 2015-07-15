#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -m ae
#PBS -j oe
#PBS -o ${templocation}/${id}.log

cd "${templocation}"
"${macslocation}/macs14" -t "${treatmentfile}" -c "${controlfile}" -n "${id}" -g "${assembly}" --bw ${fragmentlength} --bdg --single-profile

if [[ "${pospeakfile}" =~ \.txt\.gz$ ]]; then
	gzip -9 -c "${id}_peaks.xls" > "${pospeakfile}"
elif [[ "${pospeakfile}" =~ \.txt$ ]]; then
	mv "${id}_peaks.xls" "${pospeakfile}"
elif [[ "${pospeakfile}" =~ \.bed\.gz$ ]]; then
	gzip -9 -c "${id}_peaks.bed" > "${pospeakfile}"
elif [[ "${pospeakfile}" =~ \.bed$ ]]; then
	mv "${id}_peaks.bed" "${pospeakfile}"
fi

if [[ "${negpeakfile}" =~ \.txt\.gz$ ]]; then
	gzip -9 -c "${id}_negative_peaks.xls" > "${negpeakfile}"
elif [[ "${negpeakfile}" =~ \.txt$ ]]; then
	mv "${id}_negative_peaks.xls" "${negpeakfile}"
fi

if [[ "${trackfile}" =~ \.gz$ ]]; then
	mv "${id}_MACS_bedGraph/treat/${id}_treat_afterfiting_all.bdg.gz" "${trackfile}"
fi

if [ -n "${macsmodelfile}" ]; then
#	cat '' > "${macsmodelfile}"
	for symbol in p m s; do
		cat "${id}_model.r" | grep "^${symbol} <- c(\([0-9.,]*\))$" | sed -e "s/^${symbol} <- c(\([0-9.,]*\))/,\1/g" -e "s/,/;${symbol},/g" | tr ';' '\n' >> "${macsmodelfile}"
	done
fi

#rm -Rf "${templocation}"
