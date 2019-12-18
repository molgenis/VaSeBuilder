#!/bin/bash

ml argon2

echo -e "#DonorSample\tSalt\tHash\tEncodings" > "${2}"

while read -r line; do
	if [[ "${line}" == \#* ]]; then
		continue
	fi
	sample_id=$(cut -f 2 <<< "${line}")
	sample_salt=$(tr -dc '[:alnum:]' < /dev/urandom | head -c 12)
#	sample_hash=$(argon2 "${sample_salt}" -r <<< "${sample_id}")
	sample_hash=$(awk -v IFS=" " -v OFS="\t" '{print $11,$13}' <<< $(argon2 "${sample_salt}" -l 6 <<< "${sample_id}"))
	echo -e "${sample_id}\t${sample_salt}\t${sample_hash}" >> "${2}"
done < "${1}"

./python3 Locker.py "${1}" "${2}"
