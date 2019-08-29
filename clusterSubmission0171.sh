#!/bin/bash
# submit a job for each recording
for irec in {1..22}
do
	echo Submitting request $irec
	sbatch -p general -N 1 -n 36 -t 07-00:00:00 --mem=16g --wrap="matlab -nodesktop -nosplash -singleCompThread -r CSRTT_SpkPLV\($irec\) -logfile CSRTT_SpkPLV_Init$irec.out"
done