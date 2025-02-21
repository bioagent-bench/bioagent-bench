mamba create -n rnaseq
mamba activate rnaseq

mkdir data
prefetch SRR1278968 SRR1278969 SRR1278970 SRR1278971 SRR1278972 SRR1278973 -O data/ 