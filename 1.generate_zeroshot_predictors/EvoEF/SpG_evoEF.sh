evoEF=path_to/EvoEF/EvoEF
wkdir=SpCas9_evoEF
cd ${wkdir}
${evoEF} --command=BuildMutant --pdb=7ox9_Repair.pdb --mutant_file=SpG_evoEF_mutlist.txt --num_of_runs=10
PDB=($(ls 7ox9_Repair_Model_*.pdb))
for i in ${PDB[@]}
do
${evoEF} --command=ComputeStability  --pdb=${i} > ${i/.pdb/.log}
done

touch SpG_log.txt

for i in ${PDB[@]}
do
grep -e "Total                 =              " -e " was read by EvoEF" ${i/.pdb/.log} | \
tr '\n' ' ' >> SpG_log.txt
done

sed -i "s/was read by EvoEF Total                 =              //g" SpG_log.txt

sed -i "s/pdb file /\n/g" SpG_log.txt