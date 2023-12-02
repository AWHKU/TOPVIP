evoEF=path_to/bin/EvoEF/EvoEF
wkdir=SpCas9_evoEF
### make structure
cd ${wkdir}
${evoEF} --command=RepairStructure --pdb=7ox9.pdb --num_of_runs=3
${evoEF} --command=BuildMutant --pdb=7ox9_Repair.pdb --mutant_file=SpCas9_mutlist.txt --num_of_runs=10

PDB=($(ls 7ox9_Repair_Model_*.pdb))
for i in ${PDB[@]}
do
${evoEF} --command=ComputeStability  --pdb=${i} > ${i/.pdb/.log}
done

touch SpCas9_log.txt

for i in ${PDB[@]}
do
grep -e "Total                 =              " -e " was read by EvoEF" ${i/.pdb/.log} | \
tr '\n' ' ' >> SpCas9_log.txt
done

sed -i "s/was read by EvoEF Total                 =              //g" SpCas9_log.txt

sed -i "s/pdb file /\n/g" SpCas9_log.txt