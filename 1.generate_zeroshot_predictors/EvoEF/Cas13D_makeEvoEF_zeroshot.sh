evoEF=path_to/EvoEF/EvoEF
wkdir=Cas13d_evoEF
### make structure
cd ${wkdir}
${evoEF} --command=RepairStructure --pdb=RxfCas13d_Itasser.pdb --num_of_runs=3
${evoEF} --command=BuildMutant --pdb=RxfCas13d_Itasser_Repair.pdb --mutant_file=Cas13d_mutlist.txt --num_of_runs=10

PDB=($(ls RxfCas13d_Itasser_Repair_Model_*.pdb))
for i in ${PDB[@]}
do
${evoEF} --command=ComputeStability  --pdb=${i} > ${i/.pdb/.log}
done

touch Cas13_log.txt

for i in ${PDB[@]}
do
grep -e "Total                 =              " -e " was read by EvoEF" ${i/.pdb/.log} | \
tr '\n' ' ' >> Cas13_log.txt
done

sed -i "s/was read by EvoEF Total                 =              //g" Cas13_log.txt

sed -i "s/pdb file /\n/g" Cas13_log.txt