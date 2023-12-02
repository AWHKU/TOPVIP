evoEF=path_to/EvoEF/EvoEF
wkdir=SaCas9_evoEF
### make structure
cd ${wkdir}
#${evoEF} --command=RepairStructure --pdb=5cZZ_KKH.pdb --num_of_runs=3
${evoEF} --command=BuildMutant --pdb=5cZZ_KKH_Repair.pdb --mutant_file=SAV_evoef_mutlist.txt --num_of_runs=10

PDB=($(ls 5cZZ_KKH_Repair_Model_*.pdb))
for i in ${PDB[@]}
do
${evoEF} --command=ComputeStability  --pdb=${i} > ${i/.pdb/.log}
done

touch SAV_evoef_log.txt

for i in ${PDB[@]}
do
grep -e "Total                 =              " -e " was read by EvoEF" ${i/.pdb/.log} | \
tr '\n' ' ' >> SAV_evoef_log.txt
done

sed -i "s/was read by EvoEF Total                 =              //g" SAV_evoef_log.txt

sed -i "s/pdb file /\n/g" SAV_evoef_log.txt

