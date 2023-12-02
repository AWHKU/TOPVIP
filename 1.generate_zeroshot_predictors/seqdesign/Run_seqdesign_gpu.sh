export 'PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:1000'
conda activate /data/users/athchu/.conda/envs/Seqdesign
echo "set enable-bracketed-paste off" >> ~/.inputrc
program_dir=data/zeroshot_prediction/Seqdesign/src/seqdesign_pt/scripts
export PYTHONPATH="${PYTHONPATH}:/path_to/seqdesign-pytorch/src"
#LD_LIBRARY_PATH=/data/users/athchu/.conda/envs/Seqdesign/lib/python3.7/site-packages/nvidia/cublas/lib/:$LD_LIBRARY_PATH
#CUDNN_H_PATH=/usr/include/cudnn.h
cd data/zeroshot_prediction/Seqdesign/SpCas9
$program_dir/run_autoregressive_fr.py --gpu 1 \
--alphabet-type protein \
--batch-size 10 \
--num-iterations 50000 \
--dataset SPCAS9_e40_w1.fasta &
cd data/zeroshot_prediction/Seqdesign/SaCas9 
$program_dir/run_autoregressive_fr.py --gpu 2 \
--alphabet-type protein \
--batch-size 10 \
--num-iterations 50000 \
--dataset SaCas9_J7RUa5_w1.fasta &
cd data/zeroshot_prediction/Seqdesign/PABP
$program_dir/run_autoregressive_fr.py --gpu 3 \
--alphabet-type protein \
--batch-size 10 \
--num-iterations 50000 \
--dataset PABP_homologues_w1.fasta &

cd data/zeroshot_prediction/Seqdesign/PhoQ && 
$program_dir/run_autoregressive_fr.py --gpu 1 \
--alphabet-type protein \
--batch-size 10 \
--num-iterations 50000 \
--dataset PhoQ_homologues_w1.fasta &


### infer mutations effect
cd data/zeroshot_prediction/Seqdesign/GFP 
$program_dir/run_autoregressive_fr.py --gpu 2 \
--alphabet-type protein \
--batch-size 10 \
--num-iterations 50000 \
--dataset GFP_homologues_w1.fasta &


cd data/zeroshot_prediction/Seqdesign/UBE 
$program_dir/run_autoregressive_fr.py --gpu 3 \
--alphabet-type protein \
--batch-size 10 \
--num-iterations 50000 \
--dataset UBE_homologues_w1.fasta &


cd data/zeroshot_prediction/Seqdesign/SpCas9
$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/SpCas9 \
 --sess SPCAS9_e40_w1_v3-pt_channels-48_rseed-42_23Jun07_1216PM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 100 \
 --input data/zeroshot_prediction/Seqdesign/SpCas9/SpCas9_AACombo.fasta \
 --output SpCas9_mutant_100 &

cd data/zeroshot_prediction/Seqdesign/SaCas9
$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/SaCas9 \
 --sess SaCas9_J7RUa5_w1_v3-pt_channels-48_rseed-42_23Jun07_1217PM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 100 \
 --input SaCas9_887888887985986988989991.fasta \
 --output SaCas9_887888887985986988989991_mutant_100 &

 cd data/zeroshot_prediction/Seqdesign/SaCas9
$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/SaCas9 \
 --sess SaCas9_J7RUa5_w1_v3-pt_channels-48_rseed-42_23Jun07_1217PM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 100 \
 --input SaCAs9_888889909_mutants.fasta \
 --output SaCas9_888889909_mutant_100 &


$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/PhoQ \
 --sess PhoQ_homologues_w1_v3-pt_channels-48_rseed-42_23Jun06_0855PM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 10 \
 --input data/zeroshot_prediction/Seqdesign/PhoQ/PhoQ_mutants.fasta \
 --output PhoQ_mutant_100_1 &

cd data/zeroshot_prediction/Seqdesign/PABP
$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/PABP \
 --sess PABP_homologues_w1_v3-pt_channels-48_rseed-42_23Jun07_1236PM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 10 \
 --input PABP_mutants.fasta \
 --output PABP_100_mutants &

cd data/zeroshot_prediction/Seqdesign/UBE
$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/UBE \
 --sess UBE_homologues_w1_v3-pt_channels-48_rseed-42_23Jun08_0310AM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 10 \
 --input UBE_mutants.fasta \
 --output UBE_100_mutants &

cd data/zeroshot_prediction/Seqdesign/GFP
$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/GFP \
 --sess GFP_homologues_w1_v3-pt_channels-48_rseed-42_23Jun08_0309AM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 10 \
 --input GFP_mutants.fasta \
 --output GFP_100_mutants &

cd data/zeroshot_prediction/Seqdesign/evoAPOBEC
$program_dir/20230602_calc_logprobs_seqs_fr.py \
 --wkdir data/zeroshot_prediction/Seqdesign/evoAPOBEC \
 --sess EvoAPOBEC_homologues_w1_v3-pt_channels-48_rseed-42_23Jun08_1224PM \
 --checkpoint 50000 \
 --dropout-p 0.5 \
 --num-samples 100 \
 --minibatch-size 10 \
 --input evoAPROBEC_AACombo.fasta \
 --output evoAPOBEC_100_mutants1 &
