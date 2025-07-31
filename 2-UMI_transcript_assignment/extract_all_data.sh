data_dir=/export/home1/rstudio-homes/hamraoui/data/simulated_R10/
#Sicelore 2.1
bash bin/extract_data.sh -b /04a.matrices/PCR_1_isobam.bam -t Sicelore -p PCR_1 -o $data_dir
bash bin/extract_data.sh -b /04a.matrices/PCR_5_isobam.bam -t Sicelore -p PCR_5 -o $data_dir

#Sockeye
bash bin/extract_data.sh -b ~/simulated/simulated.tagged.bam -t Sockeye -p PCR_1 -o $data_dir
bash bin/extract_data.sh -b ~/simulated/simulated.tagged.bam -t Sockeye -p PCR_5 -o $data_dir

#FLAMES
bash bin/extract_data.sh -b ~/scLRPipe/realign2transcript.bam -t FLAMES -p PCR_1 -o $data_dir
bash bin/extract_data.sh -b ~/scLRPipe/realign2transcript.bam -t FLAMES -p PCR_5 -o $data_dir

#Bambu
bash bin/extract_data.sh -b ~/Bambu.demultiplexed.bam -t bambu -p PCR_1 -o $data_dir
bash bin/extract_data.sh -b ~/Bambu.demultiplexed.bam -t bambu -p PCR_5 -o $data_dir

