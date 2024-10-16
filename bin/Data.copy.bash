
DIR=/home/chenzh/My_project/SHI_Glu
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin

BDOC=$DIR/big_doc

cd $DIR

mkdir -p $TMPD/raw_data

paste <(\ls $DATA/Project_sunlu_20210813NA/C*_R1_*fastq.gz) <(\ls $DATA/Project_sunlu_20210813NA/C*_R1_*fastq.gz |cut -f 8 -d "/"|cut -f1-2 -d "_"|awk '{print "/home/chenzh/My_project/SHI_Glu/tmp_data/raw_data/"$0"_R1.fastq.gz"}')|awk '{print "ln -s "$1,$2}' > $BIN/temp.run.bash

paste <(\ls $DATA/Project_sunlu_20210813NA/C*_R2_*fastq.gz) <(\ls $DATA/Project_sunlu_20210813NA/C*_R2_*fastq.gz |cut -f 8 -d "/"|cut -f1-2 -d "_"|awk '{print "/home/chenzh/My_project/SHI_Glu/tmp_data/raw_data/"$0"_R2.fastq.gz"}')|awk '{print "ln -s "$1,$2}' >> $BIN/temp.run.bash


paste <(\ls $DATA/Project_sunlu_20210813NA/[!C]*_R1_*fastq.gz) <(\ls $DATA/Project_sunlu_20210813NA/[!C]*_R1_*fastq.gz |cut -f 8 -d "/"|cut -f1-3 -d "_"|awk '{print "/home/chenzh/My_project/SHI_Glu/tmp_data/raw_data/"$0"_R1.fastq.gz"}')|awk '{print "ln -s "$1,$2}' >> $BIN/temp.run.bash
paste <(\ls $DATA/Project_sunlu_20210813NA/[!C]*_R2_*fastq.gz) <(\ls $DATA/Project_sunlu_20210813NA/[!C]*_R2_*fastq.gz |cut -f 8 -d "/"|cut -f1-3 -d "_"|awk '{print "/home/chenzh/My_project/SHI_Glu/tmp_data/raw_data/"$0"_R2.fastq.gz"}')|awk '{print "ln -s "$1,$2}' >> $BIN/temp.run.bash


paste <(\ls $DATA/Project_sunlu_20211228NA/[!C]*_R1_*fastq.gz) <(\ls $DATA/Project_sunlu_20211228NA/[!C]*_R1_*fastq.gz |cut -f 8 -d "/"|cut -f1-2 -d "_"|awk '{print "/home/chenzh/My_project/SHI_Glu/tmp_data/raw_data/"$0"_R1.fastq.gz"}')|awk '{print "ln -s "$1,$2}' >> $BIN/temp.run.bash
paste <(\ls $DATA/Project_sunlu_20211228NA/[!C]*_R2_*fastq.gz) <(\ls $DATA/Project_sunlu_20211228NA/[!C]*_R2_*fastq.gz |cut -f 8 -d "/"|cut -f1-2 -d "_"|awk '{print "/home/chenzh/My_project/SHI_Glu/tmp_data/raw_data/"$0"_R2.fastq.gz"}')|awk '{print "ln -s "$1,$2}' >> $BIN/temp.run.bash
