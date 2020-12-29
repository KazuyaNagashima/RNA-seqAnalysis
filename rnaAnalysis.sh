#!/bin/env zsh -e
cd ..
echo 'SRAの名前を入力してください(.sraの前)'
read accession

touch ${accession}_analysis.detail

echo 'リファレンスゲノムを入力して下さい[mm9,mm10]'
read ref

data=(`grep ${accession} data.csv`)

arr=( `echo ${data} | tr -s ',' ' '`)
sra=(${arr[1]})
sample=(${arr[2]})
seq=(${arr[3]})
ends=(${arr[4]})
reserch=(${arr[5]})


pfastq-dump -v
fastp -v
star --version
rsem-calculate-expression -version

echo '\n <<sample detail>>' | tee -a ${accession}_analysis.detail
echo 'SRA:'${sra} | tee -a ${accession}_analysis.detail
echo 'SAMPLE:'${sample} | tee -a ${accession}_analysis.detail
echo 'SEQ:'${seq} | tee -a ${accession}_analysis.detail
echo 'READ:'${ends} | tee -a ${accession}_analysis.detail
echo 'RESERCH:'${reserch} | tee -a ${accession}_analysis.detail
echo 'REFERENCE:'${ref} | tee -a ${accession}_analysis.detail

echo '以上の設定でRNA-seq解析を行いますか？[y/n]'
read qestion
if [ $qestion = 'y' ]; then
  echo 'プログラムを実行します'
else
  echo 'プログラムを終了します'
  rm ${accession}_analysis.detail
  exit
fi

echo '@edu.k.u-tokyo.ac.jpのアドレスを設定すると、プログラム終了時にメールが届きます。メールを設定しますか？[y/n]'
read mail
if [ $mail = 'y' ]; then
  echo '@edu.k.u-tokyo.ac.jpの前を入力してください'
  read address
else
  echo 'メールを設定しません'
fi

dir=(`pwd`)


mkdir ${dir}/${reserch}
mkdir ${dir}/${reserch}/${sample}_${sra}
mkdir ${dir}/${reserch}/${sample}_${sra}/${ref}
cd ${dir}/${reserch}/${sample}_${sra}/${ref}


#pfastq-dumpでSRAをfastqに変換
if [ $ends = 'single' ]; then
   pfastq-dump --threads 8  --gzip -s ${dir}/sra/${sra}.sra
elif [ $ends = 'pair' ]; then
   #Pair Endsの時は--split-filesをつける
   pfastq-dump --threads 8 --gzip --split-files -s ${dir}/sra/${sra}.sra
else
   echo 'readをsingleかpairに設定してください'
   echo 'プログラムを終了します'
   exit
fi

mkdir fastq
if [ $ends = 'single' ]; then
   mv ${sra}_.fastq.gz fastq
else
   mv ${sra}_1.fastq.gz ${sra}_2.fastq.gz fastq
fi

#fastpによるQuality Check(20bp>,q>20)
mkdir fastp_reports
if [ $ends = 'single' ]; then
   fastp -i ./fastq/${sra}_.fastq.gz -o ./fastp_reports/fastp_${sra}.fastq -h ./fastp_reports/report_fastp.html -j ./fastp_reports/report_fastp.json -q 20 --length_required 20
else
   fastp -i ./fastq/${sra}_1.fastq.gz -I ./fastq/${sra}_2.fastq.gz -o ./fastp_reports/fastp_${sra}R1.fastq -O ./fastp_reports/fastp_${sra}R2.fastq -h ./fastp_reports/report_fastp.html -j ./fastp_reports/report_fastp.json -q 20 --length_required 20
fi

#STARによるマッピング
mkdir star
if [ $ends = 'single' ]; then
    STAR --runThreadN 8 \
     --runMode alignReads \
     --genomeDir /Users/shigenseigyo/Desktop/reference/Mus_musculus/UCSC/${ref}/Sequence/StarIndex \
     --quantMode TranscriptomeSAM GeneCounts \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn ./fastp_reports/fastp_${accession}.fastq \
     --outFileNamePrefix ./star/star_${accession}

    rm ./fastp_reports/fastp_${sra}.fastq
else
    STAR --runThreadN 8 \
     --runMode alignReads \
     --genomeDir /Users/shigenseigyo/Desktop/reference/Mus_musculus/UCSC/${ref}/Sequence/StarIndex \
     --quantMode TranscriptomeSAM GeneCounts \
     --outSAMtype BAM SortedByCoordinate \
     --readFilesIn ./fastp_reports/fastp_${accession}R1.fastq ./fastp_reports/fastp_${accession}R2.fastq \
     --outFileNamePrefix ./star/star_${accession}

    rm ./fastp_reports/fastp_${sra}R1.fastq ./fastp_reports/fastp_${sra}R2.fastq
fi

#RSEMによる遺伝子発現定量
mkdir rsem
cd rsem
if [ $ends = 'single' ]; then
    rsem-calculate-expression \
    --alignments \
    -p 8 \
    ../star/star_${accession}Aligned.toTranscriptome.out.bam \
    /Users/shigenseigyo/Desktop/reference/Mus_musculus/UCSC/${ref}/Sequence/RsemReference/RsemReference \
    ${accession}
else
    rsem-calculate-expression \
    --alignments --paired-end \
    -p 8 \
    ../star/star_${accession}Aligned.toTranscriptome.out.bam \
    /Users/shigenseigyo/Desktop/reference/Mus_musculus/UCSC/${ref}/Sequence/RsemReference/RsemReference \
    ${accession}
fi





#終了時メールで送信
if [ $mail = 'y' ]; then
  echo "Process done　メールを送信しました。" | mail -s "Process done" ${address}@edu.k.u-tokyo.ac.jp
else
  echo "Process done"
fi


mv ${dir}/${accession}_analysis.detail ..
