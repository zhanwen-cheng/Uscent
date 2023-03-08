#!/bin/bash

#brdu process
#basic used software basecalling---->cat fastq---->porechop---->minimap2---->DNAscent index---->DNAscent detect---->calculate depth file---->plot
show_help() {
cat << EOF
Usage: ${0##*/}

Process:
        basecalling using guppy
                trime the fastq with porechop
                align the read to the genome with minimap2
                generate depth file
                DNAscent to detect all U position
                calculte accumulated U/bases every 1 kb
                plot in python inside PDF

Example usage:
        bash brdu.bash -f fast5_file_path -g genome.fa -o output_dir

NOTICE:
        Make sure only fasta format are supported here
        !!!!YOU ARE ENCOVERAGED TO USE FULL PATH !!!!

EOF
}
while getopts "f:g:o:h" opt; do
        case "$opt" in
                h | --help)
                        show_help
                        exit 0
                        ;;
                f | --fast5)
                        Fast5=$OPTARG
                        ;;
                g | --genome)
                        Input_genome=$OPTARG
                        ;;
                o | --outdir)
                        Out_dir=$OPTARG
                        ;;
                \?)
                        echo "Invalid option: -$OPTARG" >&2
                        ;;
                '?')
                        show_help >&2
                        exit 1
                        ;;
                -?*)
                        print 'Warning: Unknown option (ignored) : %s\n' "$1" >&2
                        ;;
                :)
                        echo "Option -$OPTARG requires an argument." >&2
                        exit 1
                        ;;
                *) # default case: if no more options then break out of the loop
                        break
        esac

done
if [ -z "$Input_genome" ]
then
        echo "No input fasta, -f must be specified"
        exit
fi

if [ ! -d "$Out_dir" ]
then
        mkdir $Out_dir
        else
        echo "Output file already exist"
fi

working_path=$(pwd)
echo "Fast5_file = ${working_path}/${Fast5}"
echo "Genome     = ${working_path}/${Input_genome}"
echo "Output_dir = ${working_path}/${Out_dir}"

####guppy basecalling ####
echo '#############################################################################'
echo '##########################basecalling with guppy#############################'
echo '#############################################################################'
mkdir $Out_dir/guppy
~/software/ont-guppy/bin/guppy_basecaller -i ${working_path}/${Fast5} -c ~/software/ont-guppy/data/dna_r9.4.1_e8.1_fast.cfg -s ${working_path}/${Out_dir}/guppy -x "cuda:0,1"
cat ${working_path}/${Out_dir}/guppy/pass/*.fastq ${working_path}/${Out_dir}/guppy/fail/*.fastq > ${working_path}/${Out_dir}/guppy/merge.fastq

# ####porechop####
# echo '#############################################################################'
# echo '###############################porechop step#################################'
# echo '#############################################################################'
# porechop -i ${working_path}/${Out_dir}/guppy/merge.fastq -o ${working_path}/${Out_dir}/guppy/merge_porechop.fastq

####minimap2####
echo '#############################################################################'
echo '##########################align with minimap2################################'
echo '#############################################################################'
mkdir $Out_dir/minimap2
cd $Out_dir/minimap2
#minimap2 -ax map-ont -t 10 -o alignment.sam ${working_path}/${Input_genome} ${working_path}/${Out_dir}/guppy/merge_porechop.fastq
minimap2 -ax map-ont -t 10 -o alignment.sam ${working_path}/${Input_genome} ${working_path}/${Out_dir}/guppy/merge.fastq
samtools view -Sb -@ 10 -o alignment.bam alignment.sam
samtools sort -@ 10 alignment.bam -o alignment.sorted.bam
samtools index -@ 10 alignment.sorted.bam
samtools depth alignment.sorted.bam > depth.txt
cd $working_path

####DNAscent####
echo '#############################################################################'
echo '##########################DNAscent detetct U#################################'
echo '#############################################################################'
mkdir $Out_dir/DNAscent
cd ./$Out_dir/DNAscent
DNAscent index -f ${working_path}/${Fast5} -s ${working_path}/${Out_dir}/guppy/sequencing_summary.txt
DNAscent detect -b ${working_path}/$Out_dir/minimap2/alignment.sorted.bam -r ${working_path}/${Input_genome} -i ./index.dnascent -o ${Out_dir}.detect -t 10
DNAscent psl -d ${Out_dir}.detect -r ${working_path}/${Input_genome} -o ${working_path}/${Out_dir}/DNAscent/${Out_dir} --threshold 0.7
python ~/software/script/calculate_1kb_gc.py -r ${working_path}/${Input_genome} -o gc.csv
cd $working_path

####iRep####
echo '#############################################################################'
echo '############################iRep analysis####################################'
echo '#############################################################################'
mkdir $Out_dir/iRep
cd $Out_dir/iRep
iRep -f ${working_path}/${Input_genome} -t 10 -s ${working_path}/${Out_dir}/minimap2/alignment.sam -o ${working_path}/${Out_dir}/iRep/$Out_dir
cd $working_path

####calculate U depth and depth file####
echo '#############################################################################'
echo '########################calculate depth file#################################'
echo '#############################################################################'
mkdir $Out_dir/depth
python ~/software/script/calculate_depth.py -p ${working_path}/${Out_dir}/DNAscent/${Out_dir}.psl -d ${working_path}/${Out_dir}/minimap2/depth.txt -o ${working_path}/${Out_dir}/depth/${Out_dir}.csv
python ~/software/script/calculate_1kb_depth.py -d ${working_path}/${Out_dir}/depth/${Out_dir}.csv -o ${working_path}/${Out_dir}/depth/${Out_dir}_1kb.csv
cd $working_path

####normalization fit####
echo '#############################################################################'
echo '####################fit with normal distribution#############################'
echo '#############################################################################'
mkdir $Out_dir/nordis_fit
cd $Out_dir/nordis_fit
python ~/software/script/Bcounts_fit_normal+linear_dis2_trimmed_mean.py -d ${working_path}/${Out_dir}/depth/${Out_dir}_1kb.csv -g ${working_path}/${Out_dir}/DNAscent/gc.csv -o B_${Out_dir}.pdf
python ~/software/script/Ucounts_fit_normal+linear_dis2_trimmed_mean.py -d ${working_path}/${Out_dir}/depth/${Out_dir}_1kb.csv -g ${working_path}/${Out_dir}/DNAscent/gc.csv -o U_${Out_dir}.pdf
cd $working_path
