# 20190412 Updated to include new locations for executables
# Apr 12, 2018 Kai Xun Joshua Tay
# (1) processes smart3seq fastq with umipolymer.py
# (2) maps to a combined hg19 and ebv genome
# runs by:
# python map_hg19_and_ebv.py filelist.txt
# where filelist.txt is a tab-delimited file with <file.fastq.gz> in the first column and <sample_ID> at the 2nd column of each line


import sys  # import the necessary modules
import os.path
import time

sample_fastq = sys.argv[1]  # get file with sample names and fastq

datetime = time.strftime("%Y%m%d-%H%M")

STAR_executable = "/media/data/Josh/packages/STAR-2.7.0f/bin/Linux_x86_64/STAR"
STAR_hg19_ebv_genomeDir = "/media/data/Josh/NPC_star_indices/hg19ebv_starindex68"  # NB: these must be prepared for readlengths 68
UMI_HOMOPOLOYMER_script = "/media/data/Josh/packages/Smart-3SEQ/umi_homopolymer.py"
DEDUP_script = "/media/data/Josh/packages/Smart-3SEQ/umi-dedup-master/dedup.py"
Root = "/media/data/Josh/NPC_s3s_all/"
TEMP_Dir = Root + datetime + "_temp"
OUTPUT_Dir = Root + datetime + "_hg19andEBV"

print("Making output directories")
make_outputdir_command = "mkdir " + OUTPUT_Dir
os.system(make_outputdir_command)
make_tempdir_command = "mkdir " + TEMP_Dir
os.system(make_tempdir_command)

# load genome references for STAR

print("Loading STAR genome references")
load_genome_command = STAR_executable + " --genomeLoad LoadAndExit --genomeDir " + STAR_hg19_ebv_genomeDir
print(load_genome_command)
os.system(load_genome_command)


with open(sample_fastq, "r") as sample_fastq_file:

    for line in sample_fastq_file:
        line = line.strip()
        line_list = line.split('\t')
        sample_no = line_list[1]
        fastq_path = line_list[0]

        print(sample_no, fastq_path)

        print("Copying fastq")
        copy_fastq_command = "cp " + fastq_path + " " + TEMP_Dir
        print(copy_fastq_command)
        os.system(copy_fastq_command)


        print("Unzipping fastq")
        fastq_path_list = fastq_path.split('/')
        fastq_file_name_only = fastq_path_list[-1]
        fastq_location = TEMP_Dir + "/" + fastq_file_name_only
        unzip_fastq_command = "pigz -d " + fastq_location
        print(unzip_fastq_command)
        os.system(unzip_fastq_command)


        fastq_input_file = fastq_location[:-3]   # remove .gz from the location
        fastq_uh_file = fastq_input_file + ".uh" # location of umi_homopolymer processed fastq file

        # run UMI homopolymer

        print("Running UMI homopolymer script")
        umi_homopolymer_comamnd = "pypy " + UMI_HOMOPOLOYMER_script + " -u 5 -g 3 -p 8 -m 1 " + fastq_input_file + " " + fastq_uh_file + " 2> " + OUTPUT_Dir + "/" + sample_no + "_umihomopolymer.log"
        print(umi_homopolymer_comamnd)
        os.system(umi_homopolymer_comamnd)


        # STAR: map to human and EBV reference genome together
        # --genomeLoad LoadAndKeep --genomeDir /media/joannaprzybyl/SSD-Data/rna_npc_smart3SEQ/dbsnp147_gencode25-68 --readFilesIn testrun.fastq --runThreadN 38 --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outBAMcompression 10 --outReadsUnmapped Fastx --limitBAMsortRAM 2147483648 --outFilterMultimapNmax 1 --outFilterMismatchNmax 999 --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA --clip3pAdapterMMp 0.2


        print("Mapping to Combined genome")

        STAR_cmd_human = STAR_executable + \
            " --genomeLoad LoadAndKeep" + \
            " --genomeDir " + STAR_hg19_ebv_genomeDir + \
            " --runThreadN 38" + \
            " --readFilesIn " + fastq_uh_file + \
            " --outFileNamePrefix " + OUTPUT_Dir + "/" + sample_no + "_hg19ebv_" + \
            " --outSAMtype BAM SortedByCoordinate" + \
            " --quantMode GeneCounts" + \
            " --outBAMcompression 0" + \
            " --limitBAMsortRAM 2147483648" + \
            " --outFilterMultimapNmax 1" + \
            " --outFilterMismatchNmax 999" + \
            " --clip3pAdapterSeq AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + \
            " --clip3pAdapterMMp 0.2"

        print(STAR_cmd_human)
        os.system(STAR_cmd_human)

        # Remove temp files

        print("Removing temp files")
        remove_files_temp = "rm " + TEMP_Dir + "/*.*"
        print(remove_files_temp)
        os.system(remove_files_temp)


print("Unloading STAR genome references")
unload_genome_command = STAR_executable + " --genomeLoad Remove --genomeDir " + STAR_hg19_ebv_genomeDir
os.system(unload_genome_command)

print("Deduplicating bam files in parallel")
dedup_cmd = "ls " + OUTPUT_Dir + "/*.bam | parallel python " + DEDUP_script + """ -s '{}' '{.}'.dedup.bam ">" '{.}'.dedup.log"""
print(dedup_cmd)
os.system(dedup_cmd)


print("Creating BAM indexes in parallel")
bam_index_command = "ls " + OUTPUT_Dir + "/*dedup.bam | parallel samtools index '{}'"
print(bam_index_command)
os.system(bam_index_command)


