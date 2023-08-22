import argparse
from collections import defaultdict
import subprocess
import re

def cmd_out(cmd):
    return subprocess.check_output(cmd, shell=True, text=True).splitlines()

def write_coverage_to_file(bam_file, bed_file, output_file):
    genome_coverage = []

    for l in cmd_out(f"samtools view -b -L {bed_file} {bam_file} | bedtools coverage -a {bed_file} -b - -sorted "):
        row = l.split()
        if "Rv" in row[3]:
            gene = row[3]
            
            if "PE" in row[9]:
                matches = re.findall("PE\d+|PE_|pe_", row[9])
                if matches:
                    ppeBefore = row[9].split(";")
                    ppeBeforeAgain = ppeBefore[1].split("=")
                    ppe = ppeBeforeAgain[1]
            else :
                ppe = "."
            genomic_position = int(row[1]) + int(row[-2]) -1
            depth = float(row[-1])
            genome_coverage.append((genomic_position,ppe,gene,depth))

    with open(output_file, 'w') as f:
        for position, gene, ppe,depth in genome_coverage:
            f.write(f"{position}\t{gene}\t{ppe}\t{depth}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compute coverage for a given BAM file.')
    parser.add_argument('bam_file', help='The input BAM file')
    parser.add_argument('--bed_file', default='/mnt/storage9/joseph/miniconda3/envs/medaka/share/tbprofiler/tbdb.bed', help='The BED file to use (default: %(default)s)')
    parser.add_argument('--output_file', default='LETSLOOK.txt', help='The output file to write to (default: %(default)s)')
    
    args = parser.parse_args()
    
    write_coverage_to_file(bam_file=args.bam_file, bed_file=args.bed_file, output_file=args.output_file)
