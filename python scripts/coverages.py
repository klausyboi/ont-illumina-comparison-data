import argparse
from collections import defaultdict
import subprocess
import re

def cmd_out(cmd):
    """
    Executes a shell command and returns the output as a list of lines.

    :param cmd: A string representing the shell command to execute.
    :return: A list of strings, each representing a line of output from the command.
    """
    return subprocess.check_output(cmd, shell=True, text=True).splitlines()

def write_coverage_to_file(bam_file, bed_file, output_file):
    """
    Computes coverage for regions specified in a BED file for a given BAM file 
    and writes the coverage data to an output file.

    :param bam_file: Path to the input BAM file.
    :param bed_file: Path to the BED file specifying regions of interest.
    :param output_file: Path to the output file where the coverage data will be written.
    """
    genome_coverage = []

    ##iterate over the output lines of the coverage computation command
    for l in cmd_out(f"samtools view -b -L {bed_file} {bam_file} | bedtools coverage -a {bed_file} -b - -sorted"):
        row = l.split()
        ##process only lines containing the "Rv" marker
        if "Rv" in row[3]:
            gene = row[3]
            ##extract PPE information, if available
            if "PE" in row[9]:
                matches = re.findall("PE\d+|PE_|pe_", row[9])
                if matches:
                    ppeBefore = row[9].split(";")
                    ppeBeforeAgain = ppeBefore[1].split("=")
                    ppe = ppeBeforeAgain[1]
            else:
                ppe = "."
            genomic_position = int(row[1]) + int(row[-2]) - 1
            depth = float(row[-1])
            genome_coverage.append((genomic_position, ppe, gene, depth))

    ##write the coverage data to the output file
    with open(output_file, 'w') as f:
        for position, gene, ppe, depth in genome_coverage:
            f.write(f"{position}\t{gene}\t{ppe}\t{depth}\n")

if __name__ == "__main__":
    ##set up argument parsing for command-line execution
    parser = argparse.ArgumentParser(description='Compute coverage for a given BAM file.')
    parser.add_argument('bam_file', help='The input BAM file')
    parser.add_argument('--bed_file', default='/mnt/storage9/joseph/miniconda3/envs/medaka/share/tbprofiler/tbdb.bed', help='The BED file to use (default: %(default)s)')
    parser.add_argument('--output_file', default='LETSLOOK.txt', help='The output file to write to (default: %(default)s)')
    
    args = parser.parse_args()
    
    ##execute the coverage writing function with provided arguments
    write_coverage_to_file(bam_file=args.bam_file, bed_file=args.bed_file, output_file=args.output_file)
