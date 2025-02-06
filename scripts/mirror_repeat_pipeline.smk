import subprocess
import pandas as pd
import re
import csv
from utils import parse_fasta
import os

configfile: 'mr_extractions.yaml'

outdir = Path(config['outdir']).resolve()
outdir.mkdir(exist_ok=True)
indir = Path(config['indir']).resolve()

min_arm_length = int(config['min_arm_length'])
max_spacer_length = int(config['max_spacer_length'])
gfa='/storage/group/izg5139/default/external/quadrupia_database/g4/non-B_gfa/gfa'

PRIMATES = glob_wildcards('%s/{accession}.fna' % indir).accession
print(PRIMATES)

MINDI_FIELDS = [
                                         "seqID",
                                         "start",
                                         "end",
                                         "sequenceOfArm",
                                         "sequenceOfSpacer",
                                         "sequence",
                                         "armLength",
                                         "spacerLength",
                                         "length",
                                         "arm_a",
                                         "arm_c",
                                         "arm_g",
                                         "arm_t"
                                         ]

rule all:
    input:
        expand('%s/{accession}_MR.tsv' % outdir, accession=PRIMATES)


rule extract_mirrors:
    input:
        '%s/{accession}.fna' % indir
    output:
        '%s/{accession}_MR.tsv' % outdir,
        '%s/{accession}_MR.processed.tsv' % outdir
    run:
        # infile = f"{outdir}/{wildcards.accession}.fna"
        # if Path(infile).is_file():
        # raise ValueError("GUNZIP FILE TARGET EXISTS!")
        
        # command = f"gunzip -k -f -c {input} > {infile}"
        # subprocess.run(command, shell=True, check=True)
        print(outdir)
        prev_dir = os.getcwd()
        os.chdir(outdir)
        print(input[0])
        command = f"{gfa} -seq {input[0]} -out {wildcards.accession} -maxMRspacer 100 -minMRrep 8 -skipAPR -skipSTR -skipIR -skipDR -skipGQ -skipZ -skipSlipped -skipCruciform -skipTriplex -skipWGET"
        # command = f"{gfa} -seq {input[0]} -out {wildcards.accession} -skipAPR -skipMR -skipIR -minSTR 1 -maxSTR 9 -minSTRbp 8 -skipDR -skipGQ -skipZ -skipSlipped -skipCruciform -skipTriplex -skipWGET"
        print(command)
        subprocess.run(command, shell=True, check=True)
        os.remove(wildcards.accession + "_MR.gff")
        os.chdir(prev_dir)
        # os.remove(infile)

        to_drop = ["Source",
                   "Type",
                   "Score",
                   "Strand",
                   "Subset",
                   "Permutations",
                   "Sequence",
                   "Start",
                   "Stop"]
        nucleotides = {'a', 'g', 'c', 't'}

        FRAME_FIELDS = [
                             "seqID",
                             "start",
                             "end",
                             "sequenceOfArm",
                             "sequenceOfSpacer",
                             "sequence",
                             "armLength",
                             "spacerLength",
                             "sequenceLength",
                             "arm_a",
                             "arm_g",
                             "arm_c",
                             "arm_t",
                             "composition"
                          ]

        with open(output[1], mode="w", encoding="UTF-8") as fout:
            fout_writer = csv.DictWriter(
                                         fout,
                                         delimiter="\t",
                                         fieldnames=FRAME_FIELDS
                                         )

            fout_writer.writeheader()

            with open(output[0], mode="r", encoding="UTF-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")

                for row in reader:
                    arm_length = int(row['Repeat'])
                    spacer_length = int(row['Spacer'])
                    sequence = row['Sequence']

                    if (isinstance(min_arm_length, int) and arm_length < min_arm_length) or (isinstance(max_spacer_length, int) and spacer_length > max_spacer_length):
                        continue

                    start = int(row['Start']) - 1
                    end = int(row['Stop'])
                    sequence_length = int(row['Length'])
                    total_coordinate_length = end - start

                    if sequence_length < total_coordinate_length:
                        sequence = sequence[:sequence_length]
                        end = end - (total_coordinate_length - sequence_length)
                        print(colored(f'Invalid sequence length detected for {fn} on (start,end)=({start},{end}) with sequence {sequence}.', 'red'))
                        raise ValueError(f'Invalid sequence length detected for {fn}.')


                    # find sequence of arm
                    repeat = int(row['Repeat'])
                    sequence_of_arm = sequence[:repeat]

                    # find spacer
                    del row['Spacer']
                    true_spacer_length = sequence_length - 2 * repeat
                    right_arm = sequence[repeat+true_spacer_length:]


                    if any(n not in nucleotides for n in right_arm) or any(n not in nucleotides for n in sequence_of_arm):
                        continue

                    # spacer = sequence[repeat:repeat+spacer_length]
                    true_spacer = sequence[repeat:repeat+true_spacer_length]
                    if len(true_spacer) == 0:
                        true_spacer = "."

                    # process composition
                    composition = re.search(r"(\d+)A/(\d+)C/(\d+)G/(\d+)T", row['Composition'])

                    a_content = composition.group(1)
                    c_content = composition.group(2)
                    g_content = composition.group(3)
                    t_content = composition.group(4)

                    # skip maximum spacer length
                    if (isinstance(max_spacer_length, int) and true_spacer_length > max_spacer_length):

                        # invalid record?
                        continue

                    row.update({
                            "start": start,
                            "end": end,
                            "sequenceOfArm": sequence_of_arm,
                            "sequenceOfSpacer": true_spacer,
                            "spacer": true_spacer_length,
                            "sequence": sequence,
                            "arm_a": a_content,
                            "arm_c": c_content,
                            "arm_g": g_content,
                            "arm_t": t_content,
                        })

                    for col in to_drop:
                        del row[col]

                    row["sequenceLength"] = row.pop("Length")
                    row["seqID"] = row.pop("Sequence_name")
                    row["armLength"] = row.pop("Repeat")
                    row["spacerLength"] = row.pop("spacer")
                    row["composition"] = row.pop("Composition")
                    fout_writer.writerow(row)

        df = pd.read_table(output[1])
        for seqID, seq in parse_fasta(input[0]):
            temp = df[df['seqID'] == seqID]

            total_counts = 0
            total = temp.shape[0]
            if total == 0:
                print(f"SeqID {seqID} for input file {input[0]} is empty.")
            for _, row in temp.iterrows():
                  start = int(row['start'])
                  end = int(row['end'])
                  sequence = row['sequence']
                  sequenceOfArm = row['sequenceOfArm']
                  sequenceOfSpacer = row['sequenceOfSpacer']
                  arm_length = int(row['armLength'])
                  spacer_length = int(row['spacerLength'])
                  mirror_length = int(row['sequenceLength'])

                  original_sequence = seq[start: end]
                  assert arm_length * 2 + spacer_length == len(sequence) == mirror_length

                  if sequenceOfSpacer == '.':
                    sequenceOfSpacer = ''
                  assert original_sequence == sequence == sequenceOfArm + sequenceOfSpacer + sequenceOfArm[::-1]
                  total_counts += 1
            assert total_counts == total
            print(f"{seqID} is OK!")
