#! /usr/bin/env python

import sys
import numpy as np

helpMsg = '''

    usage: $./fasta_MSA2SITES.py <msa.fasta> <out_pref> <`all`/`incl`/`excl`> <acc_ls.txt>

        Outputs ARGweaver SITES format

'''

bases = ["A", "G", "C", "T", "N"]

def main(args):
    if len(args) != 5:    #4 argument(s)
        return helpMsg

    fasta_path = args[1]
    outPref = args[2]
    mode = args[3]

    with open(args[4], "r") as sampF:
        samp_acc = sampF.read().strip().split("\n")

    ID_ls = [] # list of identifier with corresponding order to MSA_mtx
    MSA_mtx = [] # sample by position

    with open(fasta_path, "r") as fastaF:
        cntr = 0
        buf_str = ''
        for line in fastaF:
            if line[0] == ">":
                if len(buf_str) != 0:
                    print(f"Read seq no. {cntr}: {name}, length {len(buf_str)}", flush=True)
                    ID_ls.append(name)
                    MSA_mtx.append(list(buf_str))
                    cntr += 1
                    buf_str = ''

                name, acc_no, _ = line[1:].split("|")
                in_list = (acc_no in samp_acc)

                if mode == 'all':
                    read_seq = True
                elif mode == 'incl' and in_list:
                    read_seq = True
                elif mode == 'excl' and not in_list:
                    read_seq = True
                else:
                    read_seq = False
                    print(f"Skip seq no. {cntr}: {name}", flush=True)

            elif read_seq:
                buf_str += line.strip().upper()

        if len(buf_str) != 0:
            print(f"Read seq no. {cntr}: {name}, length {len(buf_str)}", flush=True)
            ID_ls.append(name)
            MSA_mtx.append(list(buf_str))

    MSA_mtx = np.transpose(np.array(MSA_mtx)) # np matrix, position by sample

    print(f"Finished reading MSA: {MSA_mtx.shape}")

    seq_len = MSA_mtx.shape[0]

    # Cases: 0-Invariant; 1-Unknown; 2-Ambiguous/Gap; 3-Variant
    tally = [0, 0, 0, 0, 0]

    with open(outPref+".sites", "w") as sitesF:
        sitesF.write("NAMES\t")
        sitesF.write("\t".join(ID_ls)+"\n")
        sitesF.write(f"REGION\t1\t1\t{seq_len}\n")

        for pos in range(seq_len):
            gt = MSA_mtx[pos]
            if gt[0]!="N" and np.all(gt == gt[0]):
                tally[0] += 1
                continue
            if np.all(gt == "N"):
                tally[1] += 1
            else:                
                if not np.all(np.isin(gt, bases)):
                    gt[np.logical_not(np.isin(gt, bases))] = "N"
                    tally[2] += 1
                else:
                    tally[3] += 1
                # if np.sum(np.isin(np.unique(gt), ['A', 'G', 'C', 'T'])) < 2:
                #     tally[4] += 1
                #     continue

            sitesF.write(f"{pos+1}\t{''.join(gt)}\n")

    print(f"Finished writing SITES file: Invariant-{tally[0]}, Unk-{tally[1]}, Ambiguous/Gap-{tally[2]}, Variant-{tally[3]}")
#    print(f"Site omitted (non-biallelic): {tally[4]}")

    return 0

sys.exit(main(sys.argv))