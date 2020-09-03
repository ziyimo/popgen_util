#!/usr/bin/env python3

import sys, os
import numpy as np

helpMsg = '''
        usage: $./make_coal.py <popsize_file> <outPref>
'''

# GBR #
# dem = "-en 0.000120 0 0.124319 -en 0.000272 0 0.042569 -en 0.000399 0 0.031529 -en 0.000532 0 0.023182 -en 0.000665 0 0.017045 -en 0.000797 0 0.012532 -en 0.000930 0 0.009214 -en 0.001063 0 0.006576 -en 0.001224 0 0.009894 -en 0.001329 0 0.009894 -en 0.001595 0 0.009894 -en 0.001994 0 0.009910 -en 0.002713 0 0.076953 -en 0.003722 0 0.076953 -en 0.004918 0 0.076953 -en 0.006247 0 0.076892 -en 0.007870 0 0.038865 -en 0.008507 0 0.038865 -en 0.009304 0 0.038865 -en 0.010367 0 0.038865 -en 0.011962 0 0.038865 -en 0.014621 0 0.038865 -en 0.018608 0 0.038865 -en 0.023925 0 0.038865 -en 0.033229 0 0.038865 -en 0.046521 0 0.038865 -en 0.066458 0 0.038865 -en 0.132917 0 0.038865 -en 0.398750 0 0.038865"
# NE = 188088

# bird #
dem = "-en 0.0743 0 97.16 -en 3.125 0 9.79"
NE = 148000

def discoal_dem2popsizeF(discoal_dem_str, file_name, Ne): # diploid Ne
    dem_arg_ls = discoal_dem_str.split()
    with open(file_name, "w") as popsizeF:
        print(f"0\t{Ne}", file=popsizeF)
        for idx in range(0, len(dem_arg_ls), 4):
            print(f"{round(float(dem_arg_ls[idx+1])*4*Ne)}\t{round(float(dem_arg_ls[idx+3])*Ne)}", file=popsizeF)

def main(args):
    if len(args) != 3:    #2 argument(s)
        return helpMsg

    popF = args[1]
    out_pref = args[2]

    if not os.path.isfile(popF):
        discoal_dem2popsizeF(dem, popF, NE)

    gen_popsize = np.loadtxt(popF, dtype=int)

    with open(out_pref+".coal", "w") as coalF:
        coalF.write("0\n")
        coalF.write(" ".join(gen_popsize[:, 0].astype(str))+"\n")
        coalF.write("0 0 "+" ".join((0.5/gen_popsize[:, 1]).astype(str))+"\n")

    return 0

sys.exit(main(sys.argv))
