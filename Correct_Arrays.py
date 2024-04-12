import numpy as np
import pandas as pd
import argparse
import os, sys
from collections import Counter
parser = argparse.ArgumentParser(
                    prog='Correct_arrays',
                    description="Correct arrays' SNP position and Allele cding format ",
                    epilog='This pipeline includes 1. Check arrays version 2. Correct the position 3. Check the Allele cding format \
                    4. Convert to consistent version and format. Also, two functions is optional “Compare the difference of ref and alt allele after flip among versions” and “Compare the difference of SNP position among versions”.')
parser.add_argument('-w','--setWS',help="set work space default ./", default="./")
parser.add_argument('-s','--SNPchimp',help="SNPchimp SNP information file including the necessary version when you convert (Download from SNPchiMp website)")
parser.add_argument('-f', '--bfile',help="bim file prefix, example 'dryad12345'",default="")

parser.add_argument('-v', '--version',help="if you known the version of this array, please input the standard name as same as the SNPchimp file, default is no ",default="no")
parser.add_argument('--flag_99',help="change 99 to MT, default is Y",default="Y")
parser.add_argument('--inputF',help="Allele coding format of input,  default is no (means unknown, programnnes will search chip with the most consistent variants as the input Format).\
                                  option=Alleles_A_B_FORWARD, Alleles_A_B_TOP,Alleles_A_B_Affymetrix",default="no")
parser.add_argument('--outputF',help="Expected Allele coding format, default is Alleles_A_B_FORWARD. Options = Alleles_A_B_TOP,Alleles_A_B_Affymetrix",default="Alleles_A_B_FORWARD")
parser.add_argument('--funtion',help='Correct_postion:1, Correct allele coding format:2, Both 3, default is 3\n \
 Compare the difference of ref and alt allele after flip of a array version:11, Compare the difference of SNP position and  difference of ref and alt allele after flip between versions: 12',default=3,type=int)

parser.add_argument('-v1', '--version1',help="the version of this arrays you expect to compare, neccesary for function 11", default="no")
parser.add_argument('-v2', '--version2',help="the version of this arrays you expect to compare, neccesary for function 12", default="no")

args = parser.parse_args()

# Check which version the chip is
# he meta information of some publish data is not complete, for example show the chip name but no chip version.
# based on the SNP-name
class Correction:
    def __int__(self,SNPchimp, bfile, log):
        self.SNPchimp=SNPchimp
        self.bfile=bfile
        self.bim = pd.read_csv(f'{self.bfile}.bim', sep="\t", header=None)
        self.log=log
    def check_version(self,uni_chip):
        bim_SNP = self.bim[1]
        for i in uni_chip:
            chip_set = self.SNPchimp.loc[self.SNPchimp["chip_name"] == i, "SNP_name"]
            # Check if all values of bim_SNP are included in chip_set
            if set(bim_SNP).issubset(set(chip_set)):
                chip_name = i
                self.log.write(f"All SNPs in {self.bfile} are included in {i}\n")
        self.log.write(f"Regarding the version of {self.bfile} as {chip_name}\n")
        return chip_name

    def Correct_position(self,chip_name,flag):
        df1 = self.SNPchimp.loc[self.SNPchimp["chip_name"] == chip_name, :]
        if flag == "Y":
            df1.loc[df1["chromosome"] == "99", "chromosome"] = "MT"
            df1.loc[df1["chromosome"] == 99, "chromosome"] = "MT"
            self.log.write(f"Because plink couldn't handle '99', change 99 to MT\n")

        bim_new = pd.merge(self.bim, df1, left_on=1, right_on="SNP_name", how="left")[["chromosome", 1, 2, "position", 4, 5]]
        bim_new.to_csv(f'{self.bfile}_positon.bim', header=False, index=False, sep="\t")
        self.bim=bim_new
        self.log.write(f"output corrected postion .bim file\n")

    def fun_check_allele(self,x, column_name):
        if f"{x[4]}/{x[5]}" == x[column_name]:
            return 1
        elif f"{x[5]}/{x[4]}" == x[column_name]:
            return 1
        else:
            return 0

    def check_format(self,chip_name):
        df1 = self.SNPchimp.loc[self.SNPchimp["chip_name"] == chip_name, :]
        bim_df1 = pd.merge(self.bim, df1, left_on=1, right_on="SNP_name", how="left")

        bim_df1["compare_Alleles_A_B_FORWARD"] = bim_df1.apply(lambda x: self.fun_check_allele(x, "Alleles_A_B_FORWARD"),
                                                               axis=1)
        bim_df1["compare_Alleles_A_B_TOP"] = bim_df1.apply(lambda x: self.fun_check_allele(x, "Alleles_A_B_TOP"), axis=1)
        bim_df1["compare_Alleles_A_B_Affymetrix"] = bim_df1.apply(
            lambda x: self.fun_check_allele(x, "Alleles_A_B_Affymetrix"), axis=1)
        dic_number = {"Alleles_A_B_FORWARD": bim_df1["compare_Alleles_A_B_FORWARD"].sum(),
                      "Alleles_A_B_TOP": bim_df1["compare_Alleles_A_B_TOP"].sum(),
                      "Alleles_A_B_Affymetrix": bim_df1["compare_Alleles_A_B_Affymetrix"].sum()}
        for i, j in dic_number.items():
            self.log.write(self.bfile, f": the number of SNPs {i} consistent with format:  {j}")
        return max(dic_number.keys(), key=lambda x: x[1])

    def fun_flip_allele(self, x, Input, Output, missing=["0", "-"]):
        # code by Paolo
        Input = x[Input].split("/")
        Output = x[Output].split("/")

        result = []

        for allele in [x[4], x[5]]:

            if allele in missing:
                result.append("0")
            elif allele not in Input:
                result.append("0")
                self.log.write(f"{x['SNP_name']}: Allele {allele} not in {Input}")
            else:
                result.append(Output[Input.index(allele)])

        x[4], x[5] = result

        return x

    def convert(self,chip_name, origin_format, output_format):
        if origin_format == output_format:
            pass
        else:
            df1 = self.SNPchimp.loc[self.SNPchimp["chip_name"] == chip_name, :]
            bim_df1 = pd.merge(self.bim, df1, left_on=1, right_on="SNP_name", how="left")
            bim_df1 = bim_df1.apply(lambda x: self.fun_flip_allele(x, origin_format, output_format), axis=1)
            bim_new = bim_df1[[0, 1, 2, 3, 4, 5]]
            self.log.write("# After convert")
            self.check_format(chip_name)
            bim_new.to_csv(f'{self.bfile}_converted.bim', header=False, index=False, sep="\t")

class Compare:
    def __init__(self,SNPchimp):
        self.SNPchimp = SNPchimp
    def fun_check(self,x, c1, c2):
        dic = {"A": "T", "T": "A", "C": "G", "G": "C"}
        a, b = x[c1].split("/")
        if f"{a}/{b}" == x[c2]:
            return 1
        elif f"{b}/{a}" == x[c2]:
            return 1
        elif f"{dic[a]}/{dic[b]}" == x[c2]:
            return 2
        elif f"{dic[b]}/{dic[a]}" == x[c2]:
            return 2
        else:
            return 0
    #This step is for checking whether the same version's SNP ref and alt allele is consistent.
    def InVersion(self,version):
        if version=="no":
            df1=self.SNPchimp
        else:
            df1 = self.SNPchimp.loc[self.SNPchimp["chip_name"] ==version, :]
        c=Counter(df1.apply(lambda x: self.fun_check(x, "Alleles_A_B_FORWARD","Alleles_A_B_TOP"),axis=1))
        print("The difference between Alleles_A_B_FORWARD and Alleles_A_B_TOP")
        print(f"completely consistent: {c[1]}")
        print(f"consistent after flip: {c[2]}")
        print(f"others: {c[0]}")
        c=Counter(df1.apply(lambda x: self.fun_check(x, "Alleles_A_B_FORWARD", "Alleles_A_B_Affymetrix"), axis=1))
        print("The difference between Alleles_A_B_FORWARD and Alleles_A_B_TOP")
        print(f"completely consistent: {c[1]}")
        print(f"consistent after flip: {c[2]}")
        print(f"others: {c[0]}")
        c=Counter(df1.apply(lambda x: self.SNPchimp.fun_check(x, "Alleles_A_B_TOP", "Alleles_A_B_Affymetrix"), axis=1))
        print("The difference between Alleles_A_B_FORWARD and Alleles_A_B_TOP")
        print(f"completely consistent: {c[1]}")
        print(f"consistent after flip: {c[2]}")
        print(f"others: {c[0]}")

    def postion_compare(self,x):
        if x["chromosome_x"] == x["chromosome_y"]:
            if x["position_x"] == x["position_y"]:
                return 1
            else:
                return 0
        else:
            return 0
    def AmongVersion(self,version1,version2):
        li_f=["Alleles_A_B_FORWARD","Alleles_A_B_TOP","Alleles_A_B_Affymetrix"]
        df_1 = self.SNPchimp.loc[self.SNPchimp["chip_name"] == version1, ["rs","SNP_name","chromosome", "position","Alleles_A_B_FORWARD","Alleles_A_B_TOP","Alleles_A_B_Affymetrix"]]
        df_2 = self.SNPchimp.loc[self.SNPchimp["chip_name"] == version2, ["rs","SNP_name","chromosome", "position","Alleles_A_B_FORWARD","Alleles_A_B_TOP","Alleles_A_B_Affymetrix"]]

        df1=df_1[df_1["rs"].notna()]
        df2 = df_2[df_2["rs"].notna()]
        print(f"The number of not Nan value in rs of {version1} :{df1.shape[0]}; total SNP number: {df_1.shape[0]}")
        print(f"The number of not Nan value in rs of {version2} :{df2.shape[0]}; total SNP number: {df_2.shape[0]}")
        a=pd.merge(df1,df2,how="outer",on="rs")
        df1=df_1[df_1["SNP_name"].notna()]
        df2 = df_2[df_2["SNP_name"].notna()]
        print(f"The number of not Nan value in SNP name of {version1} :{df1.shape[0]}; total SNP number: {df_1.shape[0]}")
        print(f"The number of not Nan value in SNP name of {version2} :{df2.shape[0]}; total SNP number: {df_2.shape[0]}")
        b=pd.merge(df1,df2,how="outer",on="SNP_name")
        for i in range(len(li_f)):
            for j in range(len(li_f)):
                if (j==2) or (i==2):
                    df_deal=a
                else:
                    df_deal=b
                print(f"The difference of SNP position between {version1} and {version2}")
                c=Counter(df_deal.apply(self.postion_compare,axis=1))
                print(f"completely consistent: {c[1]}")
                print(f"others : {c[0]}")
                print(f"The difference of allele coding between {i} of {version1} and {j} of {version2}")
                c = Counter(
                        df_deal.apply(lambda x: self.fun_check(x, f"{i}_x", f"{j}_y"), axis=1))
                print(f"completely consistent: {c[1]}")
                print(f"consistent after flip: {c[2]}")
                print(f"others: {c[0]}")


os.chdir(args.setWS)
SNPchimp=pd.read_csv(args.SNPchimp)
uni_chip=SNPchimp["chip_name"].unique()
if args.funtion==11:
    C=Compare(SNPchimp)
    C.InVersion(args.version1)

elif args.funtion==12:
    C = Compare(SNPchimp)
    C.AmongVersion(args.version1)

elif args.funtion<=3:
    log = open(sys.bfile + '.log', 'w')
    P = Correction(SNPchimp, args.bfile)
    chip_name = args.version
    if chip_name == "no":
        if (args.funtion==1) or (args.funtion==3):
            log.write(f"# Check_version\n")
            chip_name=P.check_version(uni_chip)
        else:
            log.write(f"Error: function 2 must input array version\n")
            print(f"Error: function 2 must input array version\n")
            sys.exit()
    else:
        if chip_name in uni_chip:
            log.write(f"the version of {args.bfile} is {chip_name}\n")
        else:
            log.write(f"Error: the version of {args.bfile} doesn't include {chip_name}\n")
            print(f"Error: the version of {args.bfile} doesn't include {chip_name}\n")
            sys.exit()
    if (args.funtion == 1) or (args.funtion == 3):
        log.write(f"# Correct_position\n")
        P.Correct_position(chip_name,args.flag_99)

    if (args.funtion==2) or (args.funtion==3):
        Allele_format = args.inputF
        if Allele_format == "no":
            origin_format=P.check_format(chip_name)
        else:
            if Allele_format in ["Alleles_A_B_FORWARD", "Alleles_A_B_TOP","Alleles_A_B_Affymetrix"]:
                origin_format=Allele_format
            else:
                log.write(f"Error: Allele_format input is not in three format \n")
                print(f"Error: the version of {args.bfile} doesn't include {chip_name}\n")
                sys.exit()
        P.convert(chip_name, origin_format, args.outputF)
    log.close()
else:
    print(f"The function is {args.function}, is not existing")

