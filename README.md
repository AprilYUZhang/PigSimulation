Correct_Arrays.py 

packge: 
numpy                     1.24.3 
pandas                    2.0.3 

usage: Correct_arrays [-h] [-w SETWS] [-s SNPCHIMP] [-f BFILE] [-v VERSION]
                      [--flag_99 FLAG_99] [--inputF INPUTF]
                      [--outputF OUTPUTF] [--funtion FUNTION] [-v1 VERSION1]
                      [-v2 VERSION2]

Correct arrays' SNP position and Allele cding format

options:
  -h, --help            show this help message and exit
  -w SETWS, --setWS SETWS
                        set work space default ./
  -s SNPCHIMP, --SNPchimp SNPCHIMP
                        SNPchimp SNP information file including the necessary
                        version when you convert (Download from SNPchiMp
                        website)
  -f BFILE, --bfile BFILE
                        bim file prefix, example 'dryad12345'
  -v VERSION, --version VERSION
                        if you known the version of this array, please input
                        the standard name as same as the SNPchimp file,
                        default is no
  --flag_99 FLAG_99     change 99 to MT, default is Y
  --inputF INPUTF       Allele coding format of input, default is no (means
                        unknown, programnnes will search chip with the most
                        consistent variants as the input Format).
                        option=Alleles_A_B_FORWARD,
                        Alleles_A_B_TOP,Alleles_A_B_Affymetrix
  --outputF OUTPUTF     Expected Allele coding format, default is
                        Alleles_A_B_FORWARD. Options =
                        Alleles_A_B_TOP,Alleles_A_B_Affymetrix
  --funtion FUNTION     Correct_postion:1, Correct allele coding format:2,
                        Both 3, default is 3 Compare the difference of ref and
                        alt allele after flip of a array version:11, Compare
                        the difference of SNP position and difference of ref
                        and alt allele after flip between versions: 12
  -v1 VERSION1, --version1 VERSION1
                        the version of this arrays you expect to compare,
                        neccesary for function 11
  -v2 VERSION2, --version2 VERSION2
                        the version of this arrays you expect to compare,
                        neccesary for function 12

This pipeline includes 1. Check arrays version 2. Correct the position 3.
Check the Allele cding format 4. Convert to consistent version and format.
Also, two functions is optional “Compare the difference of ref and alt allele
after flip among versions” and “Compare the difference of SNP position among
versions”.
