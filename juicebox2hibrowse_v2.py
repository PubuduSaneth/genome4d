
# coding: utf-8

# # Read juicebox dump output and reformat to 7 column format

# ## Juicebox dump output format - Input of the program
# <pre>
# chr1:start chr2:start Normalized_interactions
# 10000	10000	311.05484
# 10000	20000	92.60087
# 20000	20000	296.0056
# 10000	30000	47.701942
# </pre>

# ## 7 column format - Output of the program
# <pre>
# chr1	start1	end1	chr2	start2	end2	Normalized_interactions
# 1	10000	20000	1	10000	20000	311.05484
# 1	10000	20000	1	20000	30000	92.60087
# 1	20000	30000	1	20000	30000	296.0056
# 1	10000	20000	1	30000	40000	47.701942
# </pre>

# ### runDump.bash: script that runs juicebox dump command
#
# #### Jucebox dump parameters:
# * Resolution: 100K
# * Normalization method: KR (normalization method used in the rao paper)
# * Output: observed notmalized values (in a sparce matrix)
#
# #### Bash script
# <pre>
# #!/bin/bash
# ## $1 ==  chr1
# ## $2 ==  chr2
# ## $3 ==  output_file
# echo -e "run juicebox dump\n\t resolution: 100000\n\t normalization method: KR\n\t output: observed notmalized sparce matrix"
# touch $3
#
# java -jar /Users/pubudu/Downloads/MacTools/juicebox_tools.7.5.jar dump observed KR /Users/pubudu/Documents/RefData/raoPaper/copy_GSE63525_HMEC_combined.hic $1 $2 BP 100000 $3
# echo -e "\#\#\#\#juicebox dump observed KR copy_GSE63525_HMEC_combined.hic $1 $2 BP 100000 $3\n$(cat $3)" > $3
# </pre>

import pandas as pd
import subprocess
from itertools import combinations

chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y'] # List of chromosomes
f_counter = 0 # Counter for the number of files created
set_path = "/projects/rrresearch/Pubudu/hic/datasets/rao/GSE63525_HMEC/hic_file/" # path of the .hic file and dump files
file_list = [] # list of all the dump files created
res = 100000 # Resolution used to run juicebox dump command

def run_Dump(arg_list):
    subprocess.call(['/projects/rrresearch/Pubudu/hic/tools/runDump.bash', str(arg_list[0]), str(arg_list[1]), ''.join([set_path, str(arg_list[2])]) ])

# Run juicebox dump to
for chroms in chr_list:
    print 'Juicebox dump run - chromosome: {}'.format(chroms)
    f_name = '{}.chr{}_chr{}.dumpOut.txt'.format(f_counter, chroms, chroms)
    arg_list = [chroms, chroms, f_name]
    file_list.append(f_name)
    run_Dump(arg_list)
    f_counter += 1

# itertools.combinations(iterable, r)
### Return r length subsequences of elements from the input iterable.
### Combinations are emitted in lexicographic sort order. So, if the input iterable is sorted, the combination tuples will be produced in sorted order.
### Elements are treated as unique based on their position, not on their value. So if the input elements are unique, there will be no repeat values in each combination.

for combo in combinations(chr_list, 2):
    print 'Juicebox dump run - chromosomes: {} and {}'.format(combo[0], combo[1])
    f_name = '{}.chr{}_chr{}.dumpOut.txt'.format(f_counter, combo[0], combo[1])
    arg_list = [combo[0], combo[1], f_name]
    file_list.append(f_name)
    run_Dump(arg_list)
    f_counter += 1

#print file_list
print len(file_list)

for f_name in file_list:
    chr_pair = f_name.split('.')[1].split('_')
    print 'Reformatting to 7 columns: {} - chr1 {} and chr2 {} '.format(''.join([set_path,f_name]), chr_pair[0], chr_pair[1])
    file = open('{}.hibrowse.txt'.format(''.join([set_path,f_name])),'w')
    with open(''.join([set_path,f_name])) as fx:
        next(fx)
        for line in fx:
            line = line.replace('\n','').replace('\r','').split('\t')
            if line[2] == 'NaN':
                intCounts = 0
            else:
                intCounts = float(line[2])
            file.write('{}\n'.format('\t'.join([chr_pair[0], line[0], str(int(line[0])+100000), chr_pair[1], line[1], str(int(line[1])+100000), str(intCounts)])))
    file.close()
