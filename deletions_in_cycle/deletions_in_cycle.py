#!/usr/bin/python

import fileinput
import re
import matplotlib.pyplot as plt
import numpy as np
from operator import truediv
def get_del_pos (cigar):
    pos=0
    pozs = []
    for cig in cigar:
        if cig[1] == 'M' or cig[1] == 'I' or cig[1] == 'S':
            pos+=int(cig[0])
        elif cig[1] == 'D':
            pozs.append(pos)
        else:
            print "another flag"
            pozs.append(-1)
    return pozs


del_positions=[]
cum_readlengths=[]
read_lengths=[]
cigar="3M2S14I2D5M"
for line in fileinput.input():
    fields=line.split("\t")
    cigar =re.findall(r'(\d+)(\w)', fields[5])
    del_positions+=get_del_pos(cigar)
    read_lengths.append(len(fields[9]))


plt.hist(del_positions,bins=300)
values, base = np.histogram(read_lengths, bins=300)

del_values, del_base = np.histogram(del_positions,bins=300)
#values*=np.max(del_values)/float(np.max(np.cumsum(values)))
#del_values *= np.max(del_values)/float(np.max(np.cumsum(del_values)))
#plt.plot(del_base[:-1],np.cumsum(del_values), c='green',linewidth=4)
plt.plot(base[:-1],(len(read_lengths) - np.cumsum(values))-np.min(len(read_lengths) - np.cumsum(values)), c='red')
plt.show()
plt.plot(del_base[:-1],map ( truediv, del_values,(len(read_lengths) - np.cumsum(values)) ) ,c='black',linewidth=5)
plt.show()
