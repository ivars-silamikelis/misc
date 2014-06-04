# -*- coding: utf-8 -*-
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import re
import glob


def extract_data(csv_file):
	data={"MISSENSE":[],"NONSENSE":[],"SILENT":[],"Transitions":[],"Transversions":[],"SNP":[],"Ts_Tv":[],"dNdS":[]}
	with open(csv_file) as content:
		for line in content:
			line=re.sub(" ","",line)
			fields=line.split(",")
			if fields[0] in data.keys():
				data[fields[0]]=int(fields[1])
	print "Transversions "+str(data["Transversions"])+str(data["Transversions"]==0)
	if data["Transversions"]==0:
		data["Ts_Tv"]=0
	else:
		data["Ts_Tv"]=data["Transitions"]/float(data["Transversions"])
	print "SILENT "+str(data["SILENT"])
	if data["SILENT"]==[]:
		data["dNdS"]=0
	else:
		print data["SILENT"]
		data["dNdS"]=data["MISSENSE"]/float(data["SILENT"])
	return data
collected_data={"MISSENSE":[],"NONSENSE":[],"SILENT":[],"Transitions":[],"Transversions":[],"SNP":[],"Ts_Tv":[],"dNdS":[],"samples":[]}
samples=glob.glob("S*.csv")
for sample in samples:
	name=sample.replace(".csv","")
	print name
	collected_data["samples"].append(name)
	dati=extract_data(sample)
	for key in dati.keys():
		collected_data[key].append(dati[key])
bin_width=0.35
index=np.arange(0,len(collected_data["samples"]))*bin_width*(len(collected_data["samples"]))
fig = plt.figure(figsize=(16, 16), dpi=200)
ax = fig.add_subplot(111)
plt.bar(index+bin_width,collected_data["SNP"],bin_width,label="Kopā".decode("utf8"),color="grey")
plt.bar(index+bin_width*2,collected_data["MISSENSE"],bin_width,label="Nesinonīmie SNP".decode("utf8"),color="#09cc0a")
plt.bar(index+bin_width*3,collected_data["SILENT"],bin_width,label="Sinonīmie SNP".decode("utf8"),color="#ab04dd")
plt.bar(index+bin_width*4,collected_data["NONSENSE"],bin_width,label="Stop SNP".decode("utf8"),color="#f6ff00")
plt.title("Ar Samtools un GATK noteikto SNP\nraksturlielumi Latvijā sastopamajiem izolātiem".decode("utf8"),fontname="Times New Roman",fontweight="bold",fontsize=32)
plt.xlabel("Izolāts".decode("utf8"),fontname="Times New Roman",fontsize=28)
plt.ylabel("SNP skaits",fontname="Times New Roman",fontsize=28 )
plt.xticks(index+bin_width*(len(collected_data["samples"])/2.0),collected_data["samples"],fontname="Times New Roman",fontsize=24)
plt.yticks(fontname="Times New Roman",fontsize=24)
for i in range(0,len(index)):
	larger_key="SILENT"
	if collected_data["MISSENSE"]>collected_data["SILENT"]:
		larger_key="MISSENSE"

	ax.annotate("dN/dS:"+"{:.2f}".format(collected_data["dNdS"][i]),fontweight="bold",xy=((index+bin_width*2)[i],collected_data[larger_key][i]+2),fontname="Times New Roman",fontsize=22)
	ax.annotate("Ts/Tv:"+"{:.2f}".format(collected_data["Ts_Tv"][i]),fontweight="bold",xy=((index+bin_width)[i],collected_data["SNP"][i]+2),fontname="Times New Roman",fontsize=22)
#ax.legend()
plt.legend(prop={"size":22},loc="upper center", ncol=4,fancybox=True, shadow=True)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.19)

plt.savefig("sample_stats.png",bbox_inches="tight")
#plt.show()
