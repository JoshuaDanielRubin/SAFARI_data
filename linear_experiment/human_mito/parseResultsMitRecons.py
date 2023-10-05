#!/usr/bin/python

import sys;
from tqdm import tqdm


MAXVAL=600
STEPS10=range(10, MAXVAL, 10)
STEPSMSA=["0","10","20","30","40","50","100","150","200","300","400","500","600"]
LENGTHS=["30","35","40","45","50","55","60","65","75","85","100","150","200","400"]
RATES=["0.00010","0.00015","0.00020","0.00025","0.0003","0.00035","0.0004","0.00045","0.0005","0.0006","0.0007","0.0008","0.0009","0.001","0.003","0.005","0.007","0.009","0.01","0.03","0.05","0.07","0.09","0.1","0.3","0.5","0.7","0.9"]
DAMAGE=["dhigh","dmid","single","none"]
ALIGNERS=["mem","aln","aln_anc","bowtie2","shrimp","vg","bb"]
ALIGNERSC=["aln_anc","bowtie2","shrimp"]
CONSENSUS=["mpileup"]
NUMFRAGS=1000000
FRAGNUM=175000

filestoopen=[];

strtoprint=[];

for step in STEPSMSA:
    for fraglen in LENGTHS:
        for dam in DAMAGE:
            for rate in RATES:
                for align in ALIGNERS:
                            strtoprint.append(""+str(step)+"\t"+str(NUMFRAGS)+"\t"+str(fraglen)+"\t"+str(dam)+"\t"+str(rate)+"\t"+str(align)+"\t");
                            filestoopen.append("/home/projects/benchmarkmito/alignments/numtS_and_gen_"+str(step)+"_n"+str(NUMFRAGS)+"_l"+str(fraglen)+"_d"+str(dam)+"_s"+str(rate)+"_"+str(align)+".stat");


for i in tqdm(range(len(filestoopen))):
#    print(i);
    fin=open(filestoopen[i],"r");
    fieldstoprint=[];
    for line in fin.readlines():
    #fout=open(sys.argv[2],w);

        line=line.rstrip()
        #print(line);
        fields=line.split();
        #print(len(fields));
        fieldstoprint.append( fields[1] );
#    fout.writelines(line+"\n");

    fin.close();
    print(strtoprint[i]+"\t"+("\t".join(fieldstoprint)));
#fout.close();

