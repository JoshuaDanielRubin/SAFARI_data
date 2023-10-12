#!/usr/bin/python

import pysam
import sys; 

pysam.set_verbosity(0)
bamInputFile  = pysam.Samfile(sys.argv[1], "rb");
#bamOutputFile = pysam.Samfile(sys.argv[2], "wb", template=bamInputFile);

mapped=0;
unmapped=0;
correctlocationmapped=0;
numt=0;
total=0;
def intersects(left1,right1,left2,right2):
    return not (right1 < left2 or left1 > right2)

for read in bamInputFile:
    if(read.is_unmapped):
        unmapped+=1;
        continue;

    
    if(read.query_name[0:10] == "generation"):
        mapped+=1;
        fields=read.query_name.split(":");
        #print(fields);
        #print(read.reference_start)
        #print(read.reference_start+read.reference_length);
        intcs=intersects(int(fields[2]),int(fields[3]),read.reference_start-50,read.reference_start+read.reference_length+50);
        #print(intcs);
        if(intcs):
            correctlocationmapped+=1;
            
    else:
        numt+=1;
    total+=1;
    
    #bamOutputFile.write(read);
                
bamInputFile.close();
#bamOutputFile.close();


print("Total:\t"+str(total));
print("unmapped:\t"+str(unmapped));
print("unmapped%:\t"+str(100*(unmapped/total)));
print("mapped:\t"+str(mapped));
print("mapped%:\t"+str(100*(mapped/total)));
print("correctmap:\t"+str(correctlocationmapped));
print("correctmap%:\t"+str(100*(correctlocationmapped/mapped)));
print("numt:\t"+str(numt));
print("numt%:\t"+str(100*(numt/total)));




