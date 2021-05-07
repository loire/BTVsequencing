import os,sys
f= os.listdir(os.getcwd())
files=[]
for i in f:
    if i[-3:]=="tsv":
        files.append(i)        
print(files)
header = open(files[0],'r').readline()
#print(header)
fout = open("merged_vcf.tsv",'w')
fout.write("SAMPLE\t"+header)
for f in files:
    if not f[-3:]=="tsv":
        continue
    sample = f.split(".")[0]
    fin = open(f,'r')
    fin.readline()
    for line in fin:
        fout.write(sample+"\t"+line)
    fin.close()
fout.close()




        
