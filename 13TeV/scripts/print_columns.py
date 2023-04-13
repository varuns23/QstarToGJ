a=[]
for line in open("xs_f0p1.txt"):
       col = line.split()
       a.append(col[2])

print ', '.join(a) 


