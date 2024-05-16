
#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

#print('Number of arguments: {}'.format(len(sys.argv)))
#print('Argument(s) passed: {}'.format(str(sys.argv)))


# In[25]:


#sys.argv[1]="../testmesh/uvsphere.obj"
#print(sys.argv[1])
control=sys.argv[2]
file1 = open(sys.argv[1], 'r') 
Lines = file1.readlines() 
count = 0
fcount = 0
val=float(sys.argv[3])
statement ='True'
if len(sys.argv) > 4 :
    statement=sys.argv[4]
#statement= 'z > 10 or z < 0'
# Strips the newline character 
select=list()
for line in Lines: 
    a=line.split()
    if(len(a)==0):
        a="X"
    if(a[0]=="v"):
            x=float(a[1])
            y=float(a[2])
            z=float(a[3])
            if(eval(statement)):
                select.append(count+1)
            count+=1
    if(a[0]=="f"):
        if int(a[1]) in select or int(a[2]) in select or int(a[3]) in select:
            print("%i %s %f" % (fcount,control,val))
        fcount+=1


# In[ ]:

