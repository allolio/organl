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
statement ='True'
if (control=="c"):
    statement=sys.argv[3]
#statement= 'z > 10 or z < 0'
# Strips the newline character 
for line in Lines: 
    a=line.split()
    if(len(a)==0):
        a="X"
    if(a[0]=="v"):
        if(control=="n" or control=="a"):
            print ("Normal %i 0" % count)
            print ("Normal %i 1" %  count)
        if(control=="a" or control=="v"):
            print ("Vertex %i 0" % count)
            print ("Vertex %i 1" %  count)
            print ("Vertex %i 2" %  count)
        if(control=="c"):
            x=float(a[1])
            y=float(a[2])
            z=float(a[3])
            if(eval(statement)):
                print ("Vertex %i 0" % count)
                print ("Vertex %i 1" %  count)
                print ("Vertex %i 2" %  count)                
                print ("Normal %i 0" % count)
                print ("Normal %i 1" %  count)
        count+=1


# In[ ]:




