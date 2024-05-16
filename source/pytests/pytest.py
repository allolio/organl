import pynagata                   
import numpy as np                
p=pynagata.FaceVertexMesh()       
p.LoadObj('testmesh/uvsphere.obj')
n=pynagata.Nagata()               
n.FromMesh(p)               
n.PrepareCoefficientVector()
a=pynagata.DoubleVector();
n.GetCoefficientVector(a);
print(a)
