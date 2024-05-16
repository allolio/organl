
#define MAX_NONL 0.994

double LinearEdge(NTriangle *nt)
{
    double penal=0, lval=0;
    // Edge 1 Unit tangent
    triple t1, t2, t3;
    triple g1,g2;
    triple d;
   // d=*nt->v[1]-*nt->v[0];
    t1=*nt->n[0]%nt->x_eta(0,0); t2=*nt->n[1]%nt->x_eta(1,0); t3=nt->normal(0.5,0)%nt->x_eta(0.5,0);
    t1=t1/t1.abs();t2=t2/t2.abs();t3=t3/t3.abs();
//    cout << "TAN"  << t1*t2 << "  " << endl;
//    cout << "NCHECK" << *nt->n[0]*nt->normal(0,0) << endl;
    t1=t1+t2;t1=t1/t1.abs();
    lval=t1*t3;
//    cout << "WTF A" << t1*t3  << endl;
 if (lval< MAX_NONL) penal+=(lval-1)*(lval-1);
        d=*nt->v[0]-*nt->v[2]; // Check
        g1=d-nt->cbase[2]-nt->cbase[5];
        g2=d+nt->cbase[2]+nt->cbase[5]*2;
//      t1=*nt->n[0]%nt->x_eta(0,0); t2=*nt->n[1]%nt->x_eta(1,0); t3=nt->normal(0.5,0)%d;
            t1=*nt->n[0]%g2; t2=*nt->n[2]%g1; t3=nt->normal(0.5,0.5)%(d-nt->cbase[5]*0.25);
    t1=t1/t1.abs();t2=t2/t2.abs();t3=t3/t3.abs();

        lval=t1*t3;
//    cout << "WTF A" << t1*t3  << endl;
 if (lval< MAX_NONL) penal+=(lval-1)*(lval-1);

    
  //  cout << "WTF B" << t1*t3  << endl;    
    d=*nt->v[2]-*nt->v[1];
         g1=d-nt->cbase[1]-nt->cbase[4];
         g2=d+nt->cbase[1]+nt->cbase[4]*2;
    t1=*nt->n[1]%g1; t2=*nt->n[2]%g2; t3=nt->normal(1.0,0.5)%(d-nt->cbase[4]*0.25);
       t1=t1/t1.abs();t2=t2/t2.abs();t3=t3/t3.abs();
  //      cout << "WTF C" << t1*t3  << endl;
        lval=t1*t3;
//    cout << "WTF A" << t1*t3  << endl;
 if (lval< MAX_NONL) penal+=(lval-1)*(lval-1);
       
    return penal;
    
}


double CheckEdgeNormals(Property *e) // E is EdgeProperty
{
    if(e->deps.size()<2) return 1.0; // Issue is irrelevant for dangling edges.
    if(e->deps[1]->ptype->collidable==false) return 1.0; // Dangling edge with boundary energy.
    FaceEnerProperty *f1;
    NTriangle *nt1=(NTriangle*) (e->deps[0])->content;
    NTriangle *nt2=(NTriangle*) (e->deps[1])->content;
    
    // Find common vertices
    int in1[2],in2[2],k=0,l=0;
    triple norm1,norm2;
    for(int i=0;i!=3;i++) 
        for(int j=0;j!=3;j++)
            if(nt1->v[i]==nt2->v[j])
            {
                in1[k]=i;in2[l]=j;
                k++;l++;if(k==2) break;
            }
     if(in1[0]==0)
     if(in1[1]==1) norm1=nt1->normal(0.5,0);
     else norm1=nt1->normal(1.0,0.5);
    else norm1=nt1->normal(0.5,0.5); 

    if(in2[0]==0)
     if(in2[1]==1) norm2=nt2->normal(0.5,0);
     else norm2=nt2->normal(1.0,0.5);
    else norm2=nt2->normal(0.5,0.5); 
    // cout << norm1*norm2 << endl;
  //  if(norm1*norm2 < 0.95 ) return false;
    return norm1*norm2;
}


double G1Diagnostic(Membrane *mmb)
{

    PropertyList *pl; string e="Edge"; double retval=0;
   mmb->props.GetPropertyList(e ,&pl);
 //  cout << "GotEdges" << endl;
   for(int i=0;i!=pl->theProps.size();i++)
       retval+=CheckEdgeNormals(pl->theProps[i]);
   return retval/(double) pl->theProps.size();
}  

bool CheckFaceQG1Property(NTriangle *, FaceVertexMesh *fvm, int index0) // Index of Face
{

// This requires index0 to be the first vertex in the Face
    // First we need to find the neighboring Triangles. This can be done via edges.
    // Edges--> Triangles
    // Compare By edge. We will compare three edges.
    // Triangle --> Edges of Triangle
    // So we need to get --> Triangles of Edge.
    // Then we evaluate Value at Edge. So we need orientation
    // then we return whether G1 property is present.
    
// What normal at what vertex:
  /*     inline triple eval(double eta,double xi)
    {
            return (*v[0]*(1-eta))+(*v[1]*(eta-xi))+*v[2]*xi-c[0]*(1-eta)
               *(eta-xi)-c[1]*((eta-xi)*xi)-c[2]*((1-eta)*xi);
    }
    */
  // [0]-[1] :  0.5 0  
  // [0]-[2] :  0.5 0.5
  // [1]-[2]:   1.0 0.5

  return true;
}
