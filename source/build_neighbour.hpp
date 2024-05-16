// These Properties are Dynamically Inserted just to Organize the Simul Face MC Move;
class FaceControlProperty : public Property
{
public:
    Deplist shallow;
    bool update(bool deep)
    {
          for(int i=0; i!=shallow.size(); i++) shallow[i]->update(false); // this will merely take care of blocks etc of original properties.
        if(deep)
            for(int i=0; i!=deps.size(); i++) deps[i]->update();
        return true;
    }
};

class FaceControlProperties : public PropertyList
{
public:
        virtual bool genpropertiesfromStructure(void *parent);
     FaceControlProperties()
    {   name="FaceC";
        szProperty=18; // Three vertex Vectors // Three Normal Vectors
        
    }
    FaceVertexMesh *mesh;
    ~FaceControlProperties()
    {free(simul);}
    double **simul;
};


bool FaceControlProperties::genpropertiesfromStructure(void* parent)
{
   // double* simul;
    this->parent=parent;            
    string prop1="Vertex";
    string prop2="Normal";
    PropertyMap *pm=(PropertyMap *) parent;
    if(!pm->HasProp(prop1)) return false;
    if(!pm->HasProp(prop2)) return false;
    PropertyList *plv,*pln;
    pm->GetPropertyList(prop1,&plv);
        pm->GetPropertyList(prop1,&pln);
    mesh=(FaceVertexMesh *) plv->parent;
    // Now iterate through mesh.
    simul=(double **) malloc(mesh->faces.size()*6*sizeof(double*));
    theProps.resize(mesh->faces.size());
    cout << theProps.size() << " " << endl;
    for(int i=0;i!=mesh->faces.size();i++)
    {
//        simul[i].resize(6);
        theProps[i]=new FaceControlProperty;
        theProps[i]->block=NULL;
        theProps[i]->ptype=this;
        theProps[i]->parent=mesh; // Need To decide on that
        theProps[i]->deps.clear();
        theProps[i]->index=i; // NOT SAFE DANGER
        // Assemble Dependents
       
        for(int j=0;j!=mesh->faces[i].vert.size();j++)
        {
            int vindex=mesh->faces[i].vert[j];
     //       cout << vindex << endl;
            // Merge Dependencies, e.g. Nagata Triangles into the Props
         //   cout << plv->theProps[vindex]->deps.size() << " " <<  pln->theProps[vindex]->deps.size() << " " << theProps[i]->deps.size();
            MergeDeplists(&(theProps[i]->deps),&(plv->theProps[vindex])->deps);
            MergeDeplists(&(theProps[i]->deps),&(pln->theProps[vindex])->deps);
            ((FaceControlProperty*) theProps[i])->shallow.push_back(plv->theProps[vindex]); // Take care of blocks and data 
          //  ((FaceControlProperty*) theProps[i])->shallow.push_back(pln->theProps[vindex]); // No Normal Pushing, its DANGER ous....
            // Now Populate Data Fields
           //  memcpy(&simul[i],plv->theProps[vindex]->content,3*sizeof(double));
       //       cout << "COPYING CONTENTS" << endl;
              simul[i*6+j]=plv->theProps[vindex]->content;
       
             // (simul[i])[j+1]=plv->theProps[vindex]->content;
          //    cout << *(simul[i])[j] << endl;
              simul[i*6+j+3]=(double*) &mesh->vertices[vindex].norm;
        }
        theProps[i]->content=(double *) &simul[i*6];
   //     for(int k=0;k!=3;k++)
  //      {triple *tp=(triple*) (((double **) theProps[i]->content)[k]);// << endl;
        //triple *tp=(triple*)  (theProps[i]->content+9);
 //       tp->print();}
    }
    return true;
};

class SimulFaceMove : public MCMove
{
public:
    SimulFaceMove()
    {
    prop="Vertex";
    tryFrac=0.5;
    maxstepsz=0.1;
    init=false;
//    vfcp.clear();
    }
    
        bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }

        int CreateParallelPlan(Simulation *sim) 
    {
        if(!init)
        {
            cout << "Injecting Control Faces" << endl;
            vfcp.resize(sim->Entities.size());
            // Now Inject Control Properties :)
            for(int i=0;i!=sim->Entities.size();i++)
            {
             if(CheckEntitySupportsMove(sim->Entities[i]))
             { vfcp[i].genpropertiesfromStructure(&(sim->Entities[i]->props));
               sim->Entities[i]->props.AddPropertyType(vfcp[i]); 
              }
            }
                        cout << "Done Injecting Control Faces" << endl;
          prop="FaceC"; init=true;
        }
       return MCMove::CreateParallelPlan(sim); // This allows us to use the default MC Planning Scheme.
    }
    
     double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
        // Backup Copy
 // cout << "EXECUTING" << endl;
#pragma omp parallel for
        
        for (int i=0;i<modprop.size();i++)
        {
            double oldall[18];             double newall[18];
//            vector <vector <double *> > *cnt;
            for(int k=0;k!=6;k++)
            {         
                  memcpy(&oldall[k*3],*(((double**) modprop[i]->content)+k),sizeof(triple));        
                //  triple *tp=(triple*) *(((double**) modprop[i]->content)+k);
                //  cout << "TRIPLE" << k << endl;
                //  tp->print();
   //               cout << oldall[k] << " " << oldall[k+1] <<" "<< oldall[k+2]  << endl;
            }
//     

            vector<double> collect(2,0);
//cout << "GOT HERE EP" << endl;
            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
  //          cout << "GOT HERE EP" << endl;
            // Normal Indices:
            FaceVertexMesh *fv=(FaceVertexMesh *) modprop[i]->parent;
            indexlist *nindex=&fv->faces[modprop[i]->index].vert;
            int retcon=0;
            retry:
           bool allowed=true;
    
            //triple center; center=center*0;
            /*
            for(int j=0;j!=modnodes[i].size();j++)
            {
            center=center+((NTriangle *) sim->n.col[modnodes[i].at(j)])->eval(0.666666,0.333333);
            }
            center=center/modnodes[i].size();
           */
            
       
            for(int k=0;k!=3;k++) { GenNewP(&oldall[k*3],&newall[k*3],3);}
            
            
            
            for(int k=0;k!=3;k++) memcpy(*(((double**) modprop[i]->content)+k) ,&newall[k*3],3*sizeof(double)); // Normals Updated through Mesh directly.
            
            //triple *nv=&((Vertex *) normal->parent)->norm;
           // Build All Normal Vectors;
             allowed=fv->BuildNormal(nindex->at(0),true) and fv->BuildNormal(nindex->at(1),true) and fv->BuildNormal(nindex->at(2),true);        
                        bool accepted=false;
//           fv->vertices[nindex->at(0)].pos.print();
 //                       fv->vertices[nindex->at(1)].pos.print();
  //          fv->vertices[nindex->at(2)].pos.print();
   //          cout<< nindex->at(0) << " " << nindex->at(1) << " " << nindex->at(2) << endl;
            
            if(!allowed)
            {
               for(int k=0;k!=6;k++)
               memcpy(*(((double**) modprop[i]->content)+k),&oldall[k*3],3*sizeof(double));             
               if(retcon<10) {retcon++; goto retry;}
            }
            else 
            {    
           // doublemod[i].first->update();
            modprop[i]->update();
            
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            //bool allowed=true;
            //if(*nv * (*nv) < 0.999) allowed=false;
            if(allowed) for(int j=0;j!=modnodes[i].size();j++)
            {
            
               allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]); //  and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], , &refn);
              if(!allowed) break;

               // sumE+=EdgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
            accepted=AcceptMove(sumE);}
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            if(!accepted)
            {
          
               for(int k=0;k!=6;k++)
               memcpy(*(((double**) modprop[i]->content)+k),&oldall[k*3],3*sizeof(double));

             // doublemod[i].first->update();

              modprop[i]->update();
             // TestReservation
           //  cout << "CUT" << endl;
           // for (int x=0;x!=cache.size();x++)  par->DeltaE(cache[x]);
            // cout << "END CUT" << endl;
            }
            else {    //    cout << "CACHE_WRITE" << endl;
                     /*     double conste=par->TotalEnergy(2); */
                      //    cout << "Vertex " << modprop[i]->index << " ";
//                           par->SumStored("Area");
  //                         par->SumStored("Volume");
                        //  for (int i=0;i!=cache.size();i++) cout << cache[i] << " " ;
                        //      cout << endl;
         //                 refn.print();
                        //   nv->print();
      
                             par->CacheWrite(cache);
                             #pragma omp critical
                          {
         
                          naccepted++;}
                          /*
                           cout << "Predicted Energy Change " << sumE << endl;
                         
                          par->TotalEnergy(4);
                          store-=par->SumStored("Energy");
                         
                         // par->SumStored("Area");
                          Property &m=par->props.GetProperty("MembG",0);
                          cout << "A (pred) " << (double) *m.content << endl;
                     //     m.update();
                          cout << "A " << (double) *m.content << endl;
                          Property &l=par->props.GetProperty("MembG",1);
                          cout << "V (pred) " << (double) *l.content << endl;
                       //   l.update();
                          cout << "V " << (double) *l.content << endl;
                                        conste-=par->TotalEnergy(2);
                            cout << sumE << " CONST " << -conste << " STORE " << -store << endl;
                 //         exit(0);*/
                 }  // Can Be Vectorized!}

            }
        }
        return naccepted / (double) modprop.size();
    }
 
 
private:
    bool init;
    vector<FaceControlProperties> vfcp;
   
};
