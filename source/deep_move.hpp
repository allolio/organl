// These Properties are Dynamically Inserted just to Organize the Simul Face MC Move;
class VertexEnvProperty : public Property
{
public:
    Deplist shallow;
    indexlist depnormals;
    bool update(bool deep)
    {
          for(int i=0; i!=shallow.size(); i++) shallow[i]->update(false); // this will merely take care of blocks etc of original properties.
        if(deep)
            for(int i=0; i!=deps.size(); i++) deps[i]->update();
        return true;
    }
};

class VertexEnvProperties : public PropertyList
{
public:
        virtual bool genpropertiesfromStructure(void *parent);
     VertexEnvProperties()
    {   name="VertexEnv";
        szProperty=6; // Just one Vertex, one norm0l.
        
    }
    FaceVertexMesh *mesh;
  //  ~VertexEnvProperties()
  //  {free(simul);}
  
};


bool VertexEnvProperties::genpropertiesfromStructure(void* parent)
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
  //  simul=(double **) malloc(mesh->faces.size()*6*sizeof(double*));
    theProps.resize(mesh->vertices.size());
    cout << theProps.size() << " " << endl;
    for(int i=0;i!=mesh->vertices.size();i++) // For each Vertex
    {
//        simul[i].resize(6);
        theProps[i]=new VertexEnvProperty;
        theProps[i]->block=NULL;
        theProps[i]->ptype=this;
        theProps[i]->parent=mesh; // Need To decide on that
        theProps[i]->deps.clear();
        theProps[i]->index=i; // NOT SAFE DANGER
        // Assemble Dependents and Co-Updated Normals;
        set<int> dep_normals;
        for(int j=0;j!=mesh->vl[i].size();j++) // Faces of A vertex:
        {
            int faceno=mesh->vl[i][j];
            for(int k=0;k!=mesh->faces[faceno].vert.size();k++)           // Collect Vertices / Normals of all the faces
            {
              int vindex=mesh->faces[faceno].vert[k];
              if(vindex!=i) dep_normals.insert(vindex); // Insert as dependent normals/vertices all vertices on the faces which are not I, make a unique set.
            }
        }
        ((VertexEnvProperty *) theProps[i])->depnormals=indexlist(dep_normals.begin(),dep_normals.end()); // Create A vector of dependent vertices on the thing
        // Now Downward-Merge Dependencies of the Vertices and Normals involved
        //dep_normals.clear();
        for(int k=0;k!=((VertexEnvProperty*) theProps[i])->depnormals.size();k++)
        {
            int vindex=((VertexEnvProperty*) theProps[i])->depnormals[k];
            // Merge Dependencies, e.g. Nagata Triangles into the Props
            MergeDeplists(&(theProps[i]->deps),&(plv->theProps[vindex])->deps);
            MergeDeplists(&(theProps[i]->deps),&(pln->theProps[vindex])->deps);
            ((VertexEnvProperty*) theProps[i])->shallow.push_back(plv->theProps[vindex]); // Take care of blocks and data 
          //  ((FaceControlProperty*) theProps[i])->shallow.push_back(pln->theProps[vindex]); // No Normal Pushing, its DANGER ous....
        }
      //  cout << i << " " << theProps[i]->deps.size()<< " " << endl;
        
        theProps[i]->content=(double *) &mesh->vertices[i].pos;
   //     for(int k=0;k!=3;k++)
  //      {triple *tp=(triple*) (((double **) theProps[i]->content)[k]);// << endl;
        //triple *tp=(triple*)  (theProps[i]->content+9);
 //       tp->print();}
    }
    return true;
};



class DeepVertexMove : public MCMove
{
public:
    DeepVertexMove()
    {
    prop="Vertex";
    tryFrac=0.5;
    maxstepsz=0.1;
    init=false;
//    vfcp.clear();
    }
    bool QuenchNormals()
    {
      VertexEnvProperty *aleprop=(VertexEnvProperty *) modprop[0];
      FaceVertexMesh *afv=(FaceVertexMesh *) aleprop->parent;
      afv->ResetNormals();
      return true;  
    }
    
        bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }

        int CreateParallelPlan(Simulation *sim) 
    {
        if(!init or cycle==-1)
        {
            cout << "Injecting Vertex Environments" << endl;
            venvp.resize(sim->Entities.size());
            // Now Inject Control Properties :)
            for(int i=0;i!=sim->Entities.size();i++)
            {
             if(CheckEntitySupportsMove(sim->Entities[i]))
             { venvp[i].genpropertiesfromStructure(&(sim->Entities[i]->props));
               sim->Entities[i]->props.AddPropertyType(venvp[i]); 
              }
            }
                        cout << "Done Injecting Vertex Environments" << endl;
          prop="VertexEnv"; init=true;
          int rv=MCMove::CreateParallelPlan(sim); 
	  cout << "Quenching Normals" << endl;
	  QuenchNormals();
          return rv;
        }
       return MCMove::CreateParallelPlan(sim); // This allows us to use the default MC Planning Scheme.
    }
    
     double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        
    //    if(QuenchNormals())
     //   {
     //       VertexEnvProperty *aleprop=(VertexEnvProperty *) modprop[0];
     //       FaceVertexMesh *afv=(FaceVertexMesh *) aleprop->parent;
     //       afv->ResetNormals();
     //   }   
        // foreach_move
        // Backup Copy
 // cout << "EXECUTING" << endl;
  //      cout << "MPS" << modprop.size() << endl; 
#pragma omp parallel for
        
        for (int i=0;i<modprop.size();i++)
        {
                         bool accepted=false;
                         int counter=0;
                                    bool first=true;
            while(!accepted and counter < 20)
            {
            VertexEnvProperty *leprop=(VertexEnvProperty *) modprop[i];
            FaceVertexMesh *fv=(FaceVertexMesh *) leprop->parent;
            vector<triple> oldall(leprop->depnormals.size());
            vector<triple> newall(leprop->depnormals.size());
            triple oldv(* (triple*)leprop->content);
            triple newv;
            triple oldnorm=fv->vertices[leprop->index].norm;
            // Generate Vertex Move
         //   GenNewP(oldv.p,newv.p,3);           
            // now backup copy normals
            for(int k=0;k!=oldall.size();k++)
            {          
                oldall[k]=fv->vertices[leprop->depnormals[k]].norm;
            }
//          TO help out with the lambdas   
            vector<double> collect(2,0);
//cout << "GOT HERE EP" << endl;
            double sumE=0;
           // Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Subtract-pre-penalty;
            for(int j=0;j!=modnodes[i].size();j++)
            {
               if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles       
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
  //          cout << "GOT HERE EP" << endl;
            // Normal Indices:
            int retcon=0;
           bool allowed=true;
           // We will now set the new vertex point
           retry:
           
//           newv=oldv;
            GenNewP((double *) &oldv, (double *) &newv,3);          
   //          newv=oldv;
            memcpy(leprop->content, newv.p, sizeof(double)*3);
            
            //triple *nv=&((Vertex *) normal->parent)->norm;
           // Build All Normal Vectors;
         //   cout << i << " "<< leprop->index << endl;
         //   oldnorm.print();
           // cout << " Oldnorm ?" << endl;
         //               cout << " Newnorm ?" << endl;
            allowed=fv->BuildNormal(leprop->index,true); //true
            //if(!allowed) cout << "NormalFail" << endl;
           // cout << " ----------------------------" << endl;
           // oldv.print();
           // newv.print();
           // cout << " ----------------------------" << endl;

          // fv->vertices[nindex->at(0)].pos.print();
 //                       fv->vertices[nindex->at(1)].pos.print();
  //          fv->vertices[nindex->at(2)].pos.print();
   //          cout<< nindex->at(0) << " " << nindex->at(1) << " " << nindex->at(2) << endl;            
            if(!allowed)
            {
             
             if(retcon<10) {retcon++; goto retry;}
             else{             
               memcpy(leprop->content,(double *) &oldv,sizeof(double)*3);
               fv->vertices[leprop->index].norm=oldnorm; first=false;counter++;
                }       
            }
            else 
            {
            // Only now we will reconstruct all the normals.
            for(int m=0;m!=leprop->depnormals.size();m++)
            {
               allowed=allowed and fv->BuildNormal(leprop->depnormals[m],true);        
             //              if(!allowed) cout << "SideNormalFail" << endl;

            }
           // doublemod[i].first->update();
            modprop[i]->update();
            
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].begin(),modnodes[i].end());
            //bool allowed=true;
            //if(*nv * (*nv) < 0.999) allowed=false;
            if(allowed) for(int j=0;j!=modnodes[i].size();j++)
            {
              if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles            
              allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]); //  and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], , &refn);
                     //     if(!allowed) cout << "AngleFail" << endl;
            //  if(!LinearEdge((NTriangle *) sim->n.col[modnodes[i].at(j)])) {allowed=false; cout << "EdgeFail" << endl;}
              if(!allowed) break;

               // sumE+=EdgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
              int index=modnodes[i].at(j);
              if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
              {index=sim->n.col[modnodes[i].at(j)]->node_id;
              sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
              }
              else{
                 index=modnodes[i].at(j)*-1;
              sumE+=par->DExtraE(index);
              cache[j]=-index;
             
              }
              //  cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
                for(int j=0;j!=modnodes[i].size();j++)
                {
                 if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                }
            accepted=AcceptMove(sumE);
                
            }
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                                   if(modnodes[i].at(k) >-1) // Edge penalty in the sense of Nagata triangles
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
//		fv->SaveObj("step1.obj"); exit(0);
            }
            if(!accepted)
            {
          
               for(int k=0;k!=leprop->depnormals.size();k++)
                 fv->vertices[leprop->depnormals[k]].norm=oldall[k];
            //       oldnorm.print();
            
             //               cout << " A---------------------------" << endl;
             //  (fv->vertices[leprop->index].pos).print();
               memcpy(leprop->content,(double *) &oldv,sizeof(double)*3);
               fv->vertices[leprop->index].norm=oldnorm;
             //  oldv.print();
               //                            cout << " E---------------------------" << endl;

             // doublemod[i].first->update();
               modprop[i]->update();
               //((triple *) leprop->content)->print(); 
                first=false;
                counter++; 
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
         if(first)
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
      } }
        return naccepted / (double) modprop.size();
    }
 
 
private:
    bool init;
    vector<VertexEnvProperties> venvp;
   
};


class ANormalMove : public MCMove
{
public:
    ANormalMove()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Normal";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.2;
    }
    
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move//
         double var=maxstepsz*maxstepsz/9.;
     //   uniform_real_distribution<> dis(-maxstepsz, maxstepsz);

        //sim->Entities[0]->DeltaE(0); return -1;
        int szProp=modprop[0]->ptype->szProperty; // Only One Property Type in This example;
    
#pragma omp parallel for
        for (int i=0;i<modprop.size();i++)
        {
                     bool accepted=false;
                                    bool first=true;
                                    int counter =0;
            while(!accepted and counter < 20)
            {
            double newp[3];newp[0]=0;newp[1]=0;newp[2]=0;
            double sumE=0;
          //  Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
                
            Property* normal=(Property*) modprop[i];
            triple *nv=&((Vertex *) normal->parent)->norm;
            //cout << sim->n.col[modnodes[i]].size() << " TOCHECK" << endl;
    
                 for(int j=0;j!=modnodes[i].size();j++)
            {
                 if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
            
            triple oldn=*nv;
            GenNewP(newp,newp,3);
            double gamma=0;
            GenNewP(&gamma,&gamma,1);
            triple newn(newp);
            newn=*nv+newn*gamma;
            *nv=newn/newn.abs();
            bool allowed=true;
   
        //    if(*nv *newn < 0) allowed=false;
            
            vector<double> collect(2,0);
            // update_coefficients.
            ((Vertex *) normal->parent)->norm=*nv;
             modprop[i]->update();
             
             if(*nv * (*nv) < 0.999) allowed=false;
	    // cout << *nv * (*nv) << endl;
        //    cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            for(int j=0;j!=modnodes[i].size();j++)
            {
                if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                {     allowed=allowed and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv,&oldn);
              if(!allowed) break;
               // if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                    allowed=allowed and allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;
                }
                else
                {
               int index=modnodes[i].at(j)*-1;
               sumE+=par->DExtraE(index);
               cache[j]=-index;
               }
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
         
            }
            if(allowed)
            {
            sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
              for(int j=0;j!=modnodes[i].size();j++)
                {
                 if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                }
            }
              accepted = AcceptMove(sumE) and allowed; //:wqg NOrmals
              //AcceptMove(sumE) and allowed;
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(modnodes[i].at(k) >-1) // Edge penalty in the sense of Nagata triangles
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            if(!accepted)
            {
                          ((Vertex *) normal->parent)->norm=oldn;
              modprop[i]->update();
              first=false;counter++;
             // TestReservation
           //  cout << "CUT" << endl;
           // for (int x=0;x!=cache.size();x++)  par->DeltaE(cache[x]);
            // cout << "END CUT" << endl;
            }
            else {    //    cout << "CACHE_WRITE" << endl;
                    //      double conste=par->TotalEnergy(2); 
                      //    cout << "Vertex " << modprop[i]->index << " ";
                      //     par->SumStored("Area");
                       //    par->SumStored("Volume");
                       // (*nv).print();
                       // cout << (*nv).abs() << endl;
                        //  for (int i=0;i!=cache.size();i++) cout << cache[i] << " " ;
                        //      cout << endl;
                        //     double store=par->TotalEnergy(4);
                         // store-=par->SumStored("Energy");
                          
                             par->CacheWrite(cache);
                             #pragma omp critical
                          {
         if(first)
                          naccepted++;}
                          
                         //  cout << "Predicted Energy Change " << sumE << endl;
                         
                        
                         
                         // par->SumStored("Area");
                        //  Property &m=par->props.GetProperty("MembG",0);
                      //    cout << "A (pred) " << (double) *m.content << endl;
                     //     m.update();
                      //    cout << "A " << (double) *m.content << endl;
                       //   Property &l=par->props.GetProperty("MembG",1);
                       //   cout << "V (pred) " << (double) *l.content << endl;
                       //   l.update();
                        //  cout << "V " << (double) *l.content << endl;
                         //               conste-=par->TotalEnergy(2);
                          //  cout << sumE << " CONST " << -conste << " STORE " << -store << "RES" << conste+store << endl;
                 //         exit(0);*/
                 }  // Can Be Vectorized!}

         
        }  }
        return naccepted / (double) modprop.size();
    }
 
protected:
  
private:
    };
