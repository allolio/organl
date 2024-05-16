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

class NormalTravelMove : public MCMove
{
public:
    NormalTravelMove()
    {
    numMoves=0;
    tryFrac=0.5;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.01;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }
    
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int CreateParallelPlan(Simulation *sim) 
    {
        modprop.clear(); modnodes.clear(); busynodes.clear(); doublemod.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            PropertyList *npl;
            string prop2="Normal";
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            sim->Entities[i]->props.GetPropertyList(prop2,&npl);
            int numProps=pl->theProps.size();
            uniform_real_distribution<> dis_u(0, 1);
#pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    Property *p2;
                    indexlist reservation;
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(select,&p);
                    npl->pofs(select,&p2);
                //    cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                        //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteDoubleReservation(p,p2,reservation); //Not
                    }
               }
            }
            // FIXME Separate plans for separate entities!
        }
       return doublemod.size(); 
    }
 
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
         
        //sim->Entities[0]->DeltaE(0); return -1;
       const int szProp=3 ; //doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp];             double oldn[2]; 

             double newno[3];newno[0]=0;newno[1]=0;newno[2]=0;
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            memcpy(oldn, doublemod[i].second->content,2*sizeof(double));

            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
           
            GenNewP(newno,newno,szProp);
            
            triple *nv=&((Vertex *) normal->parent)->norm;
            triple oldnv=*nv;
            triple newn(newno);
            newn=(newn*(1./maxstepsz/8.))+*nv;
            newn=newn/newn.abs();
            
            duple newnorm=normang(newn);
            
            double step=0, origin=0, dtheta=0, dphi=0;
            
            GenNewP(&origin,&step,1);
            
            vector<double> collect(2,0);
              triple newp(oldp); 
              newp=newp+(newn*step);
            // update_coefficients.
            memcpy(doublemod[i].first->content,newp.p,szProp*sizeof(double));
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
           // doublemod[i].first->update();
            doublemod[i].second->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            bool allowed=true;
            for(int j=0;j!=modnodes[i].size();j++)
            {
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

              allowed=allowed and (allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]) and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv,&oldnv));
              if(!allowed) break;
            }
            sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
              bool accepted = AcceptMove(sumE) and allowed;
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            if(!accepted)
            {
              memcpy(doublemod[i].first->content,oldp,szProp*sizeof(double));
              memcpy(doublemod[i].second->content,oldn,2*sizeof(double));
             // doublemod[i].first->update();

              doublemod[i].second->update();
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
        return naccepted / (double) doublemod.size();
    }
 
protected:
  
    double tryFrac;
    string prop;
    vector < pair < Property *, Property * > > doublemod; 
   
private:
      bool WriteDoubleReservation(Property *p, Property *p2, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        doublemod.push_back(make_pair(p,p2));
        modnodes.push_back(reservation);
        return true;
    }
  

          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
    
};



class FaceWalker : public MCMove
{
public:
    FaceWalker()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.3;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int CreateParallelPlan(Simulation *sim) 
    {
        modprop.clear(); modnodes.clear(); busynodes.clear(); doublemod.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            PropertyList *npl;
            string prop2="Normal";
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            sim->Entities[i]->props.GetPropertyList(prop2,&npl);
            int numProps=pl->theProps.size();
            uniform_real_distribution<> dis_u(0, 1);
#pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    Property *p2;
                    indexlist reservation;
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(select,&p);
                    npl->pofs(select,&p2);
                //    cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                        //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteDoubleReservation(p,p2,reservation); //Not
                    }
               }
            }
            // FIXME Separate plans for separate entities!
        }
       return doublemod.size(); 
    }
 
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
              double var=maxstepsz*maxstepsz/9.;
              normal_distribution<> dis(0.0, var);
       
        //sim->Entities[0]->DeltaE(0); return -1;
       const int szProp=3; //doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp];             double oldn[2]; 

             double newpl[2];newpl[0]=0;newpl[1]=0;//newno[2]=0;
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            memcpy(oldn, doublemod[i].second->content,2*sizeof(double));

            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
            
            

            //GenNewP(newno,newno,szProp);
            
            triple *nv=&((Vertex *) normal->parent)->norm;
          
            // Select Random Face
            uniform_real_distribution<> dis_u(0, modnodes[i].size());
            int face=dis_u(gen);
           // roll:
            double xi=1,eta=1;
    /*        while(xi*xi > 1 ) xi=dis(gen);
            while(eta*eta > xi*xi ) eta=dis(gen);
            
            
   */
            NTriangle *nt=(NTriangle *) sim->n.col[modnodes[i].at(face)];
            triple npos=nt->eval(0.6666666,0.33333333);
            triple newn=nt->normal(0.666666,0.3333333);
    //        double eta=dis(gen);
            triple oldnv = *nv;
            if(newn.abs()<0.99) newn=*nv;
            if(*nv *newn <0) newn=newn*-1;
            newn=newn+(*nv*0.1);
            newn=newn/newn.abs();
            
            duple newnorm=normang(newn);
            
            double step=0, origin=0, dtheta=0, dphi=0;
            
            vector<double> collect(2,0);
            // update_coefficients.
            memcpy(doublemod[i].first->content,npos.p,szProp*sizeof(double));
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
           // doublemod[i].first->update();
            doublemod[i].second->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            bool allowed=true;
            for(int j=0;j!=modnodes[i].size();j++)
            {
              allowed=allowed and (allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]) and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv, &oldnv));
              if(!allowed) break;

                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                                //((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
            if(allowed) sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
              bool accepted = AcceptMove(sumE) and allowed;
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            if(!accepted)
            {
              memcpy(doublemod[i].first->content,oldp,szProp*sizeof(double));
              memcpy(doublemod[i].second->content,oldn,2*sizeof(double));
             // doublemod[i].first->update();

              doublemod[i].second->update();
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
        return naccepted / (double) doublemod.size();
    }
 
protected:
  
    double tryFrac;
    string prop;
    vector < pair < Property *, Property * > > doublemod; 
   
private:
      bool WriteDoubleReservation(Property *p, Property *p2, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        doublemod.push_back(make_pair(p,p2));
        modnodes.push_back(reservation);
        return true;
    }
  
     
          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
    
};




class VertexBuildMove : public MCMove
{
public:
    VertexBuildMove()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.01;
     avlength=706;
     kappa=0.1*avlength;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int CreateParallelPlan(Simulation *sim) 
    {
        modprop.clear(); modnodes.clear(); busynodes.clear(); doublemod.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            PropertyList *npl;
            string prop2="Normal";
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            sim->Entities[i]->props.GetPropertyList(prop2,&npl);
            int numProps=pl->theProps.size();
            fv = (FaceVertexMesh*) pl-> parent;
            uniform_real_distribution<> dis_u(0, 1);
#pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    Property *p2;
                    indexlist reservation;
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(select,&p);
                    npl->pofs(select,&p2);
                //    cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                        //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteDoubleReservation(p,p2,reservation); //Not
                    }
               }
            }
            // FIXME Separate plans for separate entities!
        }
       return doublemod.size(); 
    }
 
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
         
        //sim->Entities[0]->DeltaE(0); return -1;
       const int szProp=3; //doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp];             double oldn[2]; 
            double newp[szProp];
             double newno[3];newno[0]=0;newno[1]=0;newno[2]=0;
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            memcpy(oldn, doublemod[i].second->content,2*sizeof(double));

            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
                             for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }

            //triple center; center=center*0;
            /*
            for(int j=0;j!=modnodes[i].size();j++)
            {
            center=center+((NTriangle *) sim->n.col[modnodes[i].at(j)])->eval(0.666666,0.333333);
            }
            center=center/modnodes[i].size();
           */
            GenNewP(oldp,newp,szProp);
            
            triple *nv=&((Vertex *) normal->parent)->norm;
        
            triple refn=*nv;
//            cout << normal->index << endl;
            bool allowed=fv->BuildNormal(normal->index,true,&refn);          
            duple newnorm=normang(*nv);
                vector<double> collect(2,0);
            
            
            // update_coefficients.
            memcpy(doublemod[i].first->content,newp,szProp*sizeof(double));
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
            doublemod[i].first->update();
           // doublemod[i].second->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            //bool allowed=true;
            if(*nv * (*nv) < 0.999) allowed=false;
            if(allowed) for(int j=0;j!=modnodes[i].size();j++)
            {
               allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]) and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv, &refn);
              if(!allowed) break;

               // sumE+=EdgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
              bool accepted=false;
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
               accepted = AcceptMove(sumE);
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            if(accepted)
            {
                for(int lo=0;lo!=fv->vertices[normal->index].neighbors.size();lo++)
                    if(fv->QualityAroundVertex(fv->vertices[normal->index].neighbors[lo]) < 1) 
                    {accepted=false; break;}
            }
            }
            if(!accepted)
            {
              memcpy(doublemod[i].first->content,oldp,szProp*sizeof(double));
              ((Vertex *) normal->parent)->norm=refn;
              memcpy(doublemod[i].second->content,oldn,2*sizeof(double));
              doublemod[i].first->update();

             // doublemod[i].second->update();
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
        return naccepted / (double) doublemod.size();
    }
 
    double avlength;
    double kappa;
    
protected:
  
    double tryFrac;
    string prop;
    vector < pair < Property *, Property * > > doublemod; 
   
private:
    FaceVertexMesh *fv;
 /*   double EdgePenalty(NTriangle *nt)
    {
        double penalty=0;
        double d=((*nt->v.at(0))-(*nt->v.at(1))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
         d=((*nt->v.at(0))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
                d=((*nt->v.at(1))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
        return penalty;
    }
    */ //Was in Trial 21
      bool WriteDoubleReservation(Property *p, Property *p2, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        doublemod.push_back(make_pair(p,p2));
        modnodes.push_back(reservation);
        return true;
    }
  

          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
    
};



class VertexBuildMidMove : public MCMove
{
public:
    VertexBuildMidMove()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.1;
     avlength=706;
     kappa=0.1*avlength;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int CreateParallelPlan(Simulation *sim) 
    {
        modprop.clear(); modnodes.clear(); busynodes.clear(); doublemod.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            PropertyList *npl;
            string prop2="Normal";
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            sim->Entities[i]->props.GetPropertyList(prop2,&npl);
            int numProps=pl->theProps.size();
            fv = (FaceVertexMesh*) pl-> parent;
            uniform_real_distribution<> dis_u(0, 1);
#pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    Property *p2;
                    indexlist reservation;
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(select,&p);
                    npl->pofs(select,&p2);
                //    cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                        //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteDoubleReservation(p,p2,reservation); //Not
                    }
               }
            }
            // FIXME Separate plans for separate entities!
        }
       return doublemod.size(); 
    }
 
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
         
        //sim->Entities[0]->DeltaE(0); return -1;
       const int szProp=3;//doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp];             double oldn[2]; 
            double newp[szProp];
             double newno[3];newno[0]=0;newno[1]=0;newno[2]=0;
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            memcpy(oldn, doublemod[i].second->content,2*sizeof(double));

            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
            triple center; center=center*0;
            
                        for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }

            
            for(int j=0;j!=modnodes[i].size();j++)
            {
            center=center+((NTriangle *) sim->n.col[modnodes[i].at(j)])->eval(0.666666,0.333333);
            }
            center=center/modnodes[i].size();
           
            GenNewP(center.p,newp,szProp);
            
            triple *nv=&((Vertex *) normal->parent)->norm;
        
            triple refn=*nv;
//            cout << normal->index << endl;
            bool allowed=fv->BuildNormal(normal->index,true,&refn);          
            duple newnorm=normang(*nv);
            
            
            vector<double> collect(2,0);
            
            
            // update_coefficients.
            memcpy(doublemod[i].first->content,newp,szProp*sizeof(double));
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
           // doublemod[i].first->update();
            doublemod[i].second->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            if(allowed)
            for(int j=0;j!=modnodes[i].size();j++)
            {
               allowed=(allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]) and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv, &refn));
              if(!allowed) break;

               // sumE+=EdgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
            bool accepted=false;
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
               accepted = AcceptMove(sumE);
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            }
            if(!accepted)
            {
              memcpy(doublemod[i].first->content,oldp,szProp*sizeof(double));
              memcpy(doublemod[i].second->content,oldn,2*sizeof(double));
             // doublemod[i].first->update();

              doublemod[i].second->update();
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
        return naccepted / (double) doublemod.size();
    }
 
    double avlength;
    double kappa;
    
protected:
  
    double tryFrac;
    string prop;
    vector < pair < Property *, Property * > > doublemod; 
   
private:
    FaceVertexMesh *fv;
    double EdgePenalty(NTriangle *nt)
    {
        double penalty=0;
        double d=((*nt->v.at(0))-(*nt->v.at(1))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
         d=((*nt->v.at(0))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
                d=((*nt->v.at(1))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
        return penalty;
    }
    
      bool WriteDoubleReservation(Property *p, Property *p2, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        doublemod.push_back(make_pair(p,p2));
        modnodes.push_back(reservation);
        return true;
    }
  

          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
    
};





class VertexFaceMidMove : public MCMove
{
public:
    VertexFaceMidMove()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.1;
     avlength=706;
     kappa=0.1*avlength;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int CreateParallelPlan(Simulation *sim) 
    {
        modprop.clear(); modnodes.clear(); busynodes.clear(); doublemod.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            PropertyList *npl;
            string prop2="Normal";
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            sim->Entities[i]->props.GetPropertyList(prop2,&npl);
            int numProps=pl->theProps.size();
            fv = (FaceVertexMesh*) pl-> parent;
            uniform_real_distribution<> dis_u(0, 1);
#pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    Property *p2;
                    indexlist reservation;
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(select,&p);
                    npl->pofs(select,&p2);
                //    cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                        //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteDoubleReservation(p,p2,reservation); //Not
                    }
               }
            }
            // FIXME Separate plans for separate entities!
        }
       return doublemod.size(); 
    }
 
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
         
        //sim->Entities[0]->DeltaE(0); return -1;
      const  int szProp=3;//doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp];             double oldn[2]; 
            double newp[szProp];
             double newno[3];newno[0]=0;newno[1]=0;newno[2]=0;
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            memcpy(oldn, doublemod[i].second->content,2*sizeof(double));

            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
            triple center; center=center*0;
            
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
            triple *nv=&((Vertex *) normal->parent)->norm;
        
            triple refn=*nv;
            
            *nv=*nv*0;
            
            for(int j=0;j!=modnodes[i].size();j++)
            {
            NTriangle *nt= (NTriangle *) sim->n.col[modnodes[i].at(j)];
            int l=1,k=2;int mye=0;
            for(int i=0;i!=3;i++) { if(nt->n[i]==nv) {mye=i; l=l%3;k=k%3; break; } l++; k++;}
            center=center+*(nt->v[l])+*(nt->v[k]);
            *nv=*nv+nt->normal(0.666666,0.333333);
            }
            *nv=*nv/nv->abs(); if(*nv *refn <0) *nv=*nv *-1;
            center=center/modnodes[i].size()/2.;
            
            GenNewP(center.p,newp,szProp);
            
       
//            cout << normal->index << endl;
            //fv->BuildNormal(normal->index,true,&refn);          
            duple newnorm=normang(*nv);
            
            
            vector<double> collect(2,0);
            
            
            // update_coefficients.
            memcpy(doublemod[i].first->content,newp,szProp*sizeof(double));
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
           // doublemod[i].first->update();
            doublemod[i].second->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            bool allowed=true;
            for(int j=0;j!=modnodes[i].size();j++)
            {
               allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]) and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv, &refn);
              if(!allowed) break;

               // sumE+=EdgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
            bool accepted=false;
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
               accepted = AcceptMove(sumE);
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            }
            if(!accepted)
            {
              memcpy(doublemod[i].first->content,oldp,szProp*sizeof(double));
              memcpy(doublemod[i].second->content,oldn,2*sizeof(double));
             // doublemod[i].first->update();

              doublemod[i].second->update();
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
        return naccepted / (double) doublemod.size();
    }
 
    double avlength;
    double kappa;
    
protected:
  
    double tryFrac;
    string prop;
    vector < pair < Property *, Property * > > doublemod; 
   
private:
    FaceVertexMesh *fv;
    double EdgePenalty(NTriangle *nt)
    {
        double penalty=0;
        double d=((*nt->v.at(0))-(*nt->v.at(1))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
         d=((*nt->v.at(0))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
                d=((*nt->v.at(1))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
        return penalty;
    }
    
      bool WriteDoubleReservation(Property *p, Property *p2, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        doublemod.push_back(make_pair(p,p2));
        modnodes.push_back(reservation);
        return true;
    }
  

          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
    
};




class VertexBuildMoveWC : public MCMove
{
public:
    VertexBuildMoveWC()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.01;
     avlength=706;
     kappa=0.1*avlength;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int CreateParallelPlan(Simulation *sim) 
    {
        modprop.clear(); modnodes.clear(); busynodes.clear(); doublemod.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            PropertyList *npl;
            string prop2="Normal";
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            sim->Entities[i]->props.GetPropertyList(prop2,&npl);
            int numProps=pl->theProps.size();
            fv = (FaceVertexMesh*) pl-> parent;
            uniform_real_distribution<> dis_u(0, 1);
#pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    Property *p2;
                    indexlist reservation;
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(select,&p);
                    npl->pofs(select,&p2);
                //    cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                        //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteDoubleReservation(p,p2,reservation); //Not
                    }
               }
            }
            // FIXME Separate plans for separate entities!
        }
       return doublemod.size(); 
    }
 
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
         
        //sim->Entities[0]->DeltaE(0); return -1;
      const  int szProp=3;//doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp];             double oldn[2]; 
            double newp[szProp];
             double newno[3];newno[0]=0;newno[1]=0;newno[2]=0;
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            memcpy(oldn, doublemod[i].second->content,2*sizeof(double));

            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
                             for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }

            //triple center; center=center*0;
            /*
            for(int j=0;j!=modnodes[i].size();j++)
            {
            center=center+((NTriangle *) sim->n.col[modnodes[i].at(j)])->eval(0.666666,0.333333);
            }
            center=center/modnodes[i].size();
           */
            GenNewP(oldp,newp,szProp);
            
            triple *nv=&((Vertex *) normal->parent)->norm;
        
            triple refn=*nv;
//            cout << normal->index << endl;
            bool allowed=fv->BuildNormal(normal->index,true,&refn);          
            duple newnorm=normang(*nv);
                vector<double> collect(2,0);
            
            
            // update_coefficients.
            memcpy(doublemod[i].first->content,newp,szProp*sizeof(double));
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
           // doublemod[i].first->update();
            doublemod[i].second->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            triple cnv;cnv=cnv*0;
            for(int j=0;j!=modnodes[i].size();j++)
            {
              triple nnv;
              // bool CornerNormal(NTriangle *nt, triple *vertex, triple *normal)
              CornerNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], (triple*) doublemod[i].first->content , &nnv);
              cnv=cnv+nnv;
            }
            *nv=cnv/cnv.abs();
            newnorm=normang(*nv);
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
            doublemod[i].second->update();
            
            indexlist cache(modnodes[i].size());
            //bool allowed=true;
            if(*nv * (*nv) != *nv * (*nv)) allowed=false;
            if(*nv * (*nv) < 0.999) allowed=false;
           // nv->print();
          //  cout << "ALLOWED " << allowed << endl;
            if(allowed) for(int j=0;j!=modnodes[i].size();j++)
            {
               allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]) and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv, &refn);
              if(!allowed) break;

               // sumE+=EdgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
              bool accepted=false;
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
               accepted = AcceptMove(sumE);
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            }
            if(!accepted)
            {
              memcpy(doublemod[i].first->content,oldp,szProp*sizeof(double));
              memcpy(doublemod[i].second->content,oldn,2*sizeof(double));
             // doublemod[i].first->update();

              doublemod[i].second->update();
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
        return naccepted / (double) doublemod.size();
    }
 
    double avlength;
    double kappa;
    
protected:
  
    double tryFrac;
    string prop;
    vector < pair < Property *, Property * > > doublemod; 
   
private:
    FaceVertexMesh *fv;
 /*   double EdgePenalty(NTriangle *nt)
    {
        double penalty=0;
        double d=((*nt->v.at(0))-(*nt->v.at(1))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
         d=((*nt->v.at(0))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
                d=((*nt->v.at(1))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
        return penalty;
    }
    */ //Was in Trial 21
      bool WriteDoubleReservation(Property *p, Property *p2, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        doublemod.push_back(make_pair(p,p2));
        modnodes.push_back(reservation);
        return true;
    }
  

          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
    
};


class VertexBuildMoveON : public MCMove
{
public:
    VertexBuildMoveON()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.01;
     avlength=706;
     kappa=0.1*avlength;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    bool CheckEntitySupportsMove(Entity *ent) 
    {
       return ent->props.HasProp(prop);
    }
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int CreateParallelPlan(Simulation *sim) 
    {
        modprop.clear(); modnodes.clear(); busynodes.clear(); doublemod.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            PropertyList *npl;
            string prop2="Normal";
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            sim->Entities[i]->props.GetPropertyList(prop2,&npl);
            int numProps=pl->theProps.size();
            fv = (FaceVertexMesh*) pl-> parent;
            uniform_real_distribution<> dis_u(0, 1);
#pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    Property *p2;
                    indexlist reservation;
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(select,&p);
                    npl->pofs(select,&p2);
                //    cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                        //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteDoubleReservation(p,p2,reservation); //Not
                    }
               }
            }
            // FIXME Separate plans for separate entities!
        }
       return doublemod.size(); 
    }
 
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
         
        //sim->Entities[0]->DeltaE(0); return -1;
      const  int szProp=3;//doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp];             double oldn[2]; 
            double newp[szProp];
             double newno[3];newno[0]=0;newno[1]=0;newno[2]=0;
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            memcpy(oldn, doublemod[i].second->content,2*sizeof(double));

            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
                             for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }

            //triple center; center=center*0;
            /*
            for(int j=0;j!=modnodes[i].size();j++)
            {
            center=center+((NTriangle *) sim->n.col[modnodes[i].at(j)])->eval(0.666666,0.333333);
            }
            center=center/modnodes[i].size();
           */
            GenNewP(oldp,newp,szProp);
            
            triple *nv=&((Vertex *) normal->parent)->norm;
        
            triple refn=*nv;
//            cout << normal->index << endl;
            bool allowed=fv->BuildNormal(normal->index,true,&refn);          
            duple newnorm=normang(*nv);
                vector<double> collect(2,0);
            
            
            // update_coefficients.
            memcpy(doublemod[i].first->content,newp,szProp*sizeof(double));
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
           // doublemod[i].first->update();
            doublemod[i].second->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            triple cnv;cnv=cnv*0;
            for(int j=0;j!=modnodes[i].size();j++)
            {
              triple nnv;
              // bool CornerNormal(NTriangle *nt, triple *vertex, triple *normal)
              OuterNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], (triple*) doublemod[i].first->content , &nnv);
              cnv=cnv+nnv;
            }
            *nv=((*nv)+(cnv/cnv.abs()));
            *nv=(*nv)/nv->abs();
            newnorm=normang(*nv);
            memcpy(doublemod[i].second->content,newnorm.p,2*sizeof(double));
            doublemod[i].second->update();
            
            indexlist cache(modnodes[i].size());
            //bool allowed=true;
            if(*nv * (*nv) != *nv * (*nv)) allowed=false;
            if(*nv * (*nv) < 0.999) allowed=false;
           // nv->print();
          //  cout << "ALLOWED " << allowed << endl;
            if(allowed) for(int j=0;j!=modnodes[i].size();j++)
            {
               allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]) and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv, &refn);
              if(!allowed) break;

               // sumE+=EdgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
                              //  triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;

            }
              bool accepted=false;
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
               accepted = AcceptMove(sumE);
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            }
            if(!accepted)
            {
              memcpy(doublemod[i].first->content,oldp,szProp*sizeof(double));
              memcpy(doublemod[i].second->content,oldn,2*sizeof(double));
             // doublemod[i].first->update();

              doublemod[i].second->update();
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
        return naccepted / (double) doublemod.size();
    }
 
    double avlength;
    double kappa;
    
protected:
  
    double tryFrac;
    string prop;
    vector < pair < Property *, Property * > > doublemod; 
   
private:
    FaceVertexMesh *fv;
 /*   double EdgePenalty(NTriangle *nt)
    {
        double penalty=0;
        double d=((*nt->v.at(0))-(*nt->v.at(1))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
         d=((*nt->v.at(0))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
                d=((*nt->v.at(1))-(*nt->v.at(2))).abs();
        penalty+=0.5*kappa*(d-avlength)*(d-avlength);
        return penalty;
    }
    */ //Was in Trial 21
      bool WriteDoubleReservation(Property *p, Property *p2, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        doublemod.push_back(make_pair(p,p2));
        modnodes.push_back(reservation);
        return true;
    }
  

          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
    
};
