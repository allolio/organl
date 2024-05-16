#define NORM_REP_LIM 0.995

      void MCMove::GenNewP(double* oldp, double* newp, int szprop)
    {
         double var=(maxstepsz*maxstepsz/9.); // Three Standard deviations is Maximum step
          normal_distribution<> dis(0.0, var);
      // uniform_real_distribution<> dis(-maxstepsz, maxstepsz);
        for (int i=0;i!=szprop;i++)
        {
            newp[i]=oldp[i]+dis(gen);
        }
    }
     bool MCMove::AcceptMove(double deltaE)
    {
     if(Beta*deltaE>100) return false; // Avoid overflows
     if(Beta*deltaE<0) return true;
    uniform_real_distribution<> dis_u(0, 1);
    return (exp(-Beta*deltaE) > dis_u(gen));
    }

    bool MCMove::Reset(Simulation *sim)
    {
      cycle=-1;cyclemnd.clear(); cyclemp.clear();
      return true;
    }
       
    bool MCMove::TestReservation(indexlist &reservation)
    {
        if(busynodes.empty()) return true;
        for (int i=0;i!=reservation.size();i++)
           if(busynodes.find(reservation[i])!=busynodes.end()) { return false;}
        return true;
    }
    
    bool MCMove::WriteReservation(Property *p, indexlist &reservation)
    {
        busynodes.insert(reservation.begin(),reservation.end());
        modprop.push_back(p);
        modnodes.push_back(reservation);
        return true;
    }
        
    bool MCMove::BuildReservation(Property *p, indexlist &reservation)
    {
        Deplist  d; d.clear(); 
        FullMergeDeplists(&d,p->dependents(),true); // Automatic Recursion here.
 //       cout << d.size() << endl;
        for(int k=0;k!=d.size();k++) 
        {
//                       cout <<  "INDEX" <<  (d.at(k))->index ;
//                       cout <<  " TYPE " << (d.at(k))->ptype->name << " ";
                      if (d.at(k)->ptype->collidable or d.at(k)->ptype->energetic) // In this case it is in the Neighbor List and its value is a collidable
                       {
                           int node_id=d.at(k)->index; // This is quite a weak identification, we will have to "cheat" here.
                           reservation.push_back(node_id); //cout << " X" ;
                       }
 //                      cout << endl;
        }
                
    return true;
    }
    
     bool MCMove::allowNormal(NTriangle *nt, triple *norm, triple *nOld)
   {
       int mye=0,j=1,k=2;
       for(int i=0;i!=3;i++) { if(nt->n[i]==norm) {mye=i; j=j%3;k=k%3; break; } j++; k++;}
       if(k>=5) return true;
    
       
   // triple y = (*nt->v[k] - *nt->v[mye])
  //      % (*nt->v[j] - *nt->v[mye]);
  //    y=y/y.abs();
//      cout << (y* *nOld) << " " <<  (y* *norm) << " " << y*y  << endl;
//      if (!( (y* *nOld) * (y* *norm) > 0.1 )) return false;
      if (*norm * *nt->n[j] < 0.01 or *norm * *nt->n[k] < 0.01) return false; // 0.01 Before
      
      // Addition: Test reproduction of normals!
            double sc1=nt->normal(0,0)* *nt->n[0];
            double sc2=nt->normal(1,0)* *nt->n[1];
            double sc3=nt->normal(1,1)* *nt->n[2];
          //  cout << "sc1" << sc1 << " " << sc2 << " " << sc3 << endl; 
            if ( sc1 < NORM_REP_LIM or sc2 < NORM_REP_LIM or sc3 < NORM_REP_LIM) return false;
      
      return true;
   }
   
       
 double MCMove::edgePenalty(NTriangle *n)
 {
   //  if(llim == NO_LIMIT and hlim==NO_LIMIT) return 0;
     double penalty=0;
     triple ta  = *(n->v[0])-*(n->v[1]);
     triple tc  = *(n->v[2])-*(n->v[0]);
     triple tb  = *(n->v[1])-*(n->v[2]);
     double kappa=n->GetProperty("kappa");
// DEBUG FOR PURE HELFRICH
     double A0=n->GetProperty("A0");
     double c0=n->GetProperty("c0");
     double Ka=n->GetProperty("Ka");
 //    c0=0; // This needs to be rethought :D
     double A=n->GetProperty("Try_Area");
    //cout << "A" << A << endl;
     double I=2.*((kappa+1)*(1+A0*c0*c0)+A0*A0*Ka);
    // Energy Intensity 
     if( A < 0.4*A0 ) penalty+=(A-A0)*(A-A0)*kpenal*I*10; // 0.5 // Somewhat  new 1e3 to function as effective BAN
     penalty+=I*100*LinearEdge(n);
     // if( A > 1.2*A0 ) penalty+=(A-A0)*(A-A0)*kpenal*I*10; // Somewhat 
   //  double KA=n->GetProperty("Ka");
 //    penalty +=KA*(A0-A)*(A0-A);
//     if( A0 < 0.8*An0 ) penalty+=(A0-An0)*(A0-An0)*kpenal;
// Read-In Angle defaults: Penalty
             double a[3]; double b[3];
             a[0]=n->GetProperty("An0");
             if(a[0]!=0)
            {   
             ta=ta/ta.abs(); tb=tb/tb.abs(); tc=tc/tc.abs();
             b[0]=acos(ta*tb); b[1]=acos(tb*tc); b[2]=acos(ta*tc);
             a[1]=n->GetProperty("An1");
             a[2]=n->GetProperty("An2");
             for(int i=0;i!=3;i++)
             {
                 double ang2=(b[i]-a[i])*(b[i]-a[i]);
                 if(ang2 > 0.1) // 0.1
                 penalty+=ang2*kpenal*I/(200*10); // 250
             }
            }
    
    return penalty;
 }
    
    bool MCMove::allowAngle(NTriangle* nagata)
{
  //  return true;
    //cout << "Pen\n";
    double cosAlpha0 = 0.98078528*0.98078528;//0.9238795325;//      // cos(pi/16)      // cos(7.5 deg)
    double step = 0.1; // Originally 0.1             // step for evaluation within nagata triangle
    triple v;
    triple a;
    triple b;

    v = nagata->eval(0, 0);
    triple ta;
    triple tb;
    double cosAlpha;
    a = nagata->eval(step, 0);
    b = nagata->eval(step, step);
    ta = (a - v);
    tb = (b - v);
    triple tc;
    
    //if(hlim != NO_LIMIT or llim != NO_LIMIT) 
    //    if(!edgelimit(nagata)) return false;
    
    cosAlpha = (ta * tb) / (ta.abs() * tb.abs());
    if ( cosAlpha*cosAlpha > cosAlpha0)
    {
        return false;
       // penalty += exp(scale*(cosAlpha - cosAlpha0));
        // cout << "Penalty: angle too low: " << cosAlpha << "\n";
    }

    v = nagata->eval(1, 0);
    a = nagata->eval(1 - step, 0);
    b = nagata->eval(1, step);
    ta = (a - v);
    tb = (b - v);
    cosAlpha = (ta * tb) / (ta.abs() * tb.abs());
    if ( cosAlpha*cosAlpha > cosAlpha0)
    {
        return false;
//        penalty += exp(scale*(cosAlpha - cosAlpha0));
        // cout << "Penalty: angle too low: " << cosAlpha << "\n";
    }

    v = nagata->eval(1, 1);
    a = nagata->eval(1, 1 - step);
    b = nagata->eval(1 - step, 1 - step);
    ta = (a - v);
    tb = (b - v);
    cosAlpha = (ta * tb) / (ta.abs() * tb.abs());
    if ( cosAlpha*cosAlpha > cosAlpha0)
    {
        return false;
//        penalty += exp(scale*(cosAlpha - cosAlpha0));
        // cout << "Penalty: angle too low: " << cosAlpha << "\n";
    }

    //if (penalty > 0)
        //cout << "\t angle penatly: " << penalty << "\n";
      // Addition: Test reproduction of normals!
            double sc1=nagata->normal(0,0)* *nagata->n[0];
            double sc2=nagata->normal(1,0)* *nagata->n[1];
            double sc3=nagata->normal(1,1)* *nagata->n[2];
    //    cout << sc1 << " " << sc2 << " " << sc3 << " " << endl;
           if ( sc1 < NORM_REP_LIM or sc2 < NORM_REP_LIM or sc3 < NORM_REP_LIM) {/*cout << "Reprofail!" << endl*/;return false;}
  //  if(!LinearEdge(nagata)) return false;

    return true;
}

 bool MCMove::CheckEntitySupportsMove(Entity *ent) 
    {
    //    cout << prop << endl;
       return ent->props.HasProp(prop);
    }
    /*  2. Create Parallel Plan Return number of planned moves
     */
    int MCMove::CreateParallelPlan(Simulation *sim) 
    {
        if(cycle==-1)
        {
            cout << "preparing precomputed steps" << endl;
        modprop.clear(); modnodes.clear(); busynodes.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {par=sim->Entities[i];
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            int numProps=pl->theProps.size();
            uniform_real_distribution<> dis_u(0, 1);
            
            deque<int> leftover(numProps);
            deque<int> removal; removal.clear();
            for(int k=0;k!=numProps;k++) leftover[k]=k;
            cout << "numProps " << numProps << endl;
            while(leftover.size()>0)
            {
                int l=0;
 #pragma omp parallel for
            for(int j=0;j<(int) (numProps*tryFrac);j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    indexlist reservation; reservation.clear();
                    int select=(int) (dis_u(gen)*leftover.size()); // select random index.
                 //  cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(leftover[select],&p);
                  //  cout << "BUILDING at" << select << endl;
                    BuildReservation(p,reservation);
                    #pragma omp critical
                    {
                    if(TestReservation(reservation))
                    {  //NOT VECTORIZABLE
                  //      cout << "RESERVING at" << select << endl;
                        WriteReservation(p,reservation); //Not
                        removal.push_back(leftover[select]);}
                    }
               }
               sort(removal.begin(), removal.end(), less<int>());
               // leftover is sorted always.
                 cout << "inside " << leftover.size() << endl; 
               /*
                */
               int i2=0;
           
               for(int i1=0;i1<leftover.size();i1++)
               {
                 if(leftover[i1]==removal[i2]) 
                 {
                   leftover[i1]=-10; i2++; //cout << i1 << " " << i2 << endl;
                 }
                 if(i2>=removal.size()) break;
               }
               leftover.erase( remove (leftover.begin(), leftover.end(), -10), leftover.end() );
                 cout << "remove" << removal.size() << endl;
    //           cout << "leftover" << leftover.size() << endl;
               removal.clear();
      //         cout << "leftover" ;
//               for(int m=0;m!=leftover.size();m++) cout << leftover[m] << " " ;
       //        cout << endl;   
         //      cout << "end leftover" << endl;
                busynodes.clear();
                cyclemp.push_back(modprop);
           cyclemnd.push_back(modnodes); modnodes.clear(); modprop.clear();

            }
            cout << "Created Plan" << endl;
            // FIXME Separate plans for separate entities!
          }
                     cycle=0;
        }}
        modprop=cyclemp[cycle];
        modnodes=cyclemnd[cycle];
        cycle++;
        if(cycle==cyclemnd.size()) cycle=0;        
 
       return cycle; 
    }
 

class VertexMove : public MCMove
{
public:
    VertexMove()
    {

    numMoves=0;
    tryFrac=0.4;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=0.1;
    par=NULL;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
      if(par==NULL)
      { int i=0;
        while(modnodes[0].at(i)<0) {i++;}
        par=sim->Entities[sim->n.parent_entity[modnodes[0].at(i)]];
      }
        // foreach_move
        // cout << "MODPROP " << modprop.size() << endl;
        //sim->Entities[0]->DeltaE(0); return -1;
        const int szProp=3;//modprop[0]->ptype->szProperty; // Only One Property Type in This example;
#pragma omp parallel for
        for (int i=0;i<modprop.size();i++)
        {
                        bool accepted=false;
                                    bool first=true;
                                    int counter=0;
            while(!accepted and counter <20)
            {
            double oldp[szProp]; 
            double newp[szProp];
            memcpy(oldp, modprop[i]->content,szProp*sizeof(double));
            double sumE=0;
            for(int j=0;j!=modnodes[i].size();j++)
            {
              if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
            //
            GenNewP(oldp,newp,szProp);
            vector<double> collect(2,0);
            // update_coefficients.
            memcpy(modprop[i]->content,newp,szProp*sizeof(double));
/*            cout << "UPDATE" << endl;
            for(int x=0;x!=modprop[i]->deps.size();x++)
            {
              cout << modprop[i]->deps[x]->ptype->name << " " << modprop[i]->deps[x]->index << endl;
            }
  
                         cout << "END UPDATE" << endl;
*/
           modprop[i]->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            bool allowed=true;
            /* Check for singular edge
             * */
            vector<triple> fnorm(modnodes[i].size());
           
            for(int j=0;j!=modnodes[i].size();j++)
            {
              int index=0;
              if(modnodes[i].at(j) >-1)
              {index=sim->n.col[modnodes[i].at(j)]->node_id;
                // Check if index belongs to same 
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]);    
              }
              else
              {
              index=modnodes[i].at(j)*-1;
              sumE+=par->DExtraE(index);
              cache[j]=-index;
              }
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
            
              if(!allowed) break;
            }
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
               accepted = AcceptMove(sumE);
            }
            if(accepted)
            {
                // Only now check Singular face you can move this at will
#ifdef EDGE_TILT_LIMIT
              for(int j=0;j!=modnodes[i].size();j++)
            {
                if(modnodes[i].at(j) >-1) // Edge penalty in the sense of Nagata triangles
                {
                fnorm[j]=(*((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[0]-*(((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[1])) % (*((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[0]-*(((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[2]));
                fnorm[j]=fnorm[j]/fnorm[j].abs();
               // if(fnorm[j]* *((NTriangle *) sim->n.col[modnodes[i].at(j)])->n[1] < 0) fnorm[j]=fnorm[j]*-1;
                for(int k=0;k<j;k++)
                    if((fnorm[k]*fnorm[j])*(fnorm[k]*fnorm[j])<EDGE_TILT_LIMIT) {accepted=false; break;}
                } // FIXME to account for neighbors that are not triangles.
                if(!accepted) break;
            }
#endif
            } // Operated this out as this should be caught by normal verifcation.
            if(accepted)   
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(modnodes[i].at(k) >-1) // Edge penalty in the sense of Nagata triangles
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
              if(accepted)
            {
          //      PropertyList *pl;
               // cout << prop << endl;
            //    sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]->props.GetPropertyList(prop,&pl);
              //  FaceVertexMesh *fv;
              //  fv = (FaceVertexMesh*) pl-> parent;
             //   cout << modprop[i]->index << endl;
           //     for(int lo=0;lo!=fv->vertices[modprop[i]->index].neighbors.size();lo++)
             //       if(fv->QualityAroundVertex(fv->vertices[modprop[i]->index].neighbors[lo]) < 1.) 
               //     {accepted=false; break;}
                  //  accepted=false;
            } // This should be caught out during normal verification!
            //}
            //}
            if(!accepted)
            {
              memcpy(modprop[i]->content,oldp,szProp*sizeof(double));
              modprop[i]->update();
              first=false;counter++;
             // TestReservation
           //  cout << "CUT" << endl;
           // for (int x=0;x!=cache.size();x++)  par->DeltaE(cache[x]);
            // cout << "END CUT" << endl;
            }
            else {     par->CacheWrite(cache);
                             #pragma omp critical
                          {
         if(first)
                          naccepted++; }
             
                 }  // Can Be Vectorized!}

            }
        }
        return naccepted / (double) modprop.size();
    }
          
};


class NormalMove : public MCMove
{
public:
    NormalMove()
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
        const int szProp=2; //modprop[0]->ptype->szProperty; // Only One Property Type in This example;
    
#pragma omp parallel for
        for (int i=0;i<modprop.size();i++)
        {
            double oldp[szProp]; 
            double newp[3];newp[0]=0;newp[1]=0;newp[2]=0;
            memcpy(oldp, modprop[i]->content,szProp*sizeof(double));
            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            GenNewP(newp,newp,3);
                
            NormalProperty* normal=(NormalProperty*) modprop[i];
            triple *nv=&((Vertex *) normal->parent)->norm;
          
            triple oldn=*nv;
            triple newn(newp);
            
            newn=newn+*nv;
            newn=newn/newn.abs();
                      bool allowed=true;
   
        //    if(*nv *newn < 0) allowed=false;
            duple newnorm=normang(newn);
            
            vector<double> collect(2,0);
            // update_coefficients.
            memcpy(modprop[i]->content,newnorm.p,szProp*sizeof(double));
             modprop[i]->update();
             
             if(*nv * (*nv) < 0.999) allowed=false;
	    // cout << *nv * (*nv) << endl;
        //    cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            for(int j=0;j!=modnodes[i].size();j++)
            {
                     allowed=allowed and allowNormal((NTriangle *) sim->n.col[modnodes[i].at(j)], nv,&oldn);
              if(!allowed) break;
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                triple t=((NTriangle *) sim->n.col[modnodes[i].at(j)])->mid;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
         
            }
            sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
              bool accepted = (sumE < 0) and allowed; // Cooling NOrmals
              //AcceptMove(sumE) and allowed;
            if(accepted)
            {
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            }
            if(!accepted)
            {
              memcpy(modprop[i]->content,oldp,szProp*sizeof(double));
              modprop[i]->update();
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

         
        }
        return naccepted / (double) modprop.size();
    }
 
protected:
  
private:
    
          duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    
   
};


class NormalFollowingMove : public MCMove
{
public:
    NormalFollowingMove()
    {

    numMoves=0;
    tryFrac=0.3;
    prop="Vertex";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=30;
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
        const int szProp=3;//doublemod[0].first->ptype->szProperty;
#pragma omp parallel for
        for (int i=0;i<doublemod.size();i++)
        {
            double oldp[szProp]; 
           // double newp[szProp];
            memcpy(oldp, doublemod[i].first->content,szProp*sizeof(double));
            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            double step=0, origin=0;
            GenNewP(&origin,&step,1);
            vector<double> collect(2,0);
            // Get Normal of Vertex.
            NormalProperty* normal=(NormalProperty*) doublemod[i].second;
            triple *nv=&((Vertex *) normal->parent)->norm;
            triple newp(oldp); newp=newp+(*nv*step);
            // update_coefficients.
            memcpy(doublemod[i].first->content,newp.p,szProp*sizeof(double));
            doublemod[i].first->update();
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
              allowed=allowed and allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]);
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
              doublemod[i].first->update();
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
  
 
};

bool CornerNormal(NTriangle *nt, triple *vertex, triple *normal)
{
    
    // Find the correct corner
     int mye=0,j=1,k=2;
       for(int i=0;i!=3;i++) { if(nt->v[i]==vertex) {mye=i; j=j%3;k=k%3; break; } j++; k++;}
      // Nagata organization:
      const double zf= 0.791;
      const double step=0.02;
      if(mye==0) *normal=nt->normal(step,0.5*step);
      else if(mye==2) *normal=nt->normal(1-step,step);
      else *normal=nt->normal(1-zf*step,zf*step);
      return true;

}


bool OuterNormal(NTriangle *nt, triple *vertex, triple *normal)
{
    
    // Find the correct corner
     int mye=0,j=1,k=2;
       for(int i=0;i!=3;i++) { if(nt->v[i]==vertex) {mye=i; j=j%3;k=k%3; break; } j++; k++;}
       triple n1=*nt->n.at(j);
       triple n2=*nt->n.at(k);
       triple d1=*(nt->v.at(mye))-*(nt->v.at(j));
       triple d2=*(nt->v.at(mye))-*(nt->v.at(k));
       *normal=n1/(d1*d1)+n2/(d2*d2);
      return true;
}
