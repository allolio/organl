

      void MCMove::GenNewP(double* oldp, double* newp, int szprop)
    {
         double var=(maxstepsz*maxstepsz/9.); // Three Standard deviations is Maximum step
            normal_distribution<> dis(0.0, var);
       
        for (int i=0;i!=szprop;i++)
        {
            newp[i]=oldp[i]+dis(gen);
        }
    }
     bool MCMove::AcceptMove(double deltaE)
    {
//    if(deltaE>100) return false;
    uniform_real_distribution<> dis_u(0, 1);
    return (exp(-deltaE) > dis_u(gen));
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
        Deplist  *d=p->dependents(); // No Automatic Recursion here.
        for(int k=0;k!=d->size();k++)
        {
                       if (d->at(k)->ptype->collidable or d->at(k)->ptype->energetic) // In this case it is in the Neighbor List and its value is a collidable
                       {
                           int node_id=d->at(k)->index;
                           reservation.push_back(node_id);
                       }
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
      if (*norm * *nt->n[j] < 0.01 or *norm * *nt->n[k] < 0.01) return false;
      return true;
   }
   
       
 double MCMove::edgePenalty(NTriangle *n)
 {
   //  if(llim == NO_LIMIT and hlim==NO_LIMIT) return 0;
     double penalty=0;
     triple ta  = *(n->v[0])-*(n->v[1]);
     triple tb  = *(n->v[0])-*(n->v[2]);
     triple tc  = *(n->v[1])-*(n->v[2]);
    //         double area= (ta % tb).abs();
     // New: Failpenal //
             
      for (int i = 0; i < 3; i ++)
        {
            int j = (i + 1) % 3;
            triple d  = *(n->v[i])-*(n->v[j]);
             double a1=(*n->n[i]*d); double a2=(*n->n[j]*d);
             double fpow=(a1*a2);
             if(fpow > 0 ) 
             {
              
              double A0=n->GetProperty("Try_Area"); double kappa=n->GetProperty("kappa");
              double da=d.abs();
              penalty+=kpenal*kappa*10*A0*(a1*a1+a2*a2)/(da*da)+1;
             }
        }
        
     double tall,tball,tcall;
  //   tall=sqrt(ta*ta), tball=sqrt(tb*tb),tcall=sqrt(tc*tc);
    /// double amin=0.5*llim*(llim); double amax=0.5*hlim*(hlim);
 /*   if( tall < llim) penalty+=kpenal*(tall-llim)*(tall-llim);
    if( tall > hlim) penalty+=kpenal*(tall-hlim)*(tall-hlim);
    if( tball <llim) penalty+=kpenal*(tball-llim)*(tball-llim);
     if( tball > hlim) penalty+=kpenal*(tball-hlim)*(tball-hlim);
     if( tcall < llim) penalty+=kpenal*(tcall-llim)*(tcall-llim);
    if( tcall > hlim) penalty+=kpenal*(tcall-hlim)*(tcall-hlim);*/ 
  //  if( area < amin) penalty+=500*kpenal*(amin-area)*(amin-area);
  // if( area > amax) penalty+=kpenal*(amax-area)*(amax-area); 
    // Temporary: Angle Penalties
    double a1=acos(ta/ta.abs()*tb/tb.abs());
    double a2=acos(tb/tb.abs()*tc/tc.abs());
    double a3=M_PI-a1-a2;
  
 //   penalty+=kpenal*(M_PI/3-a1)*(M_PI/3-a1)/25/5;
 //   penalty+=kpenal*(M_PI/3-a2)*(M_PI/3-a2)/25/5;
 //   penalty+=kpenal*(M_PI/3-a3)*(M_PI/3-a3)/25/5;
    if(a1*a1 < (M_PI/5)*(M_PI/5)) penalty+=kpenal*(a1-M_PI/5)*(a1-M_PI/5)/100;
    if(a1*a1 > (M_PI/2.5)*(M_PI/2.5) )  penalty+=kpenal*(M_PI/2.5-a1)*(M_PI/2.5-a1)/100;
    if(a2*a2 < (M_PI/5)*(M_PI/5)) penalty+=kpenal*(a2-M_PI/5)*(a2-M_PI/5)/100;
    if(a2*a2 > (M_PI/2.5)*(M_PI/2.5) )  penalty+=kpenal*(M_PI/2.5-a2)*(M_PI/2.5-a2)/100;
    if(a3*a3 < (M_PI/5)*(M_PI/5)) penalty+=kpenal*(a3-M_PI/5)*(a3-M_PI/5)/100;
    if(a3*a3 > (M_PI/2.5)*(M_PI/2.5) )  penalty+=kpenal*(M_PI/2.5-a3)*(M_PI/2.5-a3)/100;
     // penalty+=kpenal*(M_PI/3-a2)*(M_PI/3-a2)/25/5;
     // penalty+=kpenal*(M_PI/3-a3)*(M_PI/3-a3)/25/5;
    return penalty;
 }
    
    bool MCMove::allowAngle(NTriangle* nagata)
{
    //cout << "Pen\n";
    double cosAlpha0 = 0.98078528;      // cos(pi/16)      // cos(7.5 deg)
    double step =1; // Originally 0.1             // step for evaluation within nagata triangle
    bool penalty = true;
    double scale = 200.0;
    triple v;
    triple a;
    triple b;
    triple ta;
    triple tb;
    double cosAlpha;

    v = nagata->eval(0, 0);
    a = nagata->eval(step, 0);
    b = nagata->eval(step, step);
    ta = (a - v);
    tb = (b - v);
    triple tc;
    
    //if(hlim != NO_LIMIT or llim != NO_LIMIT) 
    //    if(!edgelimit(nagata)) return false;
    
    cosAlpha = (ta * tb) / (ta.abs() * tb.abs());
    //cout << cosAlpha << "\n";
    if ( cosAlpha > cosAlpha0)
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
    //cout << cosAlpha << "\n";
    if ( cosAlpha > cosAlpha0)
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
    //cout << cosAlpha << "\n";
    if ( cosAlpha > cosAlpha0)
    {
        return false;
//        penalty += exp(scale*(cosAlpha - cosAlpha0));
        // cout << "Penalty: angle too low: " << cosAlpha << "\n";
    }

    //if (penalty > 0)
        //cout << "\t angle penatly: " << penalty << "\n";


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
        int numProps;
              modprop.clear(); modnodes.clear(); busynodes.clear();
         for(int i=0;i!=sim->Entities.size();i++)
        {
       
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
         //   cout << "MOVE Supported" << endl;
            PropertyList *pl;
                               Property *p;
            sim->Entities[i]->props.GetPropertyList(prop,&pl);
            numProps=pl->theProps.size();
            uniform_real_distribution<> dis_u(0, 1);
                    indexlist reservation; reservation.clear();
                    int select=(int) (dis_u(gen)*numProps); // select random index.
                    pl->pofs(select,&p);
                                        BuildReservation(p,reservation);

                        WriteReservation(p,reservation); //Not
 
            }
        
        }
        cycle++;
        if(cycle==numProps) cycle=0;
        return cycle;
        
        if(cycle==-1)
        {
            cout << "preparing precomputed steps" << endl;
        modprop.clear(); modnodes.clear(); busynodes.clear();
        for(int i=0;i!=sim->Entities.size();i++)
        {
            if(CheckEntitySupportsMove(sim->Entities[i]))
            {
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
            for(int j=0;j<(int) 1;j++) // CAN BE VECTORIZED
               {
                    Property *p;
                    indexlist reservation; reservation.clear();
                    int select=(int) (dis_u(gen)*leftover.size()); // select random index.
                //    cout << "GET SELECT " << select << " of" << numProps<< endl;
                    pl->pofs(leftover[select],&p);
                //    cout << "BUILDING at" << select << endl;
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
               sort(removal.begin(), removal.end(), greater<int>());
               // leftover is sorted always.
               std::deque<int>::iterator leftend=leftover.end();
               std::deque<int>::iterator kill; bool killed=false;
               cout << "inside " << leftover.size() << endl; 
               for(int m=0;m!=removal.size();m++)
               {
               //    cout << "Removal " << m << " " << removal[m] << endl;
                   while(leftend>=leftover.begin() and killed==false and leftover.size()>0)
                   {
                       if(*leftend==removal[m])
                       {
                        killed=true;
                        kill=leftend;
                        }
                   
                       leftend--;
                       if(killed) { leftover.erase(kill);  }
                   }
                   leftend++;
                   if(killed==false) cout << "NOT KILLED" << endl;
                   killed=false;
               }
               cout << "remove" << removal.size() << endl;
               cout << "leftover" << leftover.size() << endl;
               removal.clear();
      //         cout << "leftover" ;
//               for(int m=0;m!=leftover.size();m++) cout << leftover[m] << " " ;
       //        cout << endl;   
         //      cout << "end leftover" << endl;
               busynodes.clear();
                cyclemp.push_back(modprop);
           cout << "MPS " << modprop.size() << endl;
           cyclemnd.push_back(modnodes); modnodes.clear(); modprop.clear();

            }
            
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
    maxstepsz=0.01;
    }
    
    /*Schema*/
        /*  1. Check if Entity Supports Move */
    
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move
//        cout << "MPS" << modprop.size()<< " " << modprop[0]->index << endl; 

       //  cout << "MODPROP " << modprop.size() << " " << modprop[0]->index << endl;
       //  cout << "MODNODES " << modnodes[0].size() << endl;
        //sim->Entities[0]->DeltaE(0); return -1;
        int szProp=modprop[0]->ptype->szProperty; // Only One Property Type in This example;
#pragma omp parallel for
        for (int i=0;i<modprop.size();i++)
        {
            double oldp[szProp]; 
            double newp[szProp];
            memcpy(oldp, modprop[i]->content,szProp*sizeof(double));
            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
            // Generate Step;
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE-=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
            //
            GenNewP(oldp,newp,szProp);
            vector<double> collect(2,0);
            // update_coefficients.
            memcpy(modprop[i]->content,newp,szProp*sizeof(double));
             modprop[i]->update();
          //  cout << "DEPENDENCY UPDATED" << endl;
            indexlist cache(modnodes[i].size());
            bool allowed=true;
            /* Check for singular edge
             * */
            vector<triple> fnorm(modnodes[i].size());
           
            for(int j=0;j!=modnodes[i].size();j++)
            {
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                // TODO Evaluate Through-Space Couplings!
                // FIXME Penalty and Reversibility
              allowed=allowAngle((NTriangle *) sim->n.col[modnodes[i].at(j)]);
              if(!allowed) break;
            }
            bool accepted=false;
            if(allowed)
            { sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            for(int j=0;j!=modnodes[i].size();j++)
            {
                 sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
            }
               accepted = AcceptMove(sumE);
            }
            if(accepted)
            {
                // Only now check Singular face you can move this at will
              for(int j=0;j!=modnodes[i].size();j++)
            {
                fnorm[j]=(*((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[0]-*(((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[1])) % (*((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[0]-*(((NTriangle *) sim->n.col[modnodes[i].at(j)])->v[2]));
                fnorm[j]=fnorm[j]/fnorm[j].abs();
               // if(fnorm[j]* *((NTriangle *) sim->n.col[modnodes[i].at(j)])->n[1] < 0) fnorm[j]=fnorm[j]*-1;
                for(int k=0;k<j;k++)
                    if((fnorm[k]*fnorm[j])*(fnorm[k]*fnorm[j])<EDGE_TILT_LIMIT) {accepted=false; break;}
                if(!accepted) break;
            }
            }
            if(accepted)   
             // Only Now Check CollisionDetection
                for(int k=0;k!=modnodes[i].size();k++) 
                   if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
              if(accepted)
            {
                PropertyList *pl;
               // cout << prop << endl;
                sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]->props.GetPropertyList(prop,&pl);
                FaceVertexMesh *fv;
                fv = (FaceVertexMesh*) pl-> parent;
             //   cout << modprop[i]->index << endl;
                for(int lo=0;lo!=fv->vertices[modprop[i]->index].neighbors.size();lo++)
                    if(fv->QualityAroundVertex(fv->vertices[modprop[i]->index].neighbors[lo]) < 1.) 
                    {accepted=false; /*break;*/}
                  //  accepted=false;
            }
            //}
            //}
            if(!accepted)
            {
              memcpy(modprop[i]->content,oldp,szProp*sizeof(double));
              modprop[i]->update();
             // TestReservation
           //  cout << "CUT" << endl;
           // for (int x=0;x!=cache.size();x++)  par->DeltaE(cache[x]);
            // cout << "END CUT" << endl;
            }
            else {     par->CacheWrite(cache);
                             #pragma omp critical
                          {
         
                          naccepted++;}
             
                 }  // Can Be Vectorized!}

           
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
        int szProp=modprop[0]->ptype->szProperty; // Only One Property Type in This example;
    
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
    
    
  
       //j=j%3;k=k%3;
       //tuple v1=(*nt->v[j])%(*nt->v[k])
       
       // Very Ugly
       
           /*
    
     *
     * checks whether the normal n is not pointing to the other side of the triangle
     * then nOld
     *
     * */

/*       bool Nagata::checkNormals(const Face & face, int vertexid, triple & n, triple & nOld)
{

    int i = ( (face.vert[0] == vertexid) ? face.vert[1] : face.vert[0] );
    int j = ( (face.vert[2] == vertexid) ? face.vert[1] : face.vert[2] );

    triple y = (fv->vertices[i].pos - fv->vertices[vertexid].pos)
        % (fv->vertices[j].pos - fv->vertices[vertexid].pos);

    return ( (y*nOld) * (y*n) > 0);
   }*/


       
   //    return ( y > 0);
   
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
        int szProp=doublemod[0].first->ptype->szProperty;
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
