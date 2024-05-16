


class AlexanderMove : public MCMove
{
public:
    AlexanderMove()
    {

    numMoves=0;
    tryFrac=0.5;
    prop="Edge";
    modprop.clear();
    modnodes.clear();
    busynodes.clear();
    maxstepsz=4;
    }
    
    // Returns Acceptance Ratio.
    double ExecuteMoves(Simulation* sim)
    {     int naccepted=0;
        // foreach_move//
         double var=maxstepsz*maxstepsz/9.;
     //   uniform_real_distribution<> dis(-maxstepsz, maxstepsz);
         //sim->Entities[0]->DeltaE(0); return -1;
        int szProp=modprop[0]->ptype->szProperty; // Only One Property Type in This example;
       if(maxstepsz < 0.1 ) maxstepsz=0.1;
        PropertyList *pl;
               // cout << prop << endl;
                par->props.GetPropertyList(prop,&pl);
                FaceVertexMesh *fvm;
                fvm = (FaceVertexMesh*) pl-> parent;
         
//                    FaceVertexMesh *fvm=(FaceVertexMesh*) modprop[0]->parent;
// #pragma omp parallel for
        for (int i=0;i<modprop.size();i++)
        {
           QuenchAll(fvm); // Fortest
    
            bool accepted=false;
            bool first=true;
            int counter =0;
            double sumE=0; // Penalty
            while(!accepted and counter < 2)
            {
              bool allowed=true;
          
            double sumE=0;
            Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
                
            Property* edge=(Property*) modprop[i];
            indexlist cache(modnodes[i].size());
            // First Get two faces:
            if(modprop[i]->deps.size()<2) { allowed=false; counter++; first=false;}
            else if(modprop[i]->deps[0]->ptype->collidable==false or modprop[i]->deps[1]->ptype->collidable==false) { allowed=false; counter++; first=false;} // First two dependencies _must be_ triangles.
            if(modprop[i]->block!=NULL) {allowed=false;counter++;first=false; } // Prevent Remeshing of Blocked Edge.
            if(allowed){ 
            NTriangle *nt1=(NTriangle *) modprop[i]->deps[0]->content;
            NTriangle *nt2=(NTriangle *) modprop[i]->deps[1]->content;
            //       (c)
            //       l|
            //  (a) Edge (b)
            //        |r
            //       (d)
            
            // -> update triangles & everything. faceids remain the same (!)
            int a=fvm->edges[edge->index].vertA;
            int b=fvm->edges[edge->index].vertB;
       //     cout << " a " << a << " b " << b << endl;
            // -> check connectivity. a and b will lose an edge.
            if(fvm->vertices[a].neighbors.size()<6 or fvm->vertices[b].neighbors.size()<6) {allowed=false;counter++;}
               
       //     cout << "PENALTY" << endl;
             //      cout << "DOING IT in EARNEST" << endl;
            // Step 1 Identify (a,b,c,d)
            int i1a,i1b,i1c;
            int i2a,i2b,i2d;
            int c=-1; int d=-1;
            int fci1=modprop[i]->deps[0]->index;
            int fci2=modprop[i]->deps[1]->index;
    //          cout << " fci1 " << fci1 << " fci2 " << fci2 << endl;
          if(allowed)
            {
            for(int j=0;j!=3;j++) 
              {
                 if(fvm->faces[fci1].vert[j]==a) i1a=j;
                 else if(fvm->faces[fci1].vert[j]==b) i1b=j;
                 else {c=fvm->faces[fci1].vert[j];i1c=j;}
                 if(fvm->faces[fci2].vert[j]==a) i2a=j;
                 else if(fvm->faces[fci2].vert[j]==b) i2b=j;
                 else {d=fvm->faces[fci2].vert[j]; i2d=j;}
              }
//            if(allowed and fvm->vertices[d].neighbors.size()>8 or fvm->vertices[c].neighbors.size()>8) {allowed=false;counter++;}
        
           if(fvm->vertices[c].neighbors.size()>7 or fvm->vertices[d].neighbors.size()>7) {allowed=false;counter++;}
          
       //   Check Hypothetical Triangle Tilt.
#define EDGE_TILT_LIMIT_X 3.141592/6.0
#ifdef EDGE_TILT_LIMIT_X
            triple cd=fvm->vertices[d].pos-fvm->vertices[c].pos;
            triple ad=fvm->vertices[d].pos-fvm->vertices[a].pos;
            triple bd=fvm->vertices[d].pos-fvm->vertices[b].pos;
            triple n1=ad%cd; n1=n1/n1.abs();
            triple n2=cd%bd; n2=n2/n2.abs();
            //cout << n1*n2 << endl;
            double tlt=n1*n2;
            if(allowed and tlt*tlt< EDGE_TILT_LIMIT_X) {allowed=false; counter++;}
#endif        
       //
            }
            if(allowed)
            {
                 sumE-=edgePenalty(nt1);
                 sumE-=edgePenalty(nt2);

              
            // Manipulate Triangles:
            nt1->v[i1b]=&fvm->vertices[d].pos; //FTF cd
            nt1->n[i1b]=&fvm->vertices[d].norm;
            nt2->v[i2a]=&fvm->vertices[c].pos;
            nt2->n[i2a]=&fvm->vertices[c].norm; 

            // Swap Properties
              nt2->updatec();
              nt1->updatec();
       
            modprop[i]->update();
           
           
            // Energy Evaluation.
            // Swap Complete: Now update
            allowed= allowAngle(nt1) and allowAngle(nt2);
            vector<double> collect(2,0);
            if(allowed)
            {
           // cout << modnodes[i].size() << endl;
            for(int j=0;j!=modnodes[i].size();j++)
            {
              if(modnodes[i].at(j)>-1)
              {
                SetAn((NTriangle *) sim->n.col[modnodes[i].at(j)],fvm,modnodes[i].at(j));
                //sumE+=edgePenalty((NTriangle *) sim->n.col[modnodes[i].at(j)]);
                int index=sim->n.col[modnodes[i].at(j)]->node_id;
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
              }
               else
              {
              int index=modnodes[i].at(j)*-1;
              sumE+=par->DExtraE(index);
              cache[j]=-index;
              }
            
                // TODO Evaluate Through-Space Couplings!
            }
            
                 sumE+=edgePenalty(nt1);
                 sumE+=edgePenalty(nt2);
            sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            }
            accepted=(sumE < 0.0) and allowed;
//          accepted= true;
         //     if(accepted)   
             // Only Now Check CollisionDetection
           //     for(int k=0;k!=modnodes[i].size();k++) 
           //        if(sim->cd.CheckCollision(modnodes[i].at(k))) {accepted=false; break;}
            if(!accepted)
            {
              
             

                // Restore
            // Manipulate Triangles:
            nt1->v[i1b]=&fvm->vertices[b].pos;
            nt1->n[i1b]=&fvm->vertices[b].norm;
          
            nt2->v[i2a]=&fvm->vertices[a].pos;
            nt2->n[i2a]=&fvm->vertices[a].norm; 
            // Unscramble Props.
              //nt2->updatec();
             // nt1->updatec();
         
              modprop[i]->update();
        
       
       
              first=false;counter++;
            }
            else {
            
            // if faceno l,r 
            // vertex b loses l, a loses r, c gains r, d gains l
            // l : v[b], n[b] now points to (d)
            // r : v[a], n[a] now points to (c)
            //Find new end iterator

            //Rebuild Neighbors
            // Only now do complicated Manipulations.
                    // a and b lose a neighbor      
#pragma omp critical
                {
//                                                                                                   cout << "Accept" << endl;
            std::vector<int>::iterator newEnd = std::remove(fvm->vertices[b].neighbors.begin(), fvm->vertices[b].neighbors.end(), a);
            fvm->vertices[b].neighbors.erase(newEnd, fvm->vertices[b].neighbors.end());
            newEnd = std::remove(fvm->vertices[a].neighbors.begin(), fvm->vertices[a].neighbors.end(), b) ;
            fvm->vertices[a].neighbors.erase(newEnd, fvm->vertices[a].neighbors.end());
       
             // d and c gain one neighbor.
            fvm->vertices[d].neighbors.push_back(c);
            fvm->vertices[c].neighbors.push_back(d);
            // Update Edges
          //  cout<< "EDGES " << fvm->edges.size() << endl;
            fvm->edges[edge->index].vertA=d; //c
            fvm->edges[edge->index].vertB=c;  //d
            
            // Update vl             
            // vertex b loses l, a loses r
            newEnd = std::remove(fvm->vl[b].begin(), fvm->vl[b].end(), fci1);
            fvm->vl[b].erase(newEnd, fvm->vl[b].end());
            newEnd = std::remove(fvm->vl[a].begin(), fvm->vl[a].end(), fci2);
            fvm->vl[a].erase(newEnd, fvm->vl[a].end());
            //  c gains r, d gains l
            fvm->vl[c].push_back(fci2);
            fvm->vl[d].push_back(fci1);
            // Update Faces
            fvm->faces[fci1].vert[i1b]=d;
            fvm->faces[fci2].vert[i2a]=c;
            int mark1=-1,mark2=-1;
            // Find Edge Indices
            for(int k=0;k!=3;k++)
            {
           //       (c)
            //    1/  l|  \(2) --> gets nt2 --> loses nt1
            //  (a) Edge (b)
            //    3\  |r  /(4) --> gets nt1 --> loses nt2
            //       (d)
         
            int edge1=fvm->faces[fci1].edges[k];
            int edge2=fvm->faces[fci2].edges[k];
            if(edge1!=modprop[i]->index)
              {
                if(fvm->edges[edge1].vertA==b or fvm->edges[edge1].vertB==b)
                  mark1=k;
              }
             if(edge2!=modprop[i]->index)
              {
                if(fvm->edges[edge2].vertA==a or fvm->edges[edge2].vertB==a)
                  mark2=k;
              }
            }
            if(mark1==-1 or mark2==-1) {cerr << "SOMETHING HORRIBLY WRONG!" << endl ; exit(0);}
            // Flip Edges of faces
            int edge1=fvm->faces[fci1].edges[mark1];
            int edge2=fvm->faces[fci2].edges[mark2];
            fvm->faces[fci1].edges[mark1]=edge2;
            fvm->faces[fci2].edges[mark2]=edge1;
            // Update Edge Dependencies
            Property* e1=pl->theProps[edge1];
            Property* e2=pl->theProps[edge2];
            int i1=0,i2=0;
            Property* e1d, *e2d; 
            if((NTriangle *) e1->deps[ 0]->content!=nt1 ) i1=1;
            if((NTriangle *) e1->deps[i1]->content!=nt1 ) cout << i1 << " NO!" << endl;

            if((NTriangle *) e2->deps[0 ]->content!=nt2 ) i2=1;
            if((NTriangle *) e2->deps[i2]->content!=nt2 ) cout << i2 << "NO!" << endl;
            e1d=e1->deps[i1];
            e2d=e2->deps[i2];
            
            e1->deps[i1]=e2d;  
            e2->deps[i2]=e1d; 
            // Actually Irreversible here.
            for(int j=0;j!=modnodes[i].size();j++)
            {
              if(modnodes[i].at(j)>-1)
              SetAn((NTriangle *) sim->n.col[modnodes[i].at(j)],fvm,modnodes[i].at(j));
            }
            
            // Now we should be consistent.
                }
        //    cout << "ACCEPTED" << endl;
                             par->CacheWrite(cache);
                             #pragma omp critical
                          {
         if(first)
                          naccepted++;}
                          
                 }}
                
            } }}  // Can Be Vectorized!}

            
          
        return true; //naccepted / (double) modprop.size();
        
    }
protected:
    void SetAn(NTriangle *nt, FaceVertexMesh *fvm, int findex, int nbadd=0)
    {
           double dv[3];
              for (int k = 0; k < 3; k ++)
        {
            dv[k]=fvm->vertices[fvm->faces[findex].vert[k]].neighbors.size();
        }
       nt->SetProperty("An0", (1-2/dv[0])*M_PI ); 
       nt->SetProperty("An1", (1-2/dv[2])*M_PI ); 
       nt->SetProperty("An2", (1-2/dv[1])*M_PI ); 
        return;
    }
    void QuenchAll(FaceVertexMesh *fvm)
    {
        for(int i=0; i!=modprop.size();i++)
            {
            NTriangle *nt1=(NTriangle *) modprop[i]->deps[0]->content;
            SetAn(nt1,fvm,modnodes[i].at(0));
            if(modprop[i]->deps.size()>1)
            if(modprop[i]->deps[1]->ptype->collidable) // To exclude edge energy.
            {NTriangle *nt2=(NTriangle *) modprop[i]->deps[1]->content;
             SetAn(nt2,fvm,modnodes[i].at(1));
            }
            }
        
    }
    void RedistributeProps(NTriangle *nt1, NTriangle nt2)
    {
        return;
    }
    };
    
