#define DIRTY_HACK_MAX_LIP 500
#define DIRTY_HACK_MIN_LIP  0.00 // 1
#define ADS_APL -1
struct Lipid
{
double APL;
double c0;
double kappa;
double charge;
double d; // Distance
double kappa_g; // Gauss Modulus
double mu; // Chemical Potential
double B;
string name;
int type;
};

class LipidDB
{
public:
    bool LoadLipids(const string &lipidfile);
    bool DisplayLipids();
    bool HasLipid(const string &lipname);
    virtual int ProcessFace(NTriangle *nt); // To be overwritten for special lipids
    bool GetNLipids(int *lipnum, int szlip);
    bool GetNLipidName(int num, bool side,string &lipname);
    bool IsNLipidP(int num, bool side);
    double GetNLipidAPL(int num, bool side);

    bool SimplifyLipids(indexlist &indic, bool side);
protected:
    std::vector <  Lipid > thelipidsdown;
    std::vector < Lipid > thelipids;  
};

bool LipidDB::GetNLipids(int *lipnum, int szlip)
{
    if(szlip!=2) return false;
    lipnum[0]= thelipids.size();
    lipnum[1]= thelipidsdown.size();
    return true;
}

bool LipidDB::SimplifyLipids(indexlist &indic, bool side)
{
    std::vector <  Lipid > replacement;
    std::vector <  Lipid > *orig;
    if(side) orig=&thelipidsdown;
    else orig=&thelipids;
    replacement.resize(indic.size());
    for(int i=0;i!=indic.size();i++)
    {
      if(indic[i]>= orig->size()) return false;
      replacement[i]=orig->at(indic[i]);  
    }
    *orig=replacement;
    return true;
}

bool LipidDB::GetNLipidName(int num, bool side, string &lipidname)
{
    if(side)
    {if(thelipidsdown.size()<=num) {return false; }
        else lipidname=thelipidsdown[num].name; return true;}
    if (thelipids.size()<=num) return false; 
    lipidname=thelipids[num].name;
    return true;
}


bool LipidDB::IsNLipidP(int num, bool side)
{
    if(side)
    return (thelipidsdown[num].APL==ADS_APL);
    else return (thelipids[num].APL==ADS_APL);
}

double LipidDB::GetNLipidAPL(int num, bool side)
{
    if(side)
    return thelipidsdown[num].APL;
    else return thelipids[num].APL;
}

int LipidDB::ProcessFace(NTriangle* nt)
{
//    cout << "PROCESSING FACE" << endl;
    double ntot=0, ntu[2];ntu[0]=0, ntu[1]=0;
    double c0=0, c0u[2];c0u[0]=0; c0u[1]=0; // total up and down
    double kapparec=0, kru[2]; kru[0]=0, kru[1]=0; 
    double A0=0, A0u[2]; A0u[0]=0; A0u[1]=0; 
    double charge=0; double kg=0;
    double S=0, Su[2]; Su[0]=0; Su[1]=0;
    //double ntot=nt->GetProperty("NTotalLipids"); 
    bool down=false;   
    bool bypass=false;
    vector <Lipid> *lin=&thelipids;
    while(true)
    {
        double adc=0;
                bypass=false;
    for(int i=0;i!=lin->size();i++)
        {
         double un;
         if(!down) un=nt->GetProperty(lin->at(i).name); 
         else un=nt->GetProperty("_"+lin->at(i).name);
    // Prefix >up  <down                
        ntu[down]+=un;
        double weight=(un*lin->at(i).APL);
        if(lin->at(i).APL==ADS_APL and un >0.9)
          {
          kg=lin->at(i).kappa_g;
          bypass=true;
          adc=lin->at(i).c0*un*(-1+2*!down);
          weight=0;
          }
          else
          {
        if(weight>0){
        kru[down]+=1./lin->at(i).kappa*weight;
        c0u[down]+=lin->at(i).c0*weight*(-1+2*!down);
        Su[down]-=un*lin->at(i).APL*log(weight);}
        A0u[down]+=lin->at(i).APL*un;
          }
              
       // else
       // {if(un>0) {bypassc=true; bypassc0=(*found).second.c0;}}
        charge+=lin->at(i).charge*un;
        
       }
       if(bypass==true)
       {
           
          c0u[down]=A0u[down]*adc;
          nt->SetProperty("kappag",kg);
       }
       else
       {
        nt->SetProperty("kappag",0.0); // Mixing terms for K_g not defined.
       }
       if(down==false)
       {down=true;
        lin=&thelipidsdown;
       }
       
       else break;
    }  
    //if(!bypassc)
    c0=0.5*(c0u[0]/A0u[0]+c0u[1]/A0u[1]); // 0.5*(c0u[0]/=A0u[0]+c0u[1]/=A0u[1])
    kapparec=1/(kru[0]/(A0u[0]))+1/(kru[1]/(A0u[1]));
//    cout << "["<< c0 << " " <<  c0u[0]/A0u[0] << " "<< c0u[1]/A0u[1] << "]"  << endl;
    Su[0]/=A0u[0], Su[1]/=A0u[1];
    Su[0]+=log(A0u[0]);Su[1]+=log(A0u[1]);
    S=Su[0]*ntu[0]+Su[1]*ntu[1]; A0=0.5*(A0u[0]+A0u[1]);
    if(S!=S) {cerr << "Entropy Crash "<<  Su[0]<< " " << ntu[0]<< "  " << A0u[0]  << " :  " << Su[1]*ntu[1]<< endl; S=0;}
    nt->SetProperty("c0",c0);
    nt->SetProperty("kappa",kapparec);
    nt->SetProperty("A0",A0);
    nt->SetProperty("Entropy",S); // Easy to Calc on the fly
    nt->SetProperty("Charge",charge);
//    nt->SetProperty("kappag",kg);
    //cout << "["<< 1/kapparec << " " <<  c0 << " "<< A0 << "]"  << endl;

    //TODO Add Area Modulus
    if((ntu[1]+ntu[0])==0) { nt->SetProperty("kappa",0); nt->SetProperty("A0",0); nt->SetProperty("Entropy",0);}
    return 1;
}

bool LipidDB::HasLipid(const std::string& lipname)
{
    for(int i=0;i!=thelipids.size();i++)
        if(thelipids[i].name==lipname) return true;
    for(int i=0;i!=thelipidsdown.size();i++)
        if(thelipidsdown[i].name==lipname) return true;
    return false;
}


bool LipidDB::DisplayLipids()
{
    cout << "TOP" << endl;
    for(int j=0;j!=thelipids.size();j++)
        {
            cout << j << " " << thelipids[j].name << " " << thelipids[j].c0 << " " << thelipids[j].kappa << " " << thelipids[j].APL  << " " << thelipids[j].d  << " C " << thelipids[j].charge << " " << thelipids[j].kappa_g<<  " " << thelipids[j].mu << " " << thelipids[j].B << endl;
        }
    cout << "BOTTOM" << endl;
    for(int j=0;j!=thelipidsdown.size();j++)
        {
            cout << j << " " << thelipidsdown[j].name << " " << thelipidsdown[j].c0 << " " << thelipidsdown[j].kappa << " " << thelipidsdown[j].APL  << " " << thelipidsdown[j].d  << " C " << thelipidsdown[j].charge << " " << thelipidsdown[j].kappa_g<<  " " << thelipidsdown[j].mu << " " << thelipidsdown[j].B << endl;
        }
    return true;
}

bool LipidDB::LoadLipids(const string &lipidfile)
{

    ifstream hfile;
    hfile.open(lipidfile.c_str());
    if(hfile.fail()) return false;
    string current;
    int Vcounter = 0;
    int Ncounter = 0;
    while ( !hfile.eof() )
    {
        vector<string> line;
        getline(hfile, current);
        int htokens=Tokenize(current,&line, " \t\r");
        if (line.size()==0 or line[0][0] == '#')
            continue;
        Lipid lip;
        int i=0;
        while(i<line.size())
        {
        switch(i)
        {
            case 0: 
                lip.name=line[0];
                break;
            case 1: 
                lip.c0=atof(line[1].c_str());
                break;
            case 2:
                lip.kappa=atof(line[2].c_str());
                break;
            case 3:
                lip.APL=atof(line[3].c_str());
                break;
            case 4:
                lip.d=atof(line[4].c_str());
                break;
            case 5:
                lip.charge=atof(line[5].c_str());
                break;
            case 6:
                lip.kappa_g=atof(line[6].c_str());
                break;
            case 7:
                lip.mu=atof(line[7].c_str());
                break;
            case 8:
                lip.B=atof(line[8].c_str());
                break;
        }
        i++;
        //lipids[line[0]]=
        }
        // Insert defaults
        while(i<8)
        {
        switch(i)
        {
            case 0: 
                lip.name="POPC";
            case 1: 
                lip.c0=0.0;
            case 2:
                lip.kappa=15;
            case 3:
                lip.APL=0.64;
            case 4:
                lip.d=1.2;
            case 5:
                lip.charge=0.0;
            case 6:
                lip.kappa_g=-0.8*lip.kappa;
            case 7:
                lip.mu=0;
            case 8:
                lip.B=-1;
        }
        i++;
        }
        thelipids.push_back(lip);
        thelipidsdown.push_back(lip);
    }
    return true;
}

class MixFaceEnerProperty;

class MixFaceEnerProperties: public PropertyList
{
public:
    MixFaceEnerProperties()
    {   name="FaceEL";
        collidable=true;
        energetic=true;
        szProperty=sizeof(void*);
    }
    virtual bool genpropertiesfromStructure(void *parent);
};

class MixFaceEnerProperty: public Property
{
public:
    MixFaceEnerProperties *ptype;
    bool update(bool deep=true) {
    //    if(block!=NULL) block->applyBlock(content); NO BLOCKS
      //  Vertex *v=(Vertex*) parent;
      
        NTriangle *nag=(NTriangle*) content;
        LipidDB* ldb=(LipidDB *) parent;
        state=ldb->ProcessFace(nag);
        // Composition to Properties.
    //    if(deep) NO DEPS
    //        for(int i=0; i!=deps.size(); i++) deps[i]->update();
        return nag->updatec();
    }
};

 bool MixFaceEnerProperties::genpropertiesfromStructure(void* parent)
    {
        double** init= (double **) parent;
        vector< NTriangle > *tri=(vector < NTriangle > *) init[0];
        theProps.resize(tri->size());
        for(int i=0;i!=tri->size();i++)
        {
            theProps[i]=new MixFaceEnerProperty;
            theProps[i]->block=NULL;
            theProps[i]->parent=init[1];
            theProps[i]->deps.clear();
            theProps[i]->index=i;
            theProps[i]->state=0;
            theProps[i]->ptype=this;
            theProps[i]->content=(double *) &tri->at(i);
            theProps[i]->update(); // Generate Things Immediately.
        }
        return true;
    }
    

class MixedMembrane : public Membrane
{
public:
     MixedMembrane(Energy &en): Membrane(en)
    {
//        MembE=new ConstraintHelfrichEnergy;
    };
    bool SetLipids(LipidDB *lipdDB);
    bool RestoreDependencies();
protected:
    bool PopulateProps()
{
    cout << "MIXED LIPID MEMBRANE" << endl;
    double *imfep[2];
    imfep[0]=(double *) &tri;
    imfep[1]=(double *) lip;
    normp.genpropertiesfromStructure(fv);
    vertp.genpropertiesfromStructure(fv);
    cout << "VERTICES" << vertp.theProps.size() << endl;
    CullLipids();
    lip->DisplayLipids();
    membp.genpropertiesfromStructure(this);
    edgep.genpropertiesfromStructure(fv);
    mfep.genpropertiesfromStructure(imfep);
    props.AddPropertyType(normp);
    props.AddPropertyType(vertp);
    props.AddPropertyType(membp);
    props.AddPropertyType(mfep);
    props.AddPropertyType(edgep);
    props.AddPropertyType(eep);
   vector<string> names;
   props.AvailProp(names);
   for(int i=0;i!=names.size();i++)
   {cout << names[i] << " ";
    PropertyList *pl;
    props.GetPropertyList(names[i],&pl);
    cout << pl->theProps.size() << endl;
   }
      /*
    * Connect face energies to vertices
    */
 
   for(int i=0;i!=fv->vertices.size();i++)
   {
    Property *n=normp.theProps[i];
    Property *v=vertp.theProps[i];
    for(int j=0;j!=fv->vl[i].size();j++)
     {
             Property *f=mfep.theProps[fv->vl[i].at(j)];
             n->deps.push_back(f);
             v->deps.push_back(f);
//             f->parent=this;
     }
   }
      /*
    * Connect face energies to Edges
    */
    for(int i=0;i!=fv->edges.size();i++)
   {
   Property *n=edgep.theProps[i];
    for(int j=0;j!=fv->vl[fv->edges[i].vertA].size();j++)
      {
          int fce=fv->vl[fv->edges[i].vertA][j];
          for(int k=0;k!=fv->faces[fce].vert.size();k++)
          {
              if(fv->faces[fce].vert[k]==fv->edges[i].vertB)
              {
                 Property *f=mfep.theProps[fce];
                 n->deps.push_back(f);
                 f->update(); // INitialize
              }
          }
      }
    }
  // CullLipids();
   return true;
    
 //   cout << "NORMPROPS " <<  normp.theProps.size() << endl;
 //   cout << "VPROPS " <<  normp.theProps.size() << endl;
}
    MixFaceEnerProperties mfep;
    LipidDB* lip;
    bool CullLipids() // This function exists to restrict lipids which are on only one side.
    {
        vector<double> counters;
        vector<double> dcounters;
        int lud[2];
      
        lip->GetNLipids(lud,2);
        counters.resize(lud[0]);
        dcounters.resize(lud[1]);
        for(int i=0;i!=tri.size();i++)
        {
          for(int j=0;j!=counters.size();j++)
          {
              string name;
              lip->GetNLipidName(j,0,name);
              counters[j]+=tri[i].GetProperty(name);
          }
                  for(int j=0;j!=dcounters.size();j++)
          {
              string name;
              lip->GetNLipidName(j,1,name);
              dcounters[j]+=tri[i].GetProperty("_"+name);
          }
        }
       indexlist s1,s2;
       for(int i=0;i!=counters.size();i++) 
           if(counters[i]!=0) s1.push_back(i);
       for(int i=0;i!=dcounters.size();i++) 
           if(dcounters[i]!=0) s2.push_back(i);
  
        lip->SimplifyLipids(s1,0);
        lip->SimplifyLipids(s2,1);
  
      return true;
    }
};

bool MixedMembrane::RestoreDependencies()
{
          PropertyList *pl,*pl2,*fep;
          string prop="Vertex";
          this->props.GetPropertyList(prop,&pl);
          prop="Normal";
          this->props.GetPropertyList(prop,&pl2);
          prop="FaceEL";
          this->props.GetPropertyList(prop,&fep);
          cout << "FEPS" << fep->theProps.size();
          FaceVertexMesh *fv;
          fv = (FaceVertexMesh*) pl-> parent;
          
     for(int i=0;i!=fv->vertices.size();i++)
   {
    Property *n=pl2->theProps[i];
    Property *v=pl->theProps[i];
    n->deps.clear();
    Deplist cache; cache.clear();
    for(int j=0;j!=v->deps.size();j++)
    {
        if(v->deps[j]->ptype->name!="FaceEL") cache.push_back(v->deps[j]);
    }
    //v->deps.clear();
    v->deps=cache;

    for(int j=0;j!=fv->vl[i].size();j++)
     {
             Property *f=fep->theProps[fv->vl[i].at(j)];
             n->deps.push_back(f);
             v->deps.push_back(f);
     //        f->parent=this;
     }
   }

        // Restore Custom Energies
   for(int i=1;i<extraterms.size();i++)
       extraterms[i]->restore();


  return true;
}


bool MixedMembrane::SetLipids(LipidDB *lipDB)
{
    if(lipDB!=NULL) lip=lipDB;
    return true;
}

class LipidMixMove : public MCMove
{
public:
    LipidMixMove()
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
        LipidDB *ldb=(LipidDB*) modprop[0]->deps[0]->parent;
        int lipno[2];
        ldb->GetNLipids(lipno,2);
     
 #pragma omp parallel for
        for (int i=0;i<modprop.size();i++)
        {
            bool accepted=false;
            bool first=true;
            int counter =0;
            double sumE=0; // Penalty
            while(!accepted and counter < 20)
            {
              bool allowed=true;
          
            double sumE=0;
           // Entity *par=sim->Entities[sim->n.parent_entity[modnodes[i].at(0)]]; // One Move one Entity
                
            Property* edge=(Property*) modprop[i];
            indexlist cache(modnodes[i].size());
            // First Get two faces:
            if(modprop[i]->deps.size()<2 or modprop[i]->block!=NULL) { allowed=false; counter++; first=false; // Prevent swapping across blocked or dangling Edges.
	    }
	                // Note that this is inflexible because the first dependencies are presumed to be faces.
	        else if(modprop[i]->deps[0]->ptype->name!="FaceEL" or modprop[i]->deps[1]->ptype->name!="FaceEL") {allowed=false;counter++;first=false;}
            if(allowed)
            {
            
            NTriangle *nt1=(NTriangle *) modprop[i]->deps[0]->content;
            NTriangle *nt2=(NTriangle *) modprop[i]->deps[1]->content;
            // Now mix the lipids;
                 sumE-=Entropy_Penalty(nt1);
                 sumE-=Entropy_Penalty(nt2);
      
            uniform_real_distribution<> l1d(0.0,1);
            double move1=l1d(gen);
            double move2=l1d(gen);
           //  move1=move2; // This is to enforce balance.
            int lip1u=l1d(gen)*lipno[0];
            int lip1d=l1d(gen)*lipno[1];
            int lip2u=l1d(gen)*lipno[0];
            int lip2d=l1d(gen)*lipno[1];
            string l1uname,l1dname,l2dname,l2uname;
            ldb->GetNLipidName(lip1u,false,l1uname);
            ldb->GetNLipidName(lip1d,true,l1dname);
            ldb->GetNLipidName(lip2u,false,l2uname);
            ldb->GetNLipidName(lip2d,true,l2dname);
            bool l1up=false, l1dp=false, l2up=false,l2dp=false;
            double l1uap, l1dap, l2uap,l2dap;
            
            // Check if is Protein
            l1uap=ldb->GetNLipidAPL(lip1u,false);
            l1dap=ldb->GetNLipidAPL(lip1d,true);
            l2uap=ldb->GetNLipidAPL(lip2u,false);
            l2dap=ldb->GetNLipidAPL(lip2d,true);
            if(l1uap==ADS_APL) l1up=true;
            if(l2uap==ADS_APL) l2up=true;
            if(l1dap==ADS_APL) l1dp=true;
            if(l2dap==ADS_APL) l2dp=true;
            
            
            l1dname="_"+l1dname;
            l2dname="_"+l2dname;
           // cout << "LIPIu " << lip1u << " " << l1uname << " LIPII u "<< lip2u << " " << l2uname <<endl;
           // cout << "LIPId " << lip1d << " " << l1dname << " LIPII d "<< lip2d << " " << l2dname <<endl;
    //        cout << move1 << " "  << move2 << endl;
            double l1upop=nt1->GetProperty(l1uname);
            double l1dpop=nt1->GetProperty(l1dname);
            double l2upop=nt2->GetProperty(l2uname);
            double l2dpop=nt2->GetProperty(l2dname);
            move1*=maxstepsz; 
           // move1= int(move1*(maxstepsz+1.0));// +0.5 why ?
            move2*=maxstepsz; 
            if(l1dp or l2dp) {move1=0;move2=0;}
    ;
           
            if(l1up) move1=1.0; //int(move1);
            if(l2up) move2=1.0; //int(move2);
            if(l1up == false and l2up==false){ move2=move1*l1uap/l2uap;} // Comment this out for area fluctuations
//            else move2=1;
           // move2=int((maxstepsz+1.0)*move2);
           //  cout << move1 << " "  << move2 << endl;
            double old1lu=nt1->GetProperty(l2uname); 
            double old1ld=nt1->GetProperty(l2dname); 
            double old2lu=nt2->GetProperty(l1uname); 
            double old2ld=nt2->GetProperty(l1dname);
           //  cout << " l1uap " << l1uap << " l1dap " << l1dap << " l2uap " << l2uap << " l2dap " <<  l2dap << endl; 
           // cout << " l1up " << l1upop << " l1d " << l1dpop << " l2u" << l2upop << " d " <<  l2dpop << endl; 
           // cout << " ol1up " << old1lu << " ol1d " << old1ld << " ol2u" << old2lu << " d " <<  old2ld << endl; 
            double dmove1=move1*l1uap/l1dap;
            double dmove2=move2*l2uap/l2dap;
            if(l1up or l2up) {dmove1=0;dmove2=0;}
          
            if(l1dp) {dmove1=1.0;}
            if(l2dp) {dmove2=1.0;}
            //cout << " move1 " << move1 << " dmove1 "<< dmove1 << endl;
            if(l1upop-move1 < DIRTY_HACK_MIN_LIP or l2upop-move2 < DIRTY_HACK_MIN_LIP or l1dpop-dmove1 < DIRTY_HACK_MIN_LIP or l2dpop-dmove2 < DIRTY_HACK_MIN_LIP) {allowed=false;first=false;counter++; /*cout << "PHAIL" << endl;*/ }
            else if(move1+old2lu > DIRTY_HACK_MAX_LIP or dmove1+old2ld > DIRTY_HACK_MAX_LIP or move2+old1lu > DIRTY_HACK_MAX_LIP or dmove2+old1ld > DIRTY_HACK_MAX_LIP ) {allowed=false;first=false;counter++; }
            else{
           // cout <<"MOVE " << move1 <<" "<< l1pop  <<" " << move2 << " " << l2pop << " L1 " << l1name << " " << lip1 << " L2 " << l2name << " " <<  lip2 << " " <<  endl; 
           // cout << "TRYING" << endl;
           // cout << " l1POP nt1 " << l1pop << "Move1 " << move1 << " Move2 " << move2 << " old21 " << old21 << " l2pop " << l2pop << " old12 "  << old12 << endl;  
            if(l1uname!=l2uname)
            {
     
            nt1->SetProperty(l1uname,l1upop-move1);
            nt2->SetProperty(l1uname,old2lu+move1); //
            nt1->SetProperty(l2uname,old1lu+move2); //
            nt2->SetProperty(l2uname,l2upop-move2);
            }
            else if(move1!=move2)
            {
            nt1->SetProperty(l1uname,l1upop-(move1-move2));
            nt2->SetProperty(l1uname,l2upop+(move1-move2));
            }

            if(l1dname!=l2dname)
            {
            nt1->SetProperty(l1dname,l1dpop-dmove1);
            nt2->SetProperty(l1dname,old2ld+dmove1);
            nt1->SetProperty(l2dname,old1ld+dmove2); //
            nt2->SetProperty(l2dname,l2dpop-dmove2);
            }
            else if(move1!=move2)
            {
            nt1->SetProperty(l1dname,l1dpop-(dmove1-dmove2));
            nt2->SetProperty(l1dname,l2dpop+(dmove1-dmove2));
            }

            // Swap Complete: Now update
            modprop[i]->update(true); // Call explicit deep update. Due to the direct call this will call the default update function!
                                      // Not the EdgeP update function.
            
            //cout << nt1->GetProperty("POPC") <<" " << nt1->GetProperty("DOPC") << " " << nt1->GetProperty("DOPE") << endl; 
            //cout << nt2->GetProperty("POPC") <<" " << nt2->GetProperty("DOPC") << " " << nt2->GetProperty("DOPE") << endl; 

            sumE+=Entropy_Penalty(nt1);
            sumE+=Entropy_Penalty(nt2);

            if(nt1->GetProperty("A0")==0 or nt2->GetProperty("A0")==0) allowed=false; //Not allowed to empty triangle.
                        vector<double> collect(2,0);
            if(allowed)
            for(int j=0;j!=modnodes[i].size();j++)
            {
                int index=modnodes[i].at(j);
                if(index>=0)
                {
                index=sim->n.col[modnodes[i].at(j)]->node_id;
                
                sumE+=par->DeltaE(index,&collect); //Needs to collect delta of constraints for whole move;
                cache[j]=index;
                }
                else 
                {
                index=index*-1;
                sumE+=par->DExtraE(index);
                cache[j]=-index;
                }
                // TODO Evaluate Through-Space Couplings!
            }
            sumE+=par->DeltaE(COMPLETE_LAMBDA,&collect); // Finish Delta Lambda
            accepted = AcceptMove(sumE) and allowed;
            //accepted= false;
            if(!accepted)
            {
               // Restore 
               nt1->SetProperty(l1uname,l1upop);
               nt1->SetProperty(l2uname,old1lu);
               nt2->SetProperty(l1uname,old2lu);
               nt2->SetProperty(l2uname,l2upop);
               nt1->SetProperty(l1dname,l1dpop);
               nt1->SetProperty(l2dname,old1ld);
               nt2->SetProperty(l1dname,old2ld);
               nt2->SetProperty(l2dname,l2dpop);

              modprop[i]->update();
              first=false;counter++;
            }
            else {        
                             par->CacheWrite(cache);
                             #pragma omp critical
                          {
         if(first)
                          naccepted++;}
                          
                 }
                
            } }}  // Can Be Vectorized!}

            
        }  
        return naccepted / (double) modprop.size();
    }
 
protected:
    // Chargepen
    double Charge_Penalty(NTriangle *n)
    {
   //  if(llim == NO_LIMIT and hlim==NO_LIMIT) return 0;
     double penalty=n->GetProperty("Charge");
     penalty*=penalty;
    return penalty;    
    };
    double Entropy_Penalty(NTriangle *n)
    {
   //  if(llim == NO_LIMIT and hlim==NO_LIMIT) return 0;
     double penalty=n->GetProperty("Entropy");
     penalty*=-1./Beta;
    return penalty;    
    };
};



 
