#pragma once
#define INVALID_RATE -61234
#define MC_INVALID -1
class SimulationProperty
{
    public:
        // constructors
        SimulationProperty(){}
        SimulationProperty(string pname)
        {
            name = pname;
        }

        // methods
        inline void setProperty(double prop)
        {
            value = prop;
        }
        
          inline void setSProperty(string prop)
        {
            svalue = prop;
        }
        inline double getProperty(void)
        {
            return value;
        }
        
        inline string getSProperty(void)
        {
            return svalue;
        }
        string getName(void)
        {
            return name;
        }

    private:
        string svalue;
        string name;
        double value;
};

class SimulationSetup
{
    public:
        SimulationSetup(){chf=NULL;Edef.clear();Elam.clear();DefRate.clear();}
        bool ReadProperties(string file);

        bool AddProperty(SimulationProperty &sp)
        {
            props.insert(make_pair(sp.getName(), sp));
            return true;
        }

        int getNumberOfMonteCarloSteps()
        {
            if (props.find("nMCSteps") != props.end())
                return (int)(props["nMCSteps"]).getProperty();
            return 20000;
        }

        int getOutputFrequency()
        {
            if (props.find("writeEvery") != props.end())
                return (int)(props["writeEvery"]).getProperty();
            return 1000;
        }
        
        string getReportFile()
        {
             if (props.find("RFile") != props.end())
                return props["RFile"].getSProperty();
            return string("report.log");
        }
        
        int getReportFreq()
        {
              if (props.find("RFreq") != props.end())
                return (int)(props["RFreq"]).getProperty();
            return -1;
            
        }
        
          int getRemeshingSession()
        {
            if (props.find("remeshEvery") != props.end())
                return (int)(props["remeshEvery"]).getProperty();
            return -1;
        }
        
           int getRemeshingIterations()
        {
            if (props.find("remeshIter") != props.end())
                return (int)(props["remeshIter"]).getProperty();
            return 100;
        }
          int getRateFrequency()
        {
            if (props.find("updateDefEvery") != props.end())
                return (int)(props["updateDefEvery"]).getProperty();
            return -1;
        }
        
        double getRate(int index)
        {
            if(index>=DefRate.size() ) return 0;
            else return DefRate[index];
        }

	bool DeepVertexOnly()
	{
	  if(props.find("DeepVertexOnly") != props.end())	return true;
	  else return false;
	}
	
	bool DomainTension()
	{
	  if(props.find("DomainTension") != props.end())	return true;
	  else return false;
	}
	
        
        bool SetDef(int index, double def)
        {
            if(index >= chf->defaults.size()) {cerr << "DEF: Value is out of range for energy: " << index << " Maximum number is  " << chf->defaults.size()-1 << endl; exit(-1);}
            chf->defaults[index]=def;
            return true;
        }

        int getRefreshNeighborsFrequency()
        {
            if (props.find("updateNghbrsEvery") != props.end())
                return (int)(props["updateNghbrsEvery"]).getProperty();
            return 200;
        }
        
         int getAutoTuneLimit()
        {
            if (props.find("autoTuneUntil") != props.end())
                return (int)(props["autoTuneUntil"]).getProperty();
            return 0;
        }

         double getBeta()
        {
            if (props.find("Beta") != props.end())
                return (double)(props["Beta"]).getProperty();
            return 1.0;
        }

        double getNeighborhoodRadius()
        {
            if (props.find("nghbRadius") != props.end())
                return (props["nghbRadius"]).getProperty();
            return 1e18;
        }
                
        bool SetEnergyParams(Energy *e)
        {
              cout << "setting params "  << endl;
            if(Elam.size() > e->lambdas.size()) 
            {   if (e->lambdas.size()>0)  cerr << "WARNING LAMBDA Index out of Energy range, values higher than "<< e->lambdas.size()-1  << " will be ignored" << endl;
                else cerr << "WARNING LAMBDAS will be ignored" << endl;
            }
            if(Elam.size() < e->lambdas.size()) cerr << "WARNING NOT Enough LAMBDA Values provided: "<< e->lambdas.size() << " required" << endl;
            if(Edef.size() < e->defaults.size()) { cerr << "ERROR NOT Enough DEF Values provided: "<< e->defaults.size() << " required" << endl; exit(-1);
            }
            for(int i=0;i!=e->lambdas.size();i++)
            {if(Elam[i]!=LAMBDA_DEFAULT) e->lambdas[i]=Elam[i];if(i<e->defaults.size()) {if(Edef[i]!=INVALID_BULK) e->defaults[i]=Edef[i];}
             else if(i>0 and i>=e->defaults.size()-1)  {cerr << "DEF: Value is out of range for energy: " << i << " Maximum number is  " << e->defaults.size()-1 << endl; exit(-1);}
            } 
            chf=e;
            cout << "params set " << endl;
            return true;
        }
        
        vector<double> Edef;
        vector<double> Step;
    private:

        unordered_map<string, SimulationProperty> props;
        Energy *chf;
        vector<double> Elam;
        vector<double> DefRate;
      
};


bool SimulationSetup::ReadProperties(std::string file)
{
    // Format: * name value
    // lines that not begin with * are ignored
    ifstream hfile; Edef.clear(); Elam.clear(); DefRate.clear(); Step.clear();
    hfile.open(file.c_str());
    if(hfile.fail()) return false;
    string current;
    while ( !hfile.eof() )
    {
        vector<string> line;
        getline(hfile, current);
        int htokens=Tokenize(current,&line, " \t\r");
        if(htokens<3)
            continue;
        //cout << line[1] << " "<< line[2] << endl;

        if(line[0]=="S")
        {
            SimulationProperty sprop(line[1]);
	    if (line.size()>=3) {
            if (isNumber(line[2])) sprop.setProperty(atof(line[2].c_str()));
            else sprop.setSProperty(line[2]);}
            AddProperty(sprop);
        }
        
        if(line[0]=="E")
        {
            if (line[1]=="DEF")
            {
                int index=atoi(line[2].c_str());
                if(Edef.size()<index+1) Edef.resize(index+1,INVALID_BULK);
                Edef[index]=atof(line[3].c_str());
         
            }
            if (line[1]=="INCR")
            {
                int index=atoi(line[2].c_str());
                if(DefRate.size()<index+1) DefRate.resize(index+1,INVALID_RATE);
                DefRate[index]=atof(line[3].c_str());
            }

            if (line[1]=="LAMBDA")
            {
                int index=atoi(line[2].c_str());
                if(Elam.size()<index+1) Elam.resize(index+1,LAMBDA_DEFAULT);
                Elam[index]=atof(line[3].c_str());
            }
        }
       if(line[0]=="M")
           if (line[1]=="STEP")
           {
               int index=atoi(line[2].c_str());
               if(Step.size()<index+1) Step.resize(index+1,MC_INVALID);
                Step[index]=atof(line[3].c_str());
           }
    }

    hfile.close();
    return true;
}

struct Simulation;
struct LogReport;
#include"mc_head.h"


struct LogReport
{
    bool Start(string fname, int tfreq,Simulation *tsim, bool textra=false)
    {
     freq=tfreq;
     if(tfreq <0 ) return false;
     fhandle.open(fname, ofstream::out | ofstream::app);
     sim=tsim; counter=1; extra=textra;
     return true;
    }
   bool Tick() ;
    bool Stop()
    {
     fhandle.close();
     return true;
    }
    ~LogReport()
    {
     fhandle.close();   
    }
protected:
   bool extra;
   ofstream fhandle;
   Simulation *sim;
   int freq;
   unsigned long int counter;
};

struct Simulation
{
  Simulation() {cd.n=&n; interaction=NULL; Remesher=0;}
  bool AddEntity(Entity *ent); // Can be refactured further most easily. Relevant interfaces are Propertied, CollidableCollection, Some Energy Interface.
  vector<Entity*> Entities;
  vector<MCMove*> Moves;
  MCMove* Remesher;
  Energy *interaction; //
  Neighborlist n;
  CollisionDetection cd;
  bool SetUp(SimulationSetup &);
  bool Init(void);
  int Run(void);
  vector<double> * GetDefaults();
protected:
  double LocalEnergy(int node, int depth);
  SimulationSetup thesim;
  LogReport logger;

};



  bool LogReport::Tick()
    {
            
            if(freq <0 ) return false;
            if(counter % freq == 0 )
            {
            if (counter == 1) fhandle << "# Step \t Entity \t Energy \t Area \t Volume \t Integrated Mean Curv. \t  V0 \t =Reduced Volume \t Constraint Energy \t [Labeled Extra Terms]" << endl; 
            //cout << "DID" << i  << " MOVES  " << Entities[j]->TotalEnergy(1) << endl;
            //cout << "REAL Energy " << Entities[j]-> TotalEnergy(0) << endl;
                fhandle << counter <<  " \t ";
                for(int j=0; j!=sim->Entities.size();j++)
                {
                fhandle << j <<  " \t ";
                fhandle <<  sim->Entities[j]->TotalEnergy(1) <<  " \t "; 
                if(sim->Entities[j]->props.HasProp("MembG")) 
                {
                double v[3];
                for(int k=0; k!=3;k++) 
                 { v[k]= (double) *(sim->Entities[j]->props.GetProperty("MembG",k,false).content) ;
                   fhandle << v[k] << " \t "  ; 
                 }
                if(extra)
                    {
                     fhandle << " V0 \t "  << 6*v[1]*pow(v[0],-3./2.)*sqrt(M_PI) << " \t " ;
                    }
                }
                fhandle << sim-> Entities[j]->TotalEnergy(2) << " \t ";
                if(extra)
                    {

                        fhandle << "G1 Diagnostic: " << G1Diagnostic(sim->Entities[j]) << " \t " ;
                    }
                      if(extra)
                    {

                        fhandle << "Edge Terms: " << sim->Entities[j]->TotalEnergy(5) << " \t " ;
                    }
                }
            fhandle << " D \t";
            vector<double > *Edef=sim->GetDefaults();
            for(int l=0;l!=Edef->size();l++) fhandle << l << " " << Edef->at(l) << " \t " ;
            fhandle << endl;
            }
            
            counter++;
            return true;
    }

bool Simulation::Init()
{
 for (int i=0;i!=Moves.size();i++)
 {
     Moves[i]->Beta=thesim.getBeta();
 }
 for (int i=0;i!=Entities.size();i++) Entities[i]->ReadyToRun();
//  for (int i=0;i!=Entities.size();i++) Entities[i]->ReadyToRun();
 n.createNearestNeighborList(thesim.getNeighborhoodRadius());
 logger.Start(thesim.getReportFile(),thesim.getReportFreq(),this, true);
  thesim.Step.resize(Moves.size(),MC_INVALID);
 for (int i=0;i!=thesim.Step.size();i++) 
 {     if( thesim.Step[i]!=MC_INVALID)
         Moves[i]->maxstepsz=thesim.Step[i];
         cout << "MOVE " << i << " STEP " << Moves[i]->maxstepsz << endl;
 }
     
 return true;
}
/**
 * This is the main simulation loop
 **/
int Simulation::Run()
{
    for(int i=1;i<=thesim.getNumberOfMonteCarloSteps();i++)
    {
        if(i%thesim.getRefreshNeighborsFrequency()==0)  { n.createNearestNeighborList(thesim.getNeighborhoodRadius()); cout << "COLL" << thesim.getNeighborhoodRadius() << endl;}
        if(thesim.getRateFrequency()!=-1)
        if(i%thesim.getRateFrequency()==0)
        {
            cout << "Def Update" << endl;
            int k=0;
            double l=thesim.getRate(k);
            while(l!=0)
            {
                if(l!=INVALID_RATE)
                thesim.Edef[k]+=l;
                k++;
                l=thesim.getRate(k);
            }
            if(k!=0)
            {
             for(int l=0;l!=k;l++) thesim.SetDef(l,thesim.Edef[l]);   
            }
        }
        if(thesim.getRemeshingSession()!=-1)
        if(i%thesim.getRemeshingSession()==0)
        {
            cout << " REMESHING SESSION" << endl;
            for(int rmsh=0; rmsh!=thesim.getRemeshingIterations();rmsh++)
            {
               int cycle=1;
               if(Remesher==NULL) cout << "Cannot Remesh - Set Remesher First!" << endl;
               
                       else    while(cycle>0) {
                cycle=Remesher->CreateParallelPlan(this);
                Remesher->ExecuteMoves(this);
                           }
                cout << "REMESH CYCLE " << rmsh << endl;
                for(int k=0;k!=Entities.size();k++)
                Entities[k]->RestoreDependencies();
                     for(int j=0;j<Moves.size();j++)
                     {
           double ratio;
//         int cprev=0;
           cycle=1;
           Moves[j]->Reset(this);
        
           while(cycle>0) {
            cycle=(Moves[j])->CreateParallelPlan(this);//->BuildParralelPlan();
       //  cout << "cycle" << cycle << "j " << endl;
           double ratio=Moves[j]->ExecuteMoves(this);
            cout<< "MOVE " << j << " CYCLE " << cycle << " STEPSIZE NOW " <<  Moves[j]->maxstepsz  << endl;  

           }

            
            }
            
          }
          n.createNearestNeighborList(thesim.getNeighborhoodRadius()); cout << "COLL" << thesim.getNeighborhoodRadius() << endl;
                      cout << " REMESHING DONE" << endl;
        }
        logger.Tick();
        if(i%thesim.getOutputFrequency()==0)
        {
            for(int j=0;j!=Entities.size();j++)
            {
            string str;
            stringstream ss, rs;
            ss << "prog" << j << "_" << i << ".vtu";
            ss >> str;
            rs << "restart" << j << "_" << i << ".obj";
            Entities[j]->WriteVtu(str,true);
            rs >> str;
            Entities[j]->SaveToFile(str); // FIXME: Needs to be generalized, completed obviously.
            cout << "DID" << i  << " MOVES  " << Entities[j]->TotalEnergy(1) << endl;
            cout << "REAL Energy " << Entities[j]-> TotalEnergy(0) << endl;
            Entities[j]->SumStored("Area");
            Entities[j]->props.GetProperty("MembG",0,true);
            Entities[j]->SumStored("Volume");
            Entities[j]->props.GetProperty("MembG",1,true);
            Entities[j]->SumStored("Curvature");
            Entities[j]->props.GetProperty("MembG",2,true);
            cout << "Constraint" << endl;
            cout << Entities[j]->TotalEnergy(2) << endl;
            cout << "Boundary and Extra Terms" << endl;
            cout << Entities[j]->TotalEnergy(5) << endl;
            cout << "Defaults" << endl;
            for(int l=0;l!=thesim.Edef.size();l++) cout << l << " " << thesim.Edef[l] << endl;
              cout << "G! Diagnostic: " << G1Diagnostic(Entities[j]) << endl;
            }
        }

     for(int j=0;j<Moves.size();j++)
     {
         double ratio;
//         int cprev=0;
         int cycle=1;
         while(cycle>0) {
         cycle=(Moves[j])->CreateParallelPlan(this);//->BuildParralelPlan();
  //       cout << "cycle" << cycle << "j " << endl;
         if(cycle==1) ratio=(Moves[j])->ExecuteMoves(this);
         else Moves[j]->ExecuteMoves(this);
        }
         if( i<thesim.getAutoTuneLimit() )
         {
		 if(Moves[j]->prop!="Normal")
		 {
             
         if(ratio < 0.25) Moves[j]->maxstepsz*=0.95;
         if(ratio > 0.4) Moves[j]->maxstepsz*=1.05;}
         cout<< "MOVE " << j << " ACCEPTANCE " << ratio << " STEPSIZE NOW " <<  Moves[j]->maxstepsz  << endl;  
         }
     }
     
    }
    return 0;
}
vector<double>* Simulation::GetDefaults()
{
    
    return &thesim.Edef;
}

double Simulation::LocalEnergy(int node, int depth=1)
{
    if(depth==1)
    {
        indexlist &neigh=n.neighbors[node];
        for (int i=0;i!=neigh.size();i++)
        {
            cout << neigh[i] << "D:  " << n.dsq[node].at(i) << "E " << Entities[n.parent_entity[ neigh[i] ] ]->IndexEnergy(n.parent_index[neigh[i]]) << "PI"<< n.parent_index[neigh[i]] <<endl;
            // TODO ADD INTERACTION THROUGH SPACE.
        }
        return 0;
    }
    else cerr << "High Depth not implemented" << endl;
    return -1;
}

bool Simulation::SetUp(SimulationSetup &stp)
{
    thesim=stp;
    // Build the Neighbor List
    
    return true;
}

bool Simulation::AddEntity(Entity *ent)
{
    // Add Entity
    Entities.push_back(ent);
    // Insert Neighbors.
    vector <Neighbor *> nbadd;
    ent->exposeNeighbors(nbadd);
    (n.col).insert((n.col).end(),nbadd.begin(), nbadd.end());
    indexlist parents(nbadd.size());
    indexlist parents_ent(nbadd.size());
    for(int i=0;i!=parents.size();i++)
     { parents[i]=i; parents_ent[i]=Entities.size()-1;}
    n.parent_index.insert(n.parent_index.end(),parents.begin(),parents.end());
    n.parent_entity.insert(n.parent_entity.end(),parents_ent.begin(),parents_ent.end());
    return true;
}



