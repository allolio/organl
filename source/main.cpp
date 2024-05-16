#include <time.h>
#include <unordered_set>
#include <map>
#include <deque>
#include "config_ng.hpp"
#include "mesh.hpp"
#include "property.hpp"
#include "collide_nb.hpp"
#include "face_eval.hpp"
#include "geocontrol.hpp"
#include "energy.hpp"
#include "nagata.hpp"
typedef Membrane Entity; // A stand-in [BF: a stand-in for what?]

#include "blockparser.hpp"
#include "checkqg1.hpp"

#include "simulation.hpp"
#include "monte_carlo.hpp"
#include "combomove.hpp"
#include "deep_move.hpp"
#include "lipid.hpp"
#include "remesh.hpp"
#include <cstring>

int main(int argc, char **argv) {

    /* ----------------------- *
     * Parse input arguments   *
     * ----------------------- */
    cout << "OrganL 1.01 release " << __DATE__ << " " <<  __TIME__ << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Cite: Allolio, Fábián, Dostalík, Biophys. J.  " << endl;
    cout << "https://doi.org/10.1016/j.bpj.2024.04.028" << endl;
    cout << "---------------------------------------------" << endl;
    
    if (argc <= 1) {
        // No arguments: run an example
        cout << "No input file specified. Usage is organl {filename.obj} [ ADE | CHF | EAP | EAH | DCH | DCG ] [ LIP ] [ BNE | TAU ] \n";
        exit(-1);
    }
    string filename = argv[1];
    cout << "From mesh " + filename + "\n";
   

    FaceVertexMesh fvm;
    Energy *dchf;
     LipidDB ldb;
    bool bne=false;
    bool tau=false;
    bool mix=false;
    Membrane *pmmb;
    if(argc>2)
    {
        if(strcmp(argv[2],"ADE")==0)
        {
           dchf=new DCADEHelfrichEnergy;
           cout << "ADE Energy " << endl;
        }
       else if(strcmp(argv[2],"CHF")==0)
        {
           dchf=new ConstraintHelfrichEnergy;
           cout << "Constrained Helfrich Energy " << endl;

        }
        else if(strcmp(argv[2],"EAH")==0)
        {
           dchf=new ElasticAreaHelfrichEnergy ;
           cout << "Elastic Area  Helfrich Energy " << endl;
        }
           else if(strcmp(argv[2],"EAP")==0)
        {
           dchf=new HelfrichElasticAreaPressure ;
           cout << "Pressure + Elastic Area Helfrich Energy " << endl;
        }
           else if(strcmp(argv[2],"DCG")==0)
        {
           dchf=new DoubleConstraintGHelfrichEnergy;    
           cout << " Doubly Constrained Helfrich Energy with Gaussian Curvature" << endl;
           cout << " Bonnet-Boundary terms are OFF by default" << endl;
        }
            else if(strcmp(argv[2],"DCH")==0)
        {
           dchf=new DoubleConstraintHelfrichEnergy ;
             cout << "Constrained Helfrich Energy with Area Elasticity" << endl;
        }
          else if(strcmp(argv[2],"LIP")==0)
        {
            mix=true;
            dchf=new DoubleConstraintHelfrichEnergy;
        }
              else if(strcmp(argv[2],"BNE")==0)
        {
            bne=true;
            dchf=new DoubleConstraintHelfrichEnergy;
        }
                      else if(strcmp(argv[2],"TAU")==0)
        {
            tau=true;
            dchf=new DoubleConstraintHelfrichEnergy;
        }

        else {
            cout << "Cannot understand Energy Functional selection exiting. " << endl;
            exit(-1);
         }
    
    }
       else {
            dchf=new DoubleConstraintHelfrichEnergy;
             cout << "Constrained Helfrich Energy with Area Elasticity" << endl;}
 
    if(argc>3)
    {
          if(strcmp(argv[3],"LIP")==0)
        {
            mix=true;
        }
         if(strcmp(argv[3],"BNE")==0)
        {
            bne=true;
        }
            if(strcmp(argv[3],"TAU")==0)
        {
            tau=true;bne=false;
        }

    }
    
   if(argc>4)
    {
          if(strcmp(argv[4],"LIP")==0)
        {
            mix=true;
        }
         if(strcmp(argv[4],"BNE")==0)
        {
            bne=true;
        }
          if(strcmp(argv[4],"TAU")==0)
        {
            tau=true;bne=false;
        }

    }
    
    
    
    
    if (!fvm.LoadObj(filename)) {
        cout << "ERROR: file " + filename + " could not be loaded!" << endl;
        return 1;
    }
     Energy &chf=*dchf;
  
      if(mix)
    {
      cout << "LIPID MIXING ENABLED" << endl;
      if(!ldb.LoadLipids("lipids.txt"))
      {cerr << "Failed to load lipid database expected in lipids.txt. File not found" << endl; exit(-1);}
      cout << "LIPIDS INCLUDED" << endl;
     ldb.DisplayLipids();
      pmmb= new MixedMembrane(chf);
    }
    else
    {
      pmmb= new Membrane(chf);
    }
    
    /* ----------------------- *
     * Set up stuff... *
     * ----------------------- */
    Membrane &mmb=*pmmb;
    mmb.FromMesh(fvm);
    if(mix)
              ((MixedMembrane*) pmmb)->SetLipids(&ldb);

  //  PrintMeshAngles(&fvm);

   
    double avglength=0;
    for(int i=0; i!=fvm.edges.size();i++)
    {
        avglength+=(fvm.vertices[fvm.edges[i].vertA].pos-fvm.vertices[fvm.edges[i].vertB].pos).abs();
    }
    avglength/=fvm.edges.size();
    cout << "INFO: Average Edge Length" << avglength<< endl;
 
    /*
     *
     */
    cout << "READING FACE PROPERTIES" << endl;
    if(!mmb.ReadFaceProperties("props.txt")) { cout << "ERROR: FAILED TO READ FACE PROPERTIES from props.txt" << endl; exit(-1);}
    cout << "------------- END FACE PROPERTIES" << endl;
    cout << " READING BLOCKS from blocks.txt" << endl;
    if(!ParseBlocks("blocks.txt",&mmb)) cout << "No blocks.txt file read. Leaving geometry free" << endl ;
    cout << "------------- END BLOCKS" << endl;
    
    if(bne)
    {
      cout << "ADDING GAUSS-BONNET GEODESIC EDGE TERMS" << endl;
      BonnetEdgeFactory bnef;
      pmmb->setBoundaryEnergy(bnef); 
    }
     if(tau)
    {
      cout << "ADDING LINE-TENSION EDGE TERMS" << endl;
      TauEdgeFactory tnef;
      pmmb->setBoundaryEnergy(tnef); 
    }

  

    Simulation sim;
    SimulationSetup setup;
    cout << "READING CONTROL PARAMS" << endl;
    if(!setup.ReadProperties("control.txt")) { cout << "ERROR: FAILED TO READ CONTROL PARAMETERS from control.txt: File not found" << endl; exit(-1);};
    bool dvmonlyflag=setup.DeepVertexOnly();
    setup.SetEnergyParams(&chf);
    if(setup.DomainTension())
    {
            if(bne)
    {
      cout << "ADDING GAUSS-BONNET GEODESIC DOMAIN TERMS" << endl;
      BonnetEdgeFactory bnef;
      pmmb->setBoundaryEnergy(bnef,true); 
    }
     if(tau)
    {
      cout << "ADDING LINE-TENSION DOMAIN TERMS" << endl;
      TauEdgeFactory tnef;
      pmmb->setBoundaryEnergy(tnef,true); 
    }

    }
    sim.SetUp(setup);
    fvm.SaveObj("ttt.obj");
    mmb.WritePlot("out.plt",10);
    sim.AddEntity(&mmb);
    ANormalMove anormal;
    VertexMove vertex; 
    DeepVertexMove dvm; 
    AlexanderMove rmsh;
    LipidMixMove lmm;
    
    // Defaults
    lmm.maxstepsz=0.1;
    dvm.maxstepsz=avglength;
    vertex.maxstepsz=avglength*1.2;
    if(dvmonlyflag)   {
	    cout << "Using DeepVertexMove only" <<endl;
	    sim.Moves.push_back(&dvm);
    }
    else {
	    sim.Moves.push_back(&vertex);
   	    sim.Moves.push_back(&anormal);
         }
    if(mix) sim.Moves.push_back(&lmm);
    sim.Remesher=&rmsh;
    sim.Init();
    if(chf.defaults.size()>1 and chf.lambdas.size()>1 ) {
        cout << " A0 " << chf.defaults[0] << " LAMBDA A " << chf.lambdas[0]<<endl;
        cout << " V0 " << chf.defaults[1] << " LAMBDA V " << chf.lambdas[1]<<endl;
    }
    else if(chf.lambdas.size()>0) {   cout << "LAMBDA 0:" << chf.lambdas[0] << endl;
    if(chf.lambdas.size()>1)    cout << "LAMBDA 1:" << chf.lambdas[1] << endl;}
    cout << "Constraint Energy" << endl;
    mmb.TotalEnergy(2);
    cout << "Initial Energy" << endl;
    mmb.TotalEnergy();
    mmb.TotalEnergy(6);

    sim.Run();
    mmb.TotalEnergy();
    mmb.TotalEnergy(1);
    mmb.WriteVtu("out.vtu", true);            // save as curved elements
    mmb.WritePlot("out.plt",10);
    fvm.SaveObj("out.obj");
    mmb.SaveFaceProperties("props_out.txt");

    return 0;
}
