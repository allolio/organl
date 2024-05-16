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
    cout << "Autoboundary " << __DATE__ << " " <<  __TIME__ << endl;
    cout << "---------------------------------------------" << endl;
    
    if (argc <= 1) {
        // No arguments: run an example
        cout << "No input file specified. Usage is auto {filename.obj}  \n";
        exit(-1);
    }
    string filename = argv[1];
    cout << "From mesh " + filename + "\n";
   

    FaceVertexMesh fvm;
     LipidDB ldb;
    bool bne=false;
    bool tau=false;
    bool mix=false;
    ElasticAreaHelfrichEnergy dchf;
    Membrane mmb(dchf);
        
    
    if (!fvm.LoadObj(filename)) {
        cout << "ERROR: file " + filename + " could not be loaded!" << endl;
        return 1;
    }
    
    
    
    /* ----------------------- *
     * Set up stuff... *
     * ----------------------- */
    mmb.FromMesh(fvm);

  //  PrintMeshAngles(&fvm);

    /*
     *
     */
    cout << "READING FACE PROPERTIES" << endl;
    if(!mmb.ReadFaceProperties("props.txt")) { cout << "ERROR: FAILED TO READ FACE PROPERTIES from props.txt" << endl; exit(-1);}
    EdgeProperties *ep; string key("Edge");
    mmb.props.GetPropertyList(key,(PropertyList **) &ep);
    for(int i=0;i!=ep->theProps.size();i++)
    {
     bool second=false; NTriangle *nt1,*nt2;
     for(int j=0;j!=ep->theProps[i]->deps.size();j++)
       {
       if(ep->theProps[i]->deps[j]->ptype->name=="FaceE")
       {
         if(!second) {nt1=(NTriangle*) ep->theProps[i]->deps[j]->content;nt2=NULL;}
         else { nt2=(NTriangle*) ep->theProps[i]->deps[j]->content;break; }
         second=true;
       }
       if(second and nt2!=NULL) break;
       }
       if(second and nt2!=NULL)
       if(nt1->GetProperty("c0") != nt2->GetProperty("c0"))
       {
         Edge *edge= (Edge*) ep->theProps[i]->content;
         cout << "Edge " << edge->vertA << " "<< edge->vertB << endl;
       }
    }
    return 0;
}
