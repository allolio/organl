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
typedef Membrane Entity; // A stand-in 

#include "blockparser.hpp"
#include "checkqg1.hpp"

#include "simulation.hpp"
#include "monte_carlo.hpp"
#include "combomove.hpp"
#include "deep_move.hpp"
#include "lipid.hpp"
//#include "mc.hpp"


int main(int argc, char **argv) {
    FaceVertexMesh fvm;
    
 //   LipidDB ldb;
 //   ldb.LoadLipids("lipids.txt");
 //   ldb.DisplayLipids();
  /*  triple n1,n2,c1,c2,v1,v2;
    v1.x=0.0177578;v1.y= 0.853939; v1.z=0.493162;
    v2.x=-0.354837,v2.y=-0.463535,v2.z= 0.733253;
    n1.x=5.20874e-17;n1.y= -0.850652;n1.z= -0.52573;
    n2.x=0.127228;n2.y=-0.864738, n2.z=-0.485841;
    Thirdorder(&v1,&v2,&n1,&n2,&c1,&c2);
    c1.print();
    c2.print();
    */
    //string filename="testmesh/uvsphere.obj";
    string outfile;
    string filename;
    if (argc <= 1)
    {
        // filename = "opt.obj";
        filename = "testmesh/kube_scaled.obj";
        outfile  = "out.vtu";
    }
    else
    {
        filename = argv[1];
        if (argc < 3)
        {
            size_t pos;
            pos = filename.find_last_of(".");
            if (pos == string::npos)
                pos = filename.size() - 1;
            outfile = filename.substr(0, pos) + ".vtu";
        }
        else
            outfile = argv[2];
    }
    if(!fvm.LoadObj(filename)) cout << "FILE COULD NOT BE LOADED" << endl;
    ConstraintHelfrichEnergy chf;
    Membrane mmb(chf);
    //mmb.SetLipids(&ldb);
    mmb.FromMesh(fvm);
   
    /*
     * Dirty
     */
    double avglength=0;
    for(int i=0; i!=fvm.edges.size();i++)
    {
        avglength+=(fvm.vertices[fvm.edges[i].vertA].pos-fvm.vertices[fvm.edges[i].vertB].pos).abs();
    }
    avglength/=fvm.edges.size();
    cout << "AVG Edge Length" << avglength<< endl;
    double llim=0.8*avglength;
    double hlim=2*avglength; //1.2
    
    /*
     * 
     */
    mmb.ReadFaceProperties("props.txt");
 //  chf.lambdas[0]=1500;
 //  chf.lambdas[1]=400;
      mmb.ReadyToRun();

    if(chf.defaults.size()>0)
    {
        cout << " A0 " << chf.defaults[0] << " LAMBDA A " << chf.lambdas[0]<<endl;
    cout << " V0 " << chf.defaults[1] << " LAMBDA V " << chf.lambdas[1]<<endl;
    if(chf.defaults.size()>2)   cout << " H0 " << chf.defaults[2] << " LAMBDA H " << chf.lambdas[2]<<endl;
 
    }
    cout << "LAMBDA 0:" << chf.lambdas[0] << endl;
    cout << "Constraint Energy" << endl;
          cout << "GOT Here" << endl;
       mmb.TotalEnergy(2);
    cout << "Initial Energy" << endl;
    mmb.TotalEnergy();
    cout << "DONE WITH RUN" << endl;
    mmb.TotalEnergy();
    mmb.TotalEnergy(1);
   // mmb.TotalArea();
   // mmb.TotalVolume();
   // mmb.WriteVtu("Man.vtu");

        mmb.WritePlot("out.plt",10);
    //nag.WriteVtu(outfile, false);       // save straight triangles
    return 0;

}
