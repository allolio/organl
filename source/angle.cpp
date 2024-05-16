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
#include "remesh.hpp"

//#include "mc.hpp"

void PrintMeshAngles(FaceVertexMesh *fvm)
{
    for (int i=0;i!=fvm->faces.size();i++)
    {
            triple dv[3];
              for (int k = 0; k < 3; k ++)
        {
            int j = (k + 1) % 3;
            triple a=fvm->vertices[fvm->faces[i].vert[k]].pos;
            triple b=fvm->vertices[fvm->faces[i].vert[j]].pos;
            dv[k]=a-b;
            dv[k].print();
            dv[k]=dv[k]/dv[k].abs();
        }
        cout << i << " An0 " << acos(dv[0]*dv[1]) << endl;
        cout << i << " An1 " << acos(dv[1]*dv[2]) << endl;
        cout << i << " An2 " << acos(dv[2]*dv[0]) << endl;
    }
}

void PrintMeshAnglesEq(FaceVertexMesh *fvm)
{
    for (int i=0;i!=fvm->faces.size();i++)
    {
            double dv[3];
              for (int k = 0; k < 3; k ++)
        {
          
            dv[k]=fvm->vertices[fvm->faces[i].vert[k]].neighbors.size();
            cout << dv[k] << endl;
       //     dv[k].print();
//            dv[k]=dv[k]/dv[k].abs();
        }
        cout << i << " An0 " << (1-2/dv[0])*M_PI << endl;
        cout << i << " An1 " << (1-2/dv[2])*M_PI << endl;
        cout << i << " An2 " << (1-2/dv[1])*M_PI << endl;
    }
}

int main(int argc, char **argv) {
    FaceVertexMesh fvm;
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
    }
    if(!fvm.LoadObj(filename)) cout << "FILE COULD NOT BE LOADED" << endl;
            if (argc < 3)
        PrintMeshAngles(&fvm);
        else PrintMeshAnglesEq(&fvm);
    return 0;
}
