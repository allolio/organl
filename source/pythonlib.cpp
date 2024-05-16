#include <time.h>
#include <unordered_set>
#include <map>
#include <deque>
#include "config_ng.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;

#include "mesh.hpp"
#include "property.hpp"
#include "collide_nb.hpp"
#include "face_eval.hpp"
#include "geocontrol.hpp"
#include "energy.hpp"
#include "nagata.hpp"


bool Nagata3Face(vector<Vertex*> &v, vector<Edge*> &e , NTriangle &nt)
{
    /* Transparency :)
     */
    /** Make Normals Consistent with Nagata Triangle intrinsic mesh orientation**/
    triple a=v[1]->pos - v[0]->pos;
    triple b=v[2]->pos - v[1]->pos;
    triple nc= a % b;
    int flip=1;
    for(int i=0;i<3;i++)
    {
        nt.v[i]=&v[i]->pos;
        // Ensure Consistency with Nagata Parametrization, Test orientation!
        if(flip==1 and nc*v[i]->norm < 0) flip=-1;
        if(flip==-1) v[i]->norm=v[i]->norm*-1;
        if(flip==-1 and nc*v[i]->norm <0) {cerr << "INCONSISTENT MESH ORIENTATION - ORIENT MESH! "<< i <<endl;}
        v[i]->norm=v[i]->norm/v[i]->norm.abs(); // Make sure normal vectors are unit!
        nt.n[i]=&v[i]->norm;
      //  nt.c[i]=&e[i]->c;
    }
    return nt.updatec();
}

triple Rotate(triple &vec, double angle, triple &x)
{
  triple a;
  a.x=1;a.y=1;a.z=1;
  double center=(1-cos(angle));
  Matrix3 rot;
  rot.row[0]=a;
  rot.row[1]=a;
  rot.row[2]=a;
  rot=rot*(double) center;
  double sa=sin(angle);
  double ca=cos(angle);
  // Build Rotation Matrix
  rot.row[0]=rot.row[0]*vec.x;
  rot.row[1]=rot.row[1]*vec.y;
  rot.row[2]=rot.row[2]*vec.z;
//
 //  cout << " c " << center << " a " << angle << " t" << rot.row[0].x<< endl; 
  
  rot.row[0].x=vec.x*rot.row[0].x+ca;
  rot.row[0].y=vec.y*rot.row[0].y-vec.z*sa; 
  rot.row[0].z=vec.z*rot.row[0].z+vec.y*sa;
  
 // cout << rot.row[0].x << " " << rot.row[0].y << " " << rot.row[0].z << "\t    | " << x.x  << endl;
  
  rot.row[1].x=vec.x*rot.row[1].x+vec.z*sa;
  rot.row[1].y=vec.y*rot.row[1].y+ca; 
  rot.row[1].z=vec.z*rot.row[1].z-vec.x*sa;

 // cout << rot.row[1].x << " " << rot.row[1].y << " " << rot.row[1].z <<  "\t    | " << x.y  << endl;

  rot.row[2].x=vec.x*rot.row[2].x-vec.y*sa;
  rot.row[2].y=vec.y*rot.row[2].y+vec.x*sa; 
  rot.row[2].z=vec.z*vec.z*(1-ca)+ca; 
 
 // cout << rot.row[2].x << " " << rot.row[2].y << " " << rot.row[2].z <<  " \t   | " << x.z  << endl;

 // cout << endl;
  
  triple res=rot*x;
//  cout << res.x <<"  " << res.y << "  " << res.z << " NORM " << res.abs() << endl; 
  return rot*x;
}


bool ApplyMove(int index, double xax, double yax, double zax, double angle, int normvec, FaceVertexMesh *fvm)
{
if(index > 5) return false;
triple ax; ax.x=xax,ax.y=yax,ax.z=zax;
if(normvec==0 or normvec==2) fvm->vertices[index].pos=    Rotate(ax,angle ,fvm->vertices[index].pos);
if(normvec==1 or normvec==2) fvm->vertices[index].norm=    Rotate(ax,angle ,fvm->vertices[index].norm);
//    fvm->vertices[3].pos.print();
//fvm.vertices[3].norm.print();
return true;
}


double EvalNerror( int aline, int bline, NTriangle &nt1, NTriangle &nt2, int nticks=30,bool back=true)
{
    double xi1inc=0,eta1inc=0,xi2inc=0,eta2inc=0;
    if(aline == 0 or aline==2 ) xi1inc=1./nticks;
    if(bline == 0 or bline==2) xi2inc=1./nticks;
    if(aline == 1 or aline==2) eta1inc=1./nticks;
    if(bline == 1 or bline==2) eta2inc=1./nticks;
    double sumtheta=0;
    cout << " DOOM" << endl;
    nt1.n[0]->print();
    nt2.n[1]->print();
    nt1.n[1]->print();
    nt2.n[0]->print();
    nt1.normal(0,0).print();
    if(*nt1.n[0]*nt1.normal(0,0) < 0.95 ) return -2;
    if(*nt1.n[1]*nt1.normal(0,1) < 0.95) return -2;
    if(*nt2.n[0]*nt2.normal(0,0) < 0.95) return -2;
    if(*nt2.n[1]*nt2.normal(0,1) < 0.95) return -2;

    nt2.normal(1,0).print();
    nt1.normal(1,0).print();
    nt2.normal(0,0).print();
    cout << "_--------------" << endl;
    for(int i=0;i!=nticks;i++)
    {
    nt1.eval(i*eta1inc,i*xi1inc).print() ;
    nt2.eval(1-i*eta2inc,i*xi2inc).print() ;
    sumtheta+= (nt1.normal(i*eta1inc,i*xi1inc) * nt2.normal(1-i*eta2inc,i*xi2inc));
    }
    cout << sumtheta/nticks << "  " << sumtheta/nticks/3.141592*180.0 << endl;
    return sumtheta/nticks; // /M_PI*180.0;
}


double EdgeError(FaceVertexMesh &fvm)
{
    GNagataTriangle g1, g2;
    for(int i=0;i!=2;i++)
    {
        vector<Vertex*> vlist(3);
        vector<Edge*> elist(3);
//      indexlist exposed(3);
        for(int j=0;j<3;j++)
        {
            vlist[j]=&fvm.vertices[fvm.faces[i].vert[j]];
            elist[j]=&fvm.edges[fvm.faces[i].edges[j]];
            //exposed[j]=fvm.faces[i].vert[j];
        }
     if(i==0) Nagata3Face(vlist,elist,g1);
     else Nagata3Face(vlist,elist,g2);
    //   tri[i].exposedv=exposed; tri[i].exposedn=exposed;tri[i].faceID=i;
    }
    return EvalNerror(1,1,g1,g2,10);
}



PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MODULE(pynagata, m)
{
  m.doc() = "Python bindings for the Nagata library";
  py::bind_vector<std::vector<double>>(m, "DoubleVector", py::module_local(false));
  py::class_<FaceVertexMesh>(m, "FaceVertexMesh")
    .def(py::init())
    .def( "LoadObj",&FaceVertexMesh::LoadObj )
    .def( "SaveObj",&FaceVertexMesh::SaveObj );

/* Binding code */
  m.def("ApplyMove",&ApplyMove);
  m.def("EdgeError", &EdgeError);
}
