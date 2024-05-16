// Nagata Patch Implementation

#pragma once

#define NAGATA_PARALLEL 0.9999
#define NAGATA_EPS1 0.0
#define NAGATA_EPS2 0.26            // \approx cos (80 deg) + cos (85 deg)

#define MEMBP_LAZY -1

#include "triangle.hpp"
//typedef AdvancedNagataTriangle NTriangle;
typedef GNagataTriangle NTriangle;
#include "meshproperties.hpp"
#include "custom_energy.hpp"

class Nagata : public Propertied, public NeighborCollection
{
        friend class FaceProperties;
public:
    Nagata() {/*frozen.clear();gradready=false;gradfirst=true;*/}
    int FromMesh(FaceVertexMesh &fvm);
    void WritePlot(string filename, int tics=5);
    bool WriteVtu(const string &filename, bool quad);
    virtual bool exposeNeighbors(vector<Neighbor*> &col);
   
protected:
    int Nagata3Face( vector<Vertex*> &v, vector<Edge*> &e, NTriangle &nt);
    FaceVertexMesh *fv;
    vector<NTriangle> tri;
private:
    GeoControl g;
    bool AddMidpoint(int f, int v1, int v2);
    bool BuildMidpoints();
};
 
class MembraneProperties: public PropertyList
{
public:
    GeoControl g;
    MembraneProperties()
    {name="MembG";
     szProperty=1;   
    }
    vector<double> lambdas;
    virtual bool genpropertiesfromStructure(void *parent);
};

class MembraneProperty;

class Membrane : public Nagata
{
    friend class MembraneProperties;
public:
    Membrane(Energy &en): MembE(en)
    {
//        MembE=new ConstraintHelfrichEnergy;
    };
    bool AddBlock(Block &b,string &key, int index);
    int FromMesh(FaceVertexMesh &fvm);
    int ReadFaceProperties(string file);
    int SaveFaceProperties(string file);
    int FromFile(string file);
    int SaveToFile(string file);
    bool setEnergy(Energy &e);
    double TotalEnergy(int);
    double SumStored(string key);
    double IndexEnergy(int i);
    double DeltaE(int i, vector<double> *coll);
    double DExtraE(int i);
    bool setBoundaryEnergy(CustomDepEnergyFactory &cde, bool forfreeze);
    bool ReadyToRun();
    bool CacheWrite(indexlist &cache);
    virtual bool RestoreDependencies();
protected:
    virtual bool PopulateProps();
    bool ChainBoundary();
    NormalProperties normp;
    VertexProperties vertp;
    MembraneProperties membp;
    FaceEnerProperties fep;
    EdgeEnerProperties eep; //These are not by default used/connected to mesh in.
    EdgeProperties edgep;
    Energy &MembE;
    vector<string> deepProps;
    vector<int> deepPropIndex;
    vector<CustomDepEnergy*> extraterms;
};


class AreaProperty: public Property
{
public:
      virtual bool update(bool deep=true)
      {
      vector<NTriangle>* tri=(vector<NTriangle>*) parent;
      double Artot=0;
      #pragma omp parallel for reduction(+: Artot)
      for(int i=0;i<tri->size();i++) Artot+=ptype->g.Area(&tri->at(i));
      *content=Artot;
      state=0;
      cout << "AREAUPDATE" << endl;
         return true;
      }
    MembraneProperties *ptype;

};

class VolumeProperty: public Property
{
public:
      virtual bool update(bool deep=true) {
      vector<NTriangle>* tri=(vector<NTriangle>*) parent;
      double Vtot=0;
      #pragma omp parallel for reduction(+: Vtot)
      for(int i=0;i<tri->size();i++) Vtot+=ptype->g.VIncr(&tri->at(i));
      *content=Vtot;
      state=0;
      cout << "VOLUMEUPDATE" << endl;
         return true;
        }
    MembraneProperties *ptype;

};

class CurvatureProperty: public Property
{
public:
      virtual bool update(bool deep=true) {
      vector<NTriangle>* tri=(vector<NTriangle>*) parent;
      double Htot=0;
      #pragma omp parallel for reduction(+: Htot)
      for(int i=0;i<tri->size();i++) Htot+=ptype->g.HInt(&tri->at(i));
      *content=Htot;
      state=0;
      cout << "CURVATUREUPDATE" << endl;
         return true;
        }
    MembraneProperties *ptype;

};


bool MembraneProperties::genpropertiesfromStructure(void* parent)
{
        // 0 = Area // 1= Volume
        Membrane *fv=(Membrane*) parent;
        theProps.resize(3);
        lambdas.resize(3);
            lambdas[0]=MEMBP_LAZY;
            lambdas[1]=MEMBP_LAZY;
            lambdas[2]=MEMBP_LAZY;
            theProps[0]=new AreaProperty;
            theProps[0]->block=nullptr;
            theProps[0]->parent=&fv->tri;
            theProps[0]->content=&lambdas[0];
            //theProps[i].deps=NULL;
            theProps[0]->state=MEMBP_LAZY;
            theProps[0]->ptype=this;
            theProps[0]->index=0;            
            theProps[1]=new VolumeProperty;
            theProps[1]->block=nullptr;
            theProps[1]->parent=&fv->tri;
            theProps[1]->content=&lambdas[1];
            //theProp[i].deps=NULL;
            theProps[1]->state=MEMBP_LAZY;
            theProps[1]->ptype=this;
            theProps[1]->index=1;
            theProps[2]=new CurvatureProperty;
            theProps[2]->block=nullptr;
            theProps[2]->parent=&fv->tri;
            theProps[2]->content=&lambdas[2];
            //theProp[i].deps=NULL;
            theProps[2]->state=MEMBP_LAZY;
            theProps[2]->ptype=this;
            theProps[2]->index=2;      
        return true;
} 

int Membrane::FromFile(string file)
{
// TODO Centralized Read-in of All Parameters
  return 0;  
}

int Membrane::SaveToFile(std::string file)
{
// TODO Save all components;
 fv->SaveObj(file);
 return 0;   
}


bool Membrane::setEnergy(Energy& e)
{
    MembE=e;
    e.deep_props(&deepProps,&deepPropIndex);
    
    return true;
}

bool Membrane::ReadyToRun()
{
    MembE.deep_props(&deepProps,&deepPropIndex); // Just in case we never set an Energy
    double Etot=0;
    membp.theProps[0]->update();
    membp.theProps[1]->update();
    membp.theProps[2]->update();
        if(!MembE.init) MembE.initializeValues(props);
        
          #pragma omp parallel for reduction(+: Etot)
    for(int i=0;i<tri.size();i++)
    {
       // tri[i].updatec(); // Weird but apparently necessary!
        Etot+=MembE.deep_eval(&tri[i]);
    }


     /*   for (int i=0;i!=tri.size();i++)
        {
//            for(int j=0;j!=3;j++)
   //         {
          //  if(tri[i].healthyTriangle)
            {
                triple ne=tri[i].normal(0,0);
                cout << ne* (*tri[i].n[0]) << endl;
                
                ne=tri[i].normal(1,0);
                cout << ne* (*tri[i].n[1]) << endl;

                ne=tri[i].normal(1,1);
                cout << ne* (*tri[i].n[2]) << endl;
            }
    //        }
        }
        
        exit(0);
       */ 
        
 //   MembE.initializeValues(props);
    MembE.E=Etot;
    return true;
}


int Membrane::FromMesh(FaceVertexMesh& fvm)
{
    int rv=Nagata::FromMesh(fvm);
    return rv;
}

int Membrane::ReadFaceProperties(std::string file)
{
    // Format facenum, name, value ; * = default for all faces;
    ifstream hfile;
    hfile.open(file.c_str());
    if(hfile.fail()) return false;
    string current;
    while ( !hfile.eof() )
    {
        vector<string> line;
        getline(hfile, current);
        int htokens=Tokenize(current,&line, " \t\r");
        if(htokens<3) { cout << "Cannot Parse" << endl; break;}
        cout << line[1] << " "<< line[2] << endl;

        if(line[0]=="*")
            for(int i=0;i!=tri.size();i++)
            {
                FaceProperty fprop(line[1]);
                fprop.setProperty(atof(line[2].c_str()));
                tri[i].AddProperty(fprop);
            }
        else
        {
            int facid=atoi(line[0].c_str());
            FaceProperty fprop(line[1]);
            fprop.setProperty(atof(line[2].c_str()));
            tri[facid].AddProperty(fprop);
        }
    }

    hfile.close();
    PopulateProps();
    return true;
}

bool Membrane::RestoreDependencies()
{
         cout << "RESTORING PROPERTIES" << endl;
          PropertyList *pl,*pl2,*fep;
          string prop="Vertex";
          this->props.GetPropertyList(prop,&pl);
          prop="Normal";
          this->props.GetPropertyList(prop,&pl2);
          prop="FaceE";
          this->props.GetPropertyList(prop,&fep);
          FaceVertexMesh *fv;
          fv = (FaceVertexMesh*) pl-> parent;
          
     for(int i=0;i!=fv->vertices.size();i++)
   {
    Property *n=pl2->theProps[i];
    Property *v=pl->theProps[i];
    n->deps.clear();
    // For vertices carry over dependencies which are not faces.
    // Because the edge terms should not be dropped and must be locked away from
    // Remeshing
    
    Deplist cache; cache.clear();
    for(int j=0;j!=v->deps.size();j++)
    {
        if(v->deps[j]->ptype->name!="FaceE") cache.push_back(v->deps[j]);
    }
    
    //v->deps.clear();
    v->deps=cache;
    for(int j=0;j!=fv->vl[i].size();j++)
     {
             Property *f=fep->theProps[fv->vl[i].at(j)];
             n->deps.push_back(f);
             v->deps.push_back(f);
             f->parent=this;
     }
   }
   // Restore Custom Energies
   for(int i=1;i<extraterms.size();i++) 
       extraterms[i]->restore();
   
   return true;
}

int Membrane::SaveFaceProperties(std::string file)
{
        // Format facenum, name, value ; * = default for all faces;
    ofstream hfile;
    hfile.open(file.c_str());
    if(hfile.fail()) return false;
    for(int i=0;i!=tri.size();i++)
    {
        vector< pair <string, double > > vds;
        tri[i].ExportProperties(&vds);
        for(int j=0;j!=vds.size();j++)
            hfile << i << " " << vds[j].first << " " << vds[j].second << endl;
            //        tri[i].GetPropertyList
    }
    
    return true;
}


bool Membrane::AddBlock(Block& b, string& key, int index)
{
    Property &p=props.GetProperty(key,index);//->theProps[index]);
    for(int i=0;i!=b.dim.size();i++)
    {
    cout << i << " " << p.content[b.dim[i]]<< endl;
    if(b.dim[i]> p.ptype->szProperty) return false;
    b.origin[i]=p.content[b.dim[i]];
    }
    if(p.block==NULL)
    {
        blocks.push_back(b);
        p.block=&(*--blocks.end());
        return true;
    }
    if(p.block!=NULL)
    {
        cout << "MERGING BLOCKS" << endl;
        p.block->merge(b);
        return true;
    }
    return false;
}


bool Membrane::PopulateProps()
{
    normp.genpropertiesfromStructure(fv);
    vertp.genpropertiesfromStructure(fv);
    cout << "VERTICES" << vertp.theProps.size() << endl;
    membp.genpropertiesfromStructure(this);
    edgep.genpropertiesfromStructure(fv);
    fep.genpropertiesfromStructure(&tri);
    props.AddPropertyType(normp);
    props.AddPropertyType(vertp);
    props.AddPropertyType(membp);
    props.AddPropertyType(fep);
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
    * Connect Edge Properties to vertices
    */
 
   for(int i=0;i!=fv->vertices.size();i++)
   {
    Property *n=normp.theProps[i];
    Property *v=vertp.theProps[i];
    for(int j=0;j!=fv->vl[i].size();j++)
     {            // VL: Faces for each vertex
         Property *f=fep.theProps[fv->vl[i].at(j)];
         n->deps.push_back(f);
         v->deps.push_back(f);
       /*  Face *fac=&fv->faces[fv->vl[i].at(j)];
         for(int k=0;k!=fac->edges.size();k++)
         {
             if(fv->edges[fac->edges[k]].vertA==i or fv->edges[fac->edges[k]].vertB==i)
             {
            // cout << "Vertex " << i << " Face " << fv->vl[i].at(j) << " Edge " << fac->edges[k] << endl; 
             Property *e=edgep.theProps[fac->edges[k]];
             
             n->deps.push_back(e);
             v->deps.push_back(e);
             e->parent=this;
             }
        }
         */
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
                 Property *f=fep.theProps[fce];
                 n->deps.push_back(f);
                 f->parent=this;
              }
          }
      }
    }
 
   return true;
    
 //   cout << "NORMPROPS " <<  normp.theProps.size() << endl;
 //   cout << "VPROPS " <<  normp.theProps.size() << endl;
}

void Nagata::WritePlot(string filename, int tics)
{
    ofstream ofile;
    ofile.open(filename.c_str());

    int ysplits=tics;
    for(int i=0;i!=tri.size();i++)
    {
        for(int x0=0;x0!=tics;x0++)
        {
            for(int y0=0;y0<=x0;y0++)
            {
                double stepx=1./tics; double stepy=1./ysplits;
                triple o=tri[i].eval(x0*stepx,y0*stepy);
 /*               if(abs(o.y)>7) {
                    cout << "!" << i <<endl;
                    nag.tri[i].c[0].print();
                    nag.tri[i].c[1].print();
                    nag.tri[i].c[2].print();
                }*/
                ofile << o.x << "\t" << o.y << "\t" << o.z << endl;// cout << endl;
            }
        }
    }
    ofile.close();
}


/*

int Nagata::PrepareCoefficientVector()
{
    coeffs.clear();
    if(reducedNormals.size()!=fv->vertices.size())
        PrepareReducedNormals();

    for(int i=0;i!=fv->vertices.size();i++)
    {
        vector<Coefficient> pervertex(5);
        for(int j=0;j!=5;j++)
        {
            pervertex[j].type=j;
            pervertex[j].vertex=i;
            if(j<3) pervertex[j].val=(fv->vertices[i].pos.p+j);
            else
                pervertex[j].val=reducedNormals[i].p+j-3;
        }
        // Really clumsy, but I don't care once is enough.
        // Type 0-2 vertex,xyz 3,4 : Normal phi,theta
        std::list<FrozenCoeff>::iterator it;
        list<int> destroy;
        for (it = frozen.begin(); it != frozen.end(); ++it)
        {
            if ( (*it).vertex==i)
            {
                if(!(*it).isnormal)
                {
                    destroy.push_back((*it).dim);
                }
                else destroy.push_back((*it).dim+3);
            }
        }
        if(destroy.size()==0)
        {
            coeffs.insert(coeffs.end(),pervertex.begin(), pervertex.end());
        }
        else
        {
            cout << "Apply Constraint at " ;
            vector<Coefficient> pv2; pv2.clear();
            std::list<int>::iterator it;
            for(int k=0;k!=5;k++)
            {
                bool ban=false;
                for (it = destroy.begin(); it != destroy.end(); ++it)
                if( (*it) == k) ban=true;
                if(!ban) pv2.push_back(pervertex[k]);
                else cout<< "Vertex " << i+1 << " Coeff "  <<  k+1<< " ";
            }
        cout << endl;
        coeffs.insert(coeffs.end(),pv2.begin(), pv2.end());
        }
    }
    return coeffs.size();
}
*/

int Nagata::Nagata3Face(vector<Vertex*> &v, vector<Edge*> &e , NTriangle &nt)
{
    /* Transparency :)
     */
    if(v.size()!=3) return -1;
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

int Nagata::FromMesh(FaceVertexMesh& fvm)
{
    tri.resize(fvm.faces.size());

    for(int i=0; i<fvm.faces.size();i++)
    {
        if (fvm.faces[i].vert.size() !=3)
        {
            cout << "So far only Triangles" << endl;
            return -1;
        }

        vector<Vertex*> vlist(3);
        vector<Edge*> elist(3);
//      indexlist exposed(3);
        for(int j=0;j<3;j++)
        {
            vlist[j]=&fvm.vertices[fvm.faces[i].vert[j]];
            elist[j]=&fvm.edges[fvm.faces[i].edges[j]];
            //exposed[j]=fvm.faces[i].vert[j];
        }
        Nagata3Face(vlist,elist,tri[i]);
    //   tri[i].exposedv=exposed; tri[i].exposedn=exposed;tri[i].faceID=i;
    }
    fv=&fvm;
    //removeInflectionPoints();
    return 0;
}

double Membrane::IndexEnergy(int i) // Does not contain constraints
{
    return MembE.eval(&tri[i]);
}


double Membrane::DeltaE(int i, vector<double> *coll) // Crucial for MC DIFF!
{
    if(i==-1)
    {
    double dElam=0;
            for(int j=0;j<deepProps.size();j++)
            {
        double Tdx=membp.lambdas[deepPropIndex[j]]-MembE.defaults[j]; // Potential for Error here!
        double lambda=MembE.lambdas[j];
        dElam+=lambda*coll->at(j)*(0.5*coll->at(j)+Tdx);
 //       cout << j << " l" << lambda << " X0 " << MembE.defaults[j] << " X " << membp.lambdas[deepPropIndex[j]] << "DELTA X"  << coll->at(j)  << endl;
            }
    return dElam;
    }
    //DEBUG 
    double newE=MembE.deep_eval(&tri[i],true); // Create Trial Energy
    if(coll->size()!=deepProps.size()) {coll->resize(deepProps.size(),0); }
    double oldE=tri[i].GetProperty("Energy");
    // cout << " OLD E" << oldE << "New E" << newE << endl;
    // Evaluate Lambda-DeltaE
    for(int j=0;j<deepProps.size();j++)
    {   
        double x0=(tri[i].GetProperty(deepProps[j]));
        double x1=(tri[i].GetProperty("Try_"+deepProps[j]));
        coll->at(j)+=x1-x0;
      //  cout << "Volume" << tri[i].GetProperty("Volume") << endl;
      //  cout << deepProps[j] << " x0" <<x0 << " x1 "<< ("Try_"+deepProps[j]) << " " << x1<<" dx " << dx << endl;
      // 
   //     dElam+=lambda*dx*((x1+x0)*0.5-dx*MembE.defaults[j]);
   //    cout << "dElam " << dElam << " dx " << dx<<  endl;

        // For each Lambda
    } 
    //   cout << " E1 " << oldE <<" E2 "<< newE <<" dEner "<< newE-oldE << " dElam" << dElam <<endl;
    return (newE-oldE);
}

double Membrane::DExtraE (int i)
{
   // cout << "EXTRATERM CALLED" << endl;
    if(i == 0 or i>=extraterms.size()) {cerr << "Extraterm Index out of bound!" << endl; return 0;}
double  oldE=extraterms[i]->stored;
double   newE=extraterms[i]->eval(true);
 return (newE-oldE);
}

bool Membrane::CacheWrite(indexlist &cache) // Nodes are in Fact Triangle Indices
{
    vector <double> dx(deepProps.size(),0);
    //double dE=0;
    for(int i=0;i<cache.size();i++)
    {
        if(cache[i]>=0)
        {
                   for(int k=0;k!=deepProps.size();k++)
           {

           double xn=tri[cache[i]].GetProperty("Try_"+deepProps[k]);
           double x0=(tri[cache[i]].GetProperty(deepProps[k]));
               tri[cache[i]].SetProperty(deepProps[k],xn);
       //    cout << "CACHE xn" << xn << " X0 " << x0 << endl; 
           dx[k]+=(xn-x0);
      //     cout << "dx" << dx[k] << endl; 
            }
           double en=tri[cache[i]].GetProperty("Try_Energy");
      //     double e0=(tri[cache[i]].GetProperty("Energy");
      //     dE+=en-e0;
           tri[cache[i]].SetProperty("Energy",en);
        }
        else extraterms[-cache[i]]->store();
    }
   /* for (int i=0;i!=modnodes->size();i++)
    {
       for(int j=0;j!=modnodes->at(i).size();j++)
       {   
        // Try_Name --> Name
           for(int k=0;k!=deepProps.size();k++)
           {
           double xn=tri[modnodes->at(i)[j]].GetProperty("Try_"+deepProps[k]);
           double x0=(tri[i].GetProperty(deepProps[k]));
           cout << "CACHE xn" << xn << " X0 " << x0 << endl; 
           dx[k]+=(xn-x0);
           tri[modnodes->at(i)[j]].SetProperty(deepProps[k],xn);
           }
       }
    }*/
#pragma omp critical
{
    for(int i=0;i!=dx.size();i++)
        membp.lambdas[deepPropIndex[i]]+=dx[i]; // Potential for Error here!
}
return true;    
}


double Membrane::TotalEnergy(int stored=0)
{
    
    double Etot=0;
    if(!MembE.init) MembE.initializeValues(props);
    if(stored==0)
    {
    #pragma omp parallel for reduction(+: Etot)
    for(int i=0;i<tri.size();i++) Etot+=MembE.eval(&tri[i]);
   // cout << "TOTAL Energy " << Etot << endl;
    }
    if(stored==1)
    {
    #pragma omp parallel for reduction(+: Etot)
    for(int i=0;i<tri.size();i++) Etot+=tri[i].GetProperty("Energy");
   // cout << "Stored TOTAL Energy " << Etot << endl;
    }
    if(stored==2)
    { 
    Etot=MembE.computeTotalEnergy();
   // cout << "Constraint Total Energy " <<  Etot   << endl;
    }
    if(stored==4)
    { 
            #pragma omp parallel for reduction(+: Etot)
    for(int i=0;i<tri.size();i++) Etot+=tri[i].GetProperty("Energy")-MembE.eval(&tri[i]);
   // cout << "TOTAL Energy Deviation" << Etot << endl;

    }
    if(stored==5)
    {
            #pragma omp parallel for reduction(+: Etot)
    for(int i=1;i<extraterms.size();i++) Etot+=extraterms[i]->stored;
   // cout << "TOTAL Energy Deviation" << Etot << endl;

    }
    if(stored==6)
    {
            #pragma omp parallel for reduction(+: Etot)
    for(int i=1;i<extraterms.size();i++) Etot+=extraterms[i]->eval(false);
   // cout << "TOTAL Energy Deviation" << Etot << endl;

    }    

    return Etot;
}

double Membrane::SumStored(string key)
{
    double Etot=0;
    #pragma omp parallel for reduction(+: Etot)
    for(int i=0;i<tri.size();i++) Etot+=tri[i].GetProperty(key);
    cout << "Stored TOTAL "<< key << " : "  << Etot << endl;
    return Etot;
}
/*

double Nagata::TotalGaussCurvature()
{
    double totCurv = 0;
    for (int i =0; i < tri.size(); i++)
        totCurv += g.KInt(&tri[i]);
    cout << "TOTAL Gauss Curvature " << totCurv << endl;

    return totCurv;
}


double Nagata::TotalMeanCurvature()
{
    double totCurv = 0;
    #pragma omp parallel for reduction(+:totCurv)
    for (int i =0; i < tri.size(); i++)
        totCurv += g.HInt(&tri[i]);
    cout << "TOTAL Mean Curvature " << totCurv << endl;

    return totCurv;
}

*/

bool Nagata::AddMidpoint(int f, int v1, int v2)
{
        bool success=false;
        int verid1,verid2;
        verid1=fv->faces[f].vert[v1];
        verid2=fv->faces[f].vert[v2];        
        Edge *e;//=fv->edges[fv->faces[f].edges[0]];
        for(int i=0;i!=fv->faces[f].edges.size();i++)
        {
            e=&fv->edges[fv->faces[f].edges[i]];
            if((e->vertA==verid1 and e->vertB==verid2) or (e->vertB==verid1 and e->vertA==verid2))
            {success=true; break;}
        }
        if(!success) {cout << "FAIL!" << endl; return false;}
        // find midpoint of the (quadratic) edge
        if (v1 == 0)
        {
            if (v2 == 1)
            {
                e->midpoint.pos = tri[f].eval(0.5, 0.0);
                if (tri[f].n.size() == 6 && !tri[f].healthyEdge[0])       // advanced nagata triangle
                    e->midpoint.norm = *tri[f].n[3]*2.0;
                else
                    e->midpoint.norm = (tri[f].normal(0.5, 0.0))*(1+!tri[f].healthyEdge[0]);
            }
            else
            {
                e->midpoint.pos = tri[f].eval(0.5, 0.5);
               if (tri[f].n.size() == 6 && !tri[f].healthyEdge[2])       // advanced nagata triangle
                      e->midpoint.norm = *tri[f].n[5]*2.0;
                else
                    e->midpoint.norm = tri[f].normal(0.5, 0.5)*(1+!tri[f].healthyEdge[2]);
            }
        }
        else
        {
            e->midpoint.pos = tri[f].eval(1.0, 0.5);
            if (tri[f].n.size() == 6 && !tri[f].healthyEdge[1])       // advanced nagata triangle
               e->midpoint.norm = *tri[f].n[4]*2.0;
            else
                e->midpoint.norm = tri[f].normal(1.0, 0.5)*(1+!tri[f].healthyEdge[1]);
        }
//        e->midpoint.pos.print();
    //    fv->edges.push_back(e);
        return true;
}

bool Membrane::setBoundaryEnergy(CustomDepEnergyFactory &cde, bool forfreeze=false)
{
    bool block;
     cout << "Detecting Mesh Boundary" << endl;
     //0. Create EdgeEnergyPropertiesList.
    if(tri.size()==0) return false;
    if(edgep.theProps.size()==0) return false;
    Property *ep;
     //1 Detect Boundary faces.
    for(int i=0;i!=edgep.theProps.size();i++)
    {
       int face=0; FaceEnerProperty *fe;
       for(int j=0;j!=edgep.theProps[i]->deps.size();j++)
       { if(edgep.theProps[i]->deps[j]->ptype->collidable) { face++; fe=(FaceEnerProperty*) edgep.theProps[i]->deps[j];
                       
                                                            }
       }
       if(edgep.theProps[i]->block!=NULL) block=true;
           else block=false;
       if((face==1 and !forfreeze) or (block and forfreeze)) // Dangling Edge Detected. 
       {
       // cout << "DANGLING FACE" << endl ;
     //1.0 Detect Edge index of Triangle
        int aindex=fv->edges[edgep.theProps[i]->index].vertA;
        triple *v1=&fv->vertices[aindex].pos;
        int bindex=fv->edges[edgep.theProps[i]->index].vertB;
        triple *v2=&fv->vertices[bindex].pos;
        int faceind=-1;
        for(int l=0;l!=3;l++)
        {
          if(tri[fe->index].v[l]==v1)
          {faceind=l; break;}
        }
        if(faceind==-1) cerr << "Edge not found" << endl;
     //1.a Generate Properies as Necessary.
    //    cout << "Generating Properties" << endl ;
           if(extraterms.size()==0) extraterms.push_back(NULL); // Zero must not be part
           ep=new Property;
           ep->parent=(void*) edgep.theProps[i];
           ep->ptype=&eep;
           ep->block=NULL;
           ep->deps=edgep.theProps[i]->deps;
           ep->index=-extraterms.size();
           ep->state=faceind;
           ep->content=edgep.theProps[i]->content; //(double*) cd;
           CustomDepEnergy *cd= cde.GetNew((void*) ep);
           extraterms.push_back(cd);
           eep.theProps.push_back(ep);
           edgep.theProps[i]->deps.push_back(ep);
           
    //2. Hook Into Evaluation Structure.
     //             for(int a=0;a!=edgep.theProps[i]->deps.size();a++)
     //       cout << a << " " << edgep.theProps[i]->deps[a]->ptype->name << endl; 
          faceind=fe->index;
        //cout << "Hook in " << faceind << " // " << fv->faces[faceind].vert.size()<< endl ;

    for(int k=0; k!= fv->faces[faceind].vert.size(); k++)
       {
        int vnum=fv->faces[faceind].vert[k];
   //     cout << vnum << " " <<endl;
        Property *v=vertp.theProps[vnum];
        Property *n=normp.theProps[vnum];
        v->deps.push_back(ep);
        n->deps.push_back(ep);
       }
       }
    }
    if(extraterms.size()>0) {
        cout << "BOUNDARY EDGE SIZE " << extraterms.size()-1 << endl;
        if(!forfreeze) return ChainBoundary();
        else return true;
    }
    else cout << "SYSTEM IS BOUNDARY FREE" << endl;
    return true;
}

bool Membrane::ChainBoundary()
{
    //1. Iterate through EEP.
    // Each Edge Index has only one partner. We track the used ones.
    // More than one cycles is possible.
    set<int> done;
    for(int i=0;i!=eep.theProps.size();i++)
    {
        int index=eep.theProps[i]->index;
        int prevB=((Edge*) (eep.theProps[i]->content))->vertB;
        if(done.find(index) == done.end())
        {
            int next=0; // Zero is a forbidden index
            while(index!=next)
            {
                // We go from index B to index A.
                for(int j=i;j!=eep.theProps.size();j++)
                {
                    Edge *e=(Edge*) eep.theProps[j]->content;
                    if(e->vertA==prevB) 
                    {
                        if(next==0) extraterms[-index]->next=(void*) eep.theProps[j];
                        else extraterms[-next]->next=(void*) eep.theProps[j];
                        next=eep.theProps[j]->index;
                        done.insert(next);
                        prevB=e->vertB;
                    }
                }
                if(done.find(index)!=done.end()) break;
            }
        }
    }
    return true;
}


bool Nagata::BuildMidpoints()
{
        /**
     *
     * This function is used for writing results in VTU format with curved elements.
     * For this we need to know the positions of the midpoints of the edges
     * (as the edges are curved). This is what this function takes care of.
     *
     */

    for (int f = 0; f < fv->faces.size(); f++)
    {
        AddMidpoint(f, 0, 1);
        AddMidpoint(f, 1, 2);
        AddMidpoint(f, 0, 2);
    }
    return true;
}

 bool Nagata::exposeNeighbors(vector<Neighbor*> &col)
{
  col.resize(tri.size());
  for(int i=0;i!=tri.size();i++) col[i]=&tri[i];
  return true;
}


bool Nagata::WriteVtu(const string &filename, bool quad = false)
{
    ofstream ofile;
    ofile.open(filename.c_str());

    // VTK header
    ofile << "<?xml version=\"1.0\"?>\n";
    ofile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    ofile << "\t<UnstructuredGrid>\n";

    if (quad)
    {
        BuildMidpoints();
        ofile << "\t\t<Piece NumberOfPoints=\"" << fv->vertices.size() + fv->edges.size() << "\" NumberOfCells=\"" << fv->faces.size() << "\">\n";
    }
    else
        ofile << "\t\t<Piece NumberOfPoints=\"" << fv->vertices.size() << "\" NumberOfCells=\"" << fv->faces.size() << "\">\n";

    // POINT section
    ofile << "\t\t\t<Points>\n";
    ofile << "\t\t\t\t<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto vert = fv->vertices.begin(); vert != fv->vertices.end(); vert++)
    {
        ofile << "\t\t\t\t\t" << vert->pos.x << " " << vert->pos.y << " " << vert->pos.z << "\n";
    }
    if (quad)
        for (auto & e : fv->edges)
          ofile << "\t\t\t\t\t" << e.midpoint.pos.x << " " << e.midpoint.pos.y << " " << e.midpoint.pos.z << "\n";

    ofile << "\t\t\t\t</DataArray>\n";
    ofile << "\t\t\t</Points>\n";

    // FACE section
    ofile << "\t\t\t<Cells>\n";

    ofile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    ofile << "\t\t\t\t\t";
    for (auto f = fv->faces.begin(); f != fv->faces.end(); f++)
    {
        for (auto v = f->vert.begin(); v != f->vert.end(); v++)
        {
            ofile << *v << " ";
        }
        if (quad)
        {
            int en[3];
            for (auto & e : f->edges)
            {
                if(fv->edges[e].vertA==f->vert[0] and fv->edges[e].vertB==f->vert[1]) en[0]=e;
                if(fv->edges[e].vertA==f->vert[1] and fv->edges[e].vertB==f->vert[0]) en[0]=e;
                if(fv->edges[e].vertA==f->vert[1] and fv->edges[e].vertB==f->vert[2]) en[1]=e;
                if(fv->edges[e].vertA==f->vert[2] and fv->edges[e].vertB==f->vert[1]) en[1]=e;
                if(fv->edges[e].vertA==f->vert[0] and fv->edges[e].vertB==f->vert[2]) en[2]=e;
                if(fv->edges[e].vertA==f->vert[2] and fv->edges[e].vertB==f->vert[0]) en[2]=e;                
            }
            ofile << fv->vertices.size() + en[0] << " ";
            ofile << fv->vertices.size() + en[1] << " ";
            ofile << fv->vertices.size() + en[2] << " ";
        }
    }
    ofile << "\n\t\t\t\t</DataArray>\n";

    ofile << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    ofile << "\t\t\t\t\t";
    int offset = 0;
    for (auto f = fv->faces.begin(); f != fv->faces.end(); f++)
    {
        if (quad)
            offset += 6;
        else
            offset += 3;
        ofile << offset << " ";
    }
    ofile << "\n\t\t\t\t</DataArray>\n";

    ofile << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    ofile << "\t\t\t\t\t";
    string code = (quad ? "69 " : "5 ");
    for (auto f = fv->faces.begin(); f != fv->faces.end(); f++)
        ofile << code;
    ofile << "\n\t\t\t\t</DataArray>\n";

    ofile << "\t\t\t</Cells>\n";

    ofile << "\t\t\t<PointData>\n";
   
    ofile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"normal\" NumberOfComponents=\"3\" "
        << "ComponentName0=\"x\" ComponentName1=\"y\" ComponentName2=\"z\" "
        << "format=\"ascii\">\n";
    for (auto v = fv->vertices.begin(); v != fv->vertices.end(); v++)
    {
        ofile << "\t\t\t\t\t" << v->norm.x << " " << v->norm.y << " " << v->norm.z << "\n";
    }
    if (quad)
        for (auto & e : fv->edges)
        {
            ofile << "\t\t\t\t\t" << e.midpoint.norm.x << " " << e.midpoint.norm.y << " " << e.midpoint.norm.z << "\n";
        }
    ofile << "\t\t\t\t</DataArray>\n";
    ofile << "\t\t\t</PointData>\n";

    // CELL DATA section
    ofile << "\t\t\t<CellData>\n";
   
    vector< pair <string, double > > vds;
    auto f=tri.begin(); f->ExportProperties(&vds);
    for(int i=0;i!=vds.size();i++)
    {
    // area
    ofile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"" << vds[i].first << "\" NumberOfComponents=\"1\" ComponentName0=\"YesWeCanEvenNameIt\" "
        << "format=\"ascii\">\n";
    for (auto f = tri.begin(); f != tri.end(); f++)
    {
        ofile << "\t\t\t\t\t" << f->GetProperty(vds[i].first) << "\n";
        // ofile << "\t\t\t\t\t" << f->area << "\n";
    }
    ofile << "\t\t\t\t</DataArray>\n";
    }
    
    // Gauss curvature per face
    ofile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"gaussCurvature\" NumberOfComponents=\"1\" ComponentName0=\"gaussCurvature\" "
        << "format=\"ascii\">\n";
    for (auto f = tri.begin(); f != tri.end(); f++)
    {
        ofile << "\t\t\t\t\t" << g.KInt(&(*f)) / g.Area(&(*f)) << "\n";
    }
    ofile << "\t\t\t\t</DataArray>\n";
    // Mean curvature per face
    ofile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"meanCurvature\" NumberOfComponents=\"1\" ComponentName0=\"meanCurvature\" "
        << "format=\"ascii\">\n";
    for (auto f = tri.begin(); f != tri.end(); f++)
    {
        ofile << "\t\t\t\t\t" << g.HInt(&(*f)) / g.Area(&(*f)) << "\n";

    }
    ofile << "\t\t\t\t</DataArray>\n";
    // Mean curvature per face
    ofile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Energy_Dens\" NumberOfComponents=\"1\" ComponentName0=\"meanCurvature\" "
        << "format=\"ascii\">\n";
    for (auto f = tri.begin(); f != tri.end(); f++)
    {
        ofile << "\t\t\t\t\t" << f->GetProperty("Energy")/ g.Area(&(*f)) << "\n";
    }
    ofile << "\t\t\t\t</DataArray>\n";

    ofile << "\t\t\t\t<DataArray type=\"Float32\" Name=\"Enorm\" NumberOfComponents=\"3\" "
    << "ComponentName0=\"x\" ComponentName1=\"y\" ComponentName2=\"z\" "
        << "format=\"ascii\">\n";
    for (auto f = tri.begin(); f != tri.end(); f++)
    {
        triple ne=(*f).normal(0.66666666666667,0.33333333333333);
        ofile << "\t\t\t\t\t" << ne.x << " " << ne.y << " " << ne.z << "\n";
    }
    ofile << "\t\t\t\t</DataArray>\n";

    ofile << "\t\t\t</CellData>\n";


    // VTK ending
    ofile << "\t\t</Piece>\n";
    ofile << "\t</UnstructuredGrid>\n";
    ofile << "</VTKFile>\n";

    return true;
}
/*
  Looks horrible, but is it ?
*/


#ifdef HAVE_GRADIENT
#include"agradient.hpp"
#endif

