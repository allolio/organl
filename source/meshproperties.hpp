

class NormalProperty;

class NormalProperties: public PropertyList
{
public:
    NormalProperties()
    {   name="Normal";
        szProperty=2;
    }
    vector<duple> reducedNormals;
    virtual bool genpropertiesfromStructure(void *parent);
    duple normang(triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    triple angnorm(duple *ang)
    {
        triple rn;
        rn.x=cos(ang->theta)*sin(ang->phi);
        rn.y=sin(ang->theta)*sin(ang->phi);
        rn.z=cos(ang->phi);
        return rn;
    }
};

class NormalProperty: public Property
{
public:
    NormalProperties *ptype;
    bool update(bool deep=true) {
        if(block!=NULL)  { block->applyBlock(content);
        Vertex *v=(Vertex*) parent;
        v->norm=ptype->angnorm((duple*) content);}
//        else content=normang(v->norm);
        if(deep)
            for(int i=0; i!=deps.size(); i++) deps[i]->update();
        return true;
    }
};

bool NormalProperties::genpropertiesfromStructure(void* parent)
{       
    this->parent=parent;
    FaceVertexMesh *fv=(FaceVertexMesh*) parent;
    theProps.resize(fv->vertices.size());
    reducedNormals.resize(fv->vertices.size());
    for(int i=0; i!=fv->vertices.size(); i++)
    {
//            NormalProperty n;
        theProps[i]=new NormalProperty;
        theProps[i]->block=NULL;
        theProps[i]->parent=&fv->vertices[i];
        theProps[i]->deps.clear();
        theProps[i]->state=0;
        theProps[i]->ptype=this;
        theProps[i]->index=i;
        reducedNormals[i]=normang(fv->vertices[i].norm);
        theProps[i]->content=reducedNormals[i].p;
    }
    return true;
}



class VertexProperties: public PropertyList
{
public:
    VertexProperties()
    {   name="Vertex";
        szProperty=3;
    };
    virtual bool genpropertiesfromStructure(void *parent)
    {
        FaceVertexMesh *fv=(FaceVertexMesh*) parent;
        this->parent=parent;
        theProps.resize(fv->vertices.size());
        for(int i=0; i!=fv->vertices.size(); i++)
        {
            theProps[i]=new Property;
            theProps[i]->block=NULL;
            theProps[i]->parent=&fv->vertices[i];
            //theProps[i].deps=NULL;
            theProps[i]->index=i;
            theProps[i]->state=0;
            theProps[i]->ptype=this;
            theProps[i]->content=fv->vertices[i].pos.p;
        }
        return true;
    }
};    

class EdgeProperty;


class EdgeProperties: public PropertyList
{
public:
    EdgeProperties()
    {   name="Edge";
        szProperty=3;
   //     collidable=false;
   //     energetic=false;
    }
    virtual bool genpropertiesfromStructure(void *parent);

};

class EdgeProperty: public Property
{
public:
    EdgeProperties *ptype;
    bool update(bool deep=false) { // Dont Dep updating by default.
//        cout << "Called Edge UPDATE" << deep << endl;
  //      if(deep)
    //        for(int i=0; i!=deps.size(); i++) { cout <<  deps[i]->ptype->name << endl; deps[i]->update();}
  //      if(block!=NULL)  { block->applyBlock(content);}
  //      if(deep)
    //        for(int i=0; i!=deps.size(); i++) deps[i]->update();
        return true;
    }
};


bool EdgeProperties::genpropertiesfromStructure(void* parent)
{
    
        this->parent=parent;
        FaceVertexMesh *fv=(FaceVertexMesh*) parent;
        theProps.resize(fv->edges.size());
        for(int i=0; i!=fv->edges.size(); i++)
        {
            theProps[i]=new Property;
            theProps[i]->block=NULL;
            theProps[i]->parent=&fv->edges[i];
            //theProps[i].deps=NULL;
            theProps[i]->index=i;
            theProps[i]->state=0;
            theProps[i]->ptype=this;
            theProps[i]->content=(double* ) &fv->edges[i];
        }
        return true;
    
}


class FaceEnerProperty;

class FaceEnerProperties: public PropertyList
{
public:
    FaceEnerProperties()
    {   name="FaceE";
        collidable=true;
        energetic=true;
        szProperty=sizeof(void*);
    }
    virtual bool genpropertiesfromStructure(void *parent);
};

class FaceEnerProperty: public Property
{
public:
    FaceEnerProperties *ptype;
    bool update(bool deep=true) {
    //    if(block!=NULL) block->applyBlock(content); NO BLOCKS
      //  Vertex *v=(Vertex*) parent;
        NTriangle *nag=(NTriangle*) content;
   //     cout << "UPLDATE" << endl;
    //    if(deep) NO DEPS
    //        for(int i=0; i!=deps.size(); i++) deps[i]->update();
        return nag->updatec();
    }
};

 bool FaceEnerProperties::genpropertiesfromStructure(void* parent)
    {
        vector< NTriangle > *tri=(vector < NTriangle > *) parent;
        theProps.resize(tri->size());
        for(int i=0;i!=tri->size();i++)
        {
            theProps[i]=new FaceEnerProperty;
            theProps[i]->block=NULL;
            theProps[i]->parent=NULL;
            theProps[i]->deps.clear();
            theProps[i]->index=i;
            theProps[i]->state=0;
            theProps[i]->ptype=this;
            theProps[i]->content=(double *) &tri->at(i);
        }
        return true;
    }
    
    
    
class EdgeEnerProperty;

class EdgeEnerProperties: public PropertyList
{
public:
    EdgeEnerProperties()
    {   name="EdgeE";
        collidable=false;
        energetic=true;
        szProperty=sizeof(void*);
    }
    virtual bool genpropertiesfromStructure(void *parent);
};

class EdgeEnerProperty: public Property
{
public:
    EdgeEnerProperties *ptype;
    bool update(bool deep=true) {
    //    if(block!=NULL) block->applyBlock(content); NO BLOCKS
      //  Vertex *v=(Vertex*) parent;
        return true;
    }
};

 bool EdgeEnerProperties::genpropertiesfromStructure(void* parent)
    {
        // No automatic generation
        return true;
    }

