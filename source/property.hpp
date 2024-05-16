typedef vector<int> indexlist;
struct EvalStruc 
{
    
};

struct Property;
class PropertyList;
typedef vector<Property*> Deplist;


struct Block
{
  void applyBlock(double *content)
  {
      for(int i=0;i!=dim.size();i++)
           content[dim[i]]=origin[i];    
  }
  bool merge(Block &b)
  {
     bool overlap=false;            
     indexlist overlist(0);
     for(int i=0;i!=dim.size();i++)
      {
            for(int j=0;j!=b.dim.size();j++)
                if(dim[i]==b.dim[j]) {overlap=true; overlist.push_back(b.dim[j]);}            
      }
     if(!overlap) {origin.insert(origin.end(),b.origin.begin(),b.origin.end());dim.insert(dim.end(),b.dim.begin(),b.dim.end()); return true;}
     for(int i=0;i!=b.dim.size();i++)
     {         overlap=false;
               for(int j=0;j!=overlist.size();j++) if(overlist[j]==b.dim[i]) overlap=true;
               if(!overlap) {dim.push_back(b.dim[i]);origin.push_back(b.origin[i]);}
     }             
     return true;
  }
  vector<double> origin;
  indexlist dim;
};


class PropertyList // : Structure;
{
public:
  PropertyList()
  {collidable=false;energetic=false;}
  int szProperty;
  string name;
  bool collidable;
  bool energetic;
  virtual bool genpropertiesfromStructure(void *parent) {return true;}
  //  -->
  // --> Access Property_of_Structure.
  int pofs(int index, Property **p);
    // --> Access Structure of Property.
  int sofp(int index,void **struc);
  int coeffs()
  {
      return 1;
  };
  void *parent;
  vector<Property*> theProps;
  ~PropertyList()
  {
  for(int i=0;i!=theProps.size();i++) delete theProps[i];
  theProps.clear();
  }
  // --> 
  //
/*protected:
    PropertyList()
    {
    }
*/
};

struct Property
{   
    Property(){deps.clear();}
    PropertyList *ptype; // Type Info
    void *parent; // Link to Parent Object
    int type;
    int index;
    Deplist *dependents() { return &deps;}
    Deplist deps;
    virtual bool update(bool deep=true) {
    //    cout << "DEFAULT" << endl;
        if(block!=NULL) block->applyBlock(content); 
        if(deep ) 
         for(int i=0;i!=deps.size();i++) deps[i]->update();
        
    return true;
    }
    bool collapse(vector<double> &ref)
       {
          // if(block==NULL)
          memcpy(&ref, content, (ptype->szProperty)*sizeof(double));
          /* else
           {
           int edim=(ptype->szProperty)-block->dim.size();
           vector<double> temp(edim);
           int j=0,k=0;
           for(int i=0;i!=ptype->szProperty;i++)
           {
            if(block->dim[j]!=i) {temp[k]=*(content+i);k++;}
            else j++;
           }
           memcpy(&ref, temp, ptype->szProperty*sizeof(double));
           }*/
           return true;
       }
    double *content;
    Block *block;
    short int state;
};

int PropertyList::pofs(int index, Property** p)
{
      *p=theProps[index];
      return 0;
}

int PropertyList::sofp(int index, void** struc)
{
      *struc=theProps[index]->parent;
      return 0;
}


// Structure

// class PropertiedMesh : Mesh, PropertyStructure;



class PropertyMap
{
public:
 bool AddPropertyType(PropertyList &pl);
 bool GetPropertyList(string &key, PropertyList **pl);
 Property &GetProperty(const string &key, int index, bool update=false);
 map<string, PropertyList*> PropertyL;
 bool AvailProp(vector <string> &keys);
 bool HasProp(const string &pname);
};

typedef list<Block> Blocklist;
 
struct Propertied
{
    PropertyMap props;
    Blocklist blocks;
};

bool PropertyMap::AddPropertyType(PropertyList &pl)
{
    auto rv=PropertyL.insert(make_pair(pl.name,&pl));
     if(rv.second==false)
    {
        cout << "PropertyList existed - Overwriting it!" << endl;
        PropertyL.erase(rv.first);
        PropertyL.insert(make_pair(pl.name,&pl));
    }
   
    return true;
}

bool PropertyMap::GetPropertyList(string &key, PropertyList **pl)
{
    *pl=(PropertyL[key]);
    return true;
}

Property &PropertyMap::GetProperty(const string &key,int index, bool update )
{
   if(update) PropertyL[key]->theProps[index]->update();
 //  cout << *PropertyL[key]->theProps[index]->content << endl; 
   return *(PropertyL[key]->theProps[index]);
 //  return true;
}

bool PropertyMap::HasProp(const string &key)
{
   if (PropertyL.find(key)==PropertyL.end()) return false;
   return true;
}

bool PropertyMap::AvailProp(vector<std::string> &keys)
{
    keys.resize(PropertyL.size());
    int i=0;
    for(map<string, PropertyList*>::iterator it=PropertyL.begin();it!=PropertyL.end();it++)
    {
    keys[i]=it->first;
    i++;
    }
    return true;
}

// Merges Deplist2 into Deplist 1 - Does not go beyond first level! You can do it recursively.

bool MergeDeplists(Deplist *list1, Deplist* list2)
{
    set<Property*> bag(list1->begin(),list1->end());
 
    std::copy(list2->begin(), list2->end(), std::inserter(bag, bag.end()));
    list1->assign(bag.begin(), bag.end());
//    std::copy(list2->begin(), list2->end(), std::inserter(bag, bag.end()));
    return true;
}

bool FullMergeDeplists(Deplist *goal, Deplist *traverse, bool uniq=false)
{
    if(traverse==NULL) return false;
    for(int j=0;j!=traverse->size();j++)
    {
        Deplist inter; inter.clear();
        //if(traverse->at(j).dependents!=0)
        if(FullMergeDeplists(&inter,&(traverse->at(j))->deps,false))
            MergeDeplists(goal,&inter);
    }
    MergeDeplists(goal,traverse);
    if(uniq)
    {
        sort( goal->begin(), goal->end() );
       goal->erase( std::unique( goal->begin(), goal->end() ), goal->end() );
    }
    return true;
}




/*
 * Property: Get
 * 
 * 
 * 
 * 
 */
/*
 * Energy = Sum Int(Face) + Sum(Particles) + Sum Sum (Face/Particles) + lambda*Const(Sum(Face)+Sum Sum)
 * Int(Face(Vert,Mix,Norm,Props)) 
 */
/// Gradient Stuct:
// flowest(h1, h2... hn) -->    Fhigh(f

/// Updating scheme
// Mesh (scale) --> Vertex (scale,update) --> Face (update)
// Chain Rule:
// --> Elements.  
/// Derivative Scheme
