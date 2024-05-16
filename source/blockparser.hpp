

bool ParseBlocks(string blockfile, Entity *ent)
{
    ifstream hfile;
    hfile.open(blockfile.c_str());
    if(hfile.fail()) return false;
    string current;
    while ( !hfile.eof() )
    {
        vector<string> line;
        getline(hfile, current);
        int htokens=Tokenize(current,&line, " \t\r");
        if(htokens<3)
            continue;
        int index=atoi(line[1].c_str());
        if(line[0]=="Edge")
        {
          // There will be a special format Edge Blocking has no dimension, index is too hard to guess,
          // So the Format is Edge AtomA AtomB
           PropertyList *edgelist;
           if(ent->props.GetPropertyList(line[0],&edgelist))
           {
           bool found=false;
           for(int i=0;i!=edgelist->theProps.size();i++)
           {
           Edge* e= (Edge*) edgelist->theProps[i]->content ;   
           if(e->vertA==atoi(line[1].c_str()) && e->vertB==atoi(line[2].c_str())) found=true;
           if(e->vertA==atoi(line[2].c_str()) && e->vertB==atoi(line[1].c_str())) found=true;
           if(found) { index=i;cout << "BLOCK: Translated Vertices " << e->vertA << " " << e->vertB << " to Index of Edge: " << index << endl; break;}
           }
           }
           else cerr << "BLOCK Error: Entity has no Edges" << endl;
        }
        // Format Name, Index Prop, Index_Field;
        //cout << line[1] << " "<< line[2] << endl;
        Block b;b.dim.resize(1);
        if(line[0]!="Edge")
        {
        if(line[2]=="*") 
        {
            b.dim.resize(3); for(int i=0;i!=3;i++) b.dim[i]=i;
        }
        else b.dim[0]=atoi(line[2].c_str());
        }
        b.origin.resize(b.dim.size());
        cout << "Adding Block to " << line[0] << " Index " << index << " ID " << b.dim[0] << endl;
        ent->AddBlock(b,line[0],index);
    }
    hfile.close();
    return true;
}
