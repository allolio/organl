//#define EDGE_TILT_LIMIT 0.25//0.0 //0.25! // 
typedef vector<int> indexlist;

struct Vertex {
    public:
        triple pos;
        triple norm;
        indexlist neighbors;
};

struct Edge {
    public:
        int vertA;
        int vertB;
        triple c;
        Vertex midpoint;
};

struct Face {
    indexlist vert;
    indexlist edges;
};



class FaceVertexMesh {
    public:
        vector<Vertex> vertices;
        vector<Face> faces;
        vector<indexlist> vl;               // List of faces for a vertex
        vector<Edge> edges;
        bool LoadObj(const string &filename);
        bool SaveObj(const string &filename);
        bool BuildNormals();
        bool ResetNormals();
        bool BuildEdges();
        bool BuildNormal(int center,bool orient=false, triple* refvec=NULL);
        double QualityAroundVertex(int center);

    private:
        void AddEdge(int f, int v1, int v2);
        bool consistentNormals(int);
        bool vlFromFaces();
        bool inlFromvl();
};

void FaceVertexMesh::AddEdge(int f, int v1, int v2)
{
    bool newEdge = true;
    for (int i = 0; i < edges.size() and newEdge; i++)
    {
        if ((edges[i].vertA == faces[f].vert[v1] and edges[i].vertB == faces[f].vert[v2])
                or (edges[i].vertA == faces[f].vert[v2] and edges[i].vertB == faces[f].vert[v1]))
        {
            // edge already exists
            newEdge = false;
            faces[f].edges.push_back(i);
        }
    }
    if (newEdge)
    {
        faces[f].edges.push_back(edges.size());
        Edge e;
        e.vertA = faces[f].vert[v1];
        e.vertB = faces[f].vert[v2];
        edges.push_back(e);
    }
}

bool FaceVertexMesh::BuildEdges()
{

    edges.clear();
    for (int f = 0; f < faces.size(); f++)
    {
        faces[f].edges.clear();
        AddEdge(f, 0, 1);
        AddEdge(f, 1, 2);
        AddEdge(f, 0, 2);
    }

    return true;
}


bool FaceVertexMesh::vlFromFaces() {
    vl.clear();
    vl.resize(vertices.size());

    for(int i=0; i != faces.size(); i++) {
        for (int j=0;j!=faces[i].vert.size();j++) {
            vl.at(faces[i].vert[j]).push_back(i);
        }
    }

    return true;
}


bool FaceVertexMesh::inlFromvl() {
    for (int i=0;i!=vertices.size();i++) {
        set<int> nli;

        for (int j=0; j!=vl[i].size();j++) {
            for (int k=0;k!=faces[vl[i].at(j)].vert.size();k++) {
                if(faces[vl[i][j]].vert[k]!=i)
                nli.insert(faces[vl[i][j]].vert[k]);
            }
        }
        vertices[i].neighbors.insert(vertices[i].neighbors.begin(), nli.begin(),nli.end());
    }
    return true;
}

bool FaceVertexMesh::consistentNormals(int b=0)
{
    indexlist done;
    done.resize(vertices.size());
    // Built Upon the Assumption that adjacent normals cannot be more than 90°
    int i=b;
    done[b]=true;
    int last=b;
    auto glamb = [&](){int a=0; for (int x:done) a+=x; return a;};
    while( done.size()-glamb() > 0 ) {
        int change=0;
        for (int j=0; j!=vl[i].size(); j++) {
            for (int k=0; k!=faces[vl[i][j]].vert.size(); k++) {
                int nb;
                nb=faces[vl[i][j]].vert[k];
                if(nb!=i && done[nb]!=true) {
                    double d=vertices[nb].norm*vertices[i].norm;
                    if(d<0) {
                        vertices[nb].norm=vertices[nb].norm*-1;
                    }
                    done[nb]=true;
                    last=nb; change++;
                   // cout << " MOD " << nb <<  " " << glamb()<< endl;
                }
             }
        }

        i=last;
        if(change==0)
        {
            int nb=-1;
            for (int l=0;l!=done.size();l++) {
                if(done[l]==false) {
                    for (int j=0;j!=vl[l].size();j++) {
                        for (int k=0;k!=faces[vl[l][j]].vert.size();k++) {
                             int nob= faces[vl[l][j]].vert[k];
                             if (done[nob]==true) {i=nob; nb=nob; break;}
                        }
                        if(nb!=-1) break;
                    }
                }

                if(nb!=-1) break;

                if(l==done.size()-1) {
                    cout << "Disconnected Mesh! Guessing." << endl;
                    for (int x=0; x!=done.size(); x++) {
                        if(!done[x]) {
                            done[x]=true;nb=x;i=x; break;
                        }
                    }
                }
            }
        }
    }
    return true;
}

bool FaceVertexMesh::BuildNormals() {
    for (int i=0; i!=vertices.size(); i++) {
        BuildNormal(i);
    }
    return true;
}

bool FaceVertexMesh::SaveObj(const string &filename) {
    ofstream ofile;
    ofile.open(filename.c_str());

    for (int i=0; i!=vertices.size(); i++) {
        ofile << "v " ;
        ofile << vertices[i].pos.x << " " << vertices[i].pos.y<< " " << vertices[i].pos.z << " " << endl;
     //     cout << fvm.vertices[i].norm.x << " " << fvm.vertices[i].norm.y<< " " << fvm.vertices[i].norm.z << " " << endl;
     //      for (int j=0;j!=fvm.vl[i].size();j++) cout << fvm.vl[i][j]+1 <<" "; cout << endl;
     //        cout << " " << fvm.vl.size() << " " << endl;
    }

    for (int i=0; i!=vertices.size(); i++) {
        ofile << "vn " ;
        ofile << vertices[i].norm.x << " " << vertices[i].norm.y<< " " << vertices[i].norm.z << " " << endl;
     //     cout << fvm.vertices[i].norm.x << " " << fvm.vertices[i].norm.y<< " " << fvm.vertices[i].norm.z << " " << endl;
     //      for (int j=0;j!=fvm.vl[i].size();j++) cout << fvm.vl[i][j]+1 <<" "; cout << endl;
     //        cout << " " << fvm.vl.size() << " " << endl;
    }

    for (int i=0; i!=faces.size(); i++) {
        ofile << "f " ;
        for (int j=0; j!=faces[i].vert.size(); j++) {
            ofile << faces[i].vert[j]+1 << " ";
        }
        ofile << endl;
    }

    ofile.close();
    return 1;
}

double FaceVertexMesh::QualityAroundVertex(int center) {
/*
 * We now check for faces which are at angles of more than 90° angle, as this will result in a failure.
 */
#ifdef EDGE_TILT_LIMIT
    int numberOfFaces = vl[center].size();  // number of Faces around the central vertex 'center'
    vector<triple> store(numberOfFaces);

    for (int j=0; j!=numberOfFaces; j++) {

        int pop=0;
        triple a,b;

        for (int k=0; k!=faces[vl[center][j]].vert.size(); k++) {
            int cf = faces[vl[center][j]].vert[k];
            if(cf!=center) {
                if(pop==0) {
                    a=vertices[cf].pos;
                    pop++;
                }
                else {
                    b=vertices[cf].pos;
                    pop++;
                    break;
                }
            }
        }

        a = a - vertices[center].pos;
        b = b - vertices[center].pos;
        a = a;  /// a.abs();
        b = b;  /// b.abs(); // EDGE WEIGHTING
        triple n1 = a % b;
//      n1 = n1 / n1.abs(); // EQUAL WEIGHTING
        //double cosA = (a * b);
        store[j]=n1/n1.abs();
    }

    for (int i=0; i!=numberOfFaces; i++) {
        for (int j=i+1; j<numberOfFaces; j++) {
            double ab=store[i]*store[j];
            if(ab*ab<EDGE_TILT_LIMIT) return 0.0; //0.25 0.15
        }
    }
#endif
    return 1.0;
}

bool FaceVertexMesh::BuildNormal(int center,bool orient, triple* refvec) {
/**
 *
 * There are plenty of possible weightings
 *
 * The Visual Computer (2005) 21:71–82Digital Object Identifier (DOI) 10.1007/s00371-004-0271-1
 *
 * They claim that MWA approach is the best
 *
 */
    triple nvec;
    nvec=nvec*0;
    double orn=0;
    int j=0;
    /*
     * We now check for faces which are at angles of more than 90° angle this will result in a failure.
        */
    int numberOfFaces = vl[center].size();
    vector<triple> store(numberOfFaces);

    for (j=0; j!=numberOfFaces; j++) {
        int pop = 0;
        triple a,b;

        for (int k=0; k!=faces[vl[center][j]].vert.size(); k++) {
            int cf = faces[vl[center][j]].vert[k];
            if (cf != center) {
                if (pop == 0) {
                    a=vertices[cf].pos;
                    pop++;
                }
                else {
                    b=vertices[cf].pos;
                    pop++;
                    break;
                }
            }
        }

        a = a - vertices[center].pos;
        b = b - vertices[center].pos;
        a = a/a.abs();
        b = b/b.abs(); // EDGE WEIGHTING
        triple n1 = a % b;
        //n1 = n1 / n1.abs(); // EQUAL WEIGHTING
        //double cosA = (a * b);
        double w = 1;//((a*a)+(b*b));//acos(cosA);
        if (orient and refvec!=NULL)
            orn = n1 * (*refvec);
        if (orient and refvec==NULL)
            orn = n1 * vertices[center].norm;
        else
            orn = n1 * nvec;

        if (orn < 0) {
            n1 = n1 * -1;
        }
#ifdef EDGE_TILT_LIMIT
        store[j] = n1 / n1.abs();
#endif
        nvec = nvec + n1*w;
    }

    nvec = nvec / nvec.abs();
    if(orient && refvec!=NULL) {
        if(*refvec * nvec < 0) {
            nvec = nvec * -1;
        }
    }

    vertices[center].norm = nvec;
#ifdef EDGE_TILT_LIMIT
    for (int i=0; i!=numberOfFaces; i++) {
        for (int j=i+1; j<numberOfFaces; j++) {
            double ab = store[i]*store[j];
            if(ab*ab < EDGE_TILT_LIMIT) return false; //0.25 0.15
        }
    }
#endif
    return true;
}


/* Load OBJ files
 * Currently only vertex, vertex normal and face is supported
 */
bool FaceVertexMesh::LoadObj(const string &filename) {

    string currentLine;
    ifstream fileStream;
    int VertexCounter = 0;
    int NormalCounter = 0;
    bool hasNormals = false;
    int LineCounter = 0;

    // Open file
    fileStream.open( filename.c_str() );
    if ( fileStream.fail() ) {
        // Could not open file
        return false;
    }

    // Parse file contents
    while ( !fileStream.eof() ) {
        vector<string> line;
        getline(fileStream, currentLine);
        LineCounter++;
        int nTokens = Tokenize(currentLine,&line, " \t\r");

        // Skip comments and empty lines in input
        if (line.size()==0 or line[0][0] == '#') {
            continue;
        }

        // Every line must have at least 3 tokens (except comments)
        if(nTokens < 3) {
            cout << "Cannot parse broken mesh! Incorrect number of tokens in line " << LineCounter << endl;
            break;
        }

        // Vertex
        if(line[0] == "v") {
            triple vertexPosition(atof(line[1].c_str()), atof(line[2].c_str()), atof(line[3].c_str()));

            if (vertices.size() == VertexCounter) {
                // Add new vertex to fvm
                Vertex newVertex;
                newVertex.pos=vertexPosition;
                vertices.push_back(newVertex);
            }
            else {
                vertices[VertexCounter].pos = vertexPosition;
            }

            VertexCounter++;
        }

        // Vertex normal
        if(line[0] == "vn") {
            hasNormals=true;
            triple normal(atof(line[1].c_str()),atof(line[2].c_str()),atof(line[3].c_str()));

            if (vertices.size() == NormalCounter) {
                Vertex v;
                v.norm = normal;
                vertices.push_back(v);
            }
            else {
                vertices[NormalCounter].norm=normal;
            }

            NormalCounter++;
        }

        // Face
        if(line[0] == "f") {
            indexlist tmp;
            tmp.resize(nTokens-2);

            for (int i = 0; i != tmp.size(); i++) {
                tmp[i]=atoi(line[i+1].c_str())-1;
            }

            Face f1;
            f1.vert=tmp;
            faces.push_back(f1);
        }
    }

    // Close file
    fileStream.close();

    // Create additional index lists
    vlFromFaces();
    inlFromvl();

    // Evaluate mesh quality for every vertex
    for (int i=0; i!=vertices.size(); i++) {
        if( QualityAroundVertex(i) < 1 ) {
            cerr << "Bad Mesh, will not propagate properly, fix vertex " << i << endl;
        }
    }

    // If no normals are present in the input, build them
    if(!hasNormals) {
        BuildNormals();
        consistentNormals();
    }

    BuildEdges();

    return true;
}

bool FaceVertexMesh::ResetNormals()
{
     for (int i=0; i!=vertices.size(); i++) {
        triple oldn=vertices[i].norm;    
        BuildNormal(i,true,&oldn);
    }
    return true; 
}
