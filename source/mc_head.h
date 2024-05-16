#define COMPLETE_LAMBDA -1
#define NO_LIMIT -1
class MCMove
{
public:
    MCMove()
    {
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    gen.seed(rd()); //Standard mersenne_twister_engine seeded with rd()
    llim=NO_LIMIT;
    hlim=NO_LIMIT;
    prop="UNDEFINED";
    kpenal=500;
    cycle=-1;
    Beta=1.0;
    par=NULL;
    }
    virtual bool CheckEntitySupportsMove(Entity *ent);
    virtual int CreateParallelPlan(Simulation *sim);
    virtual double ExecuteMoves(Simulation *sim)=0;
    virtual bool Reset(Simulation *sim);
    int numMoves;
    double maxstepsz;
    double llim;
    double hlim;
    double Beta;

        string prop;
protected:
    double kpenal;
      std::mt19937 gen;
    set <int> busynodes;
    vector <Property *> modprop;
    deque <indexlist> modnodes; 
    void GenNewP(double *oldp, double *newp, int szprop);
    bool AcceptMove(double deltaE);
    bool TestReservation(indexlist &reservation);
    bool WriteReservation(Property *p, indexlist &reservation);
    bool BuildReservation(Property *p, indexlist &reservation);
    bool allowNormal(NTriangle *nt, triple *norm, triple *nOld);
    double edgePenalty(NTriangle *n);
    bool allowAngle(NTriangle* nagata);
    double tryFrac;
    vector <vector < Property*> > cyclemp;
    vector < deque <indexlist > > cyclemnd;
    int cycle;
    Entity *par; // One Move one Entity
};
