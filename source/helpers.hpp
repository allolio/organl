class SimulationProperty
{
    public:
        // constructors
        SimulationProperty(){}
        SimulationProperty(string pname)
        {
            name = pname;
        }

        // methods
        inline void setProperty(double prop)
        {
            value = prop;
        }

        inline double getProperty(void)
        {
            return value;
        }

        string getName(void)
        {
            return name;
        }

    private:
        string name;
        double value;
};

class SimulationSetup
{
    public:

        bool ReadProperties(string file);

        bool AddProperty(SimulationProperty &sp)
        {
            props.insert(make_pair(sp.getName(), sp));
            return true;
        }

        int getNumberOfMonteCarloSteps()
        {
            if (props.find("nMCSteps") != props.end())
                return (int)(props["nMCSteps"]).getProperty();
            return 20000;
        }

        int getOutputFrequency()
        {
            if (props.find("writeEvery") != props.end())
                return (int)(props["writeEvery"]).getProperty();
            return 1000;
        }

        int getRefreshNeighborsFrequency()
        {
            if (props.find("updateNghbrsEvery") != props.end())
                return (int)(props["updateNghbrsEvery"]).getProperty();
            return 200;
        }

        double getNeighboroodRadius()
        {
            if (props.find("nghbRadius") != props.end())
                return (props["nghbRadius"]).getProperty();
            return -1.0;
        }

    private:

        unordered_map<string, SimulationProperty> props;
};


bool SimulationSetup::ReadProperties(std::string file)
{
    // Format: * name value
    // lines that not begin with * are ignored
    ifstream hfile;
    hfile.open(file.c_str());
    if(hfile.fail()) return false;
    string current;
    while ( !hfile.eof() )
    {
        vector<string> line;
        getline(hfile, current);
        int htokens=Tokenize(current,&line, " \t");
        if(htokens<3)
            continue;
        //cout << line[1] << " "<< line[2] << endl;

        if(line[0]=="*")
        {
            SimulationProperty sprop(line[1]);
            sprop.setProperty(atof(line[2].c_str()));
            AddProperty(sprop);
        }
    }

    hfile.close();
    return true;
}
