class Funct {
    public:
        virtual double doit(double eta, double xi) = 0;
        virtual ~Funct() = 0;
};


inline Funct::~Funct() { }


struct FaceIntegrator
{
    double Gauss7Point(Funct* fu);
    double Gauss7PointSub(Funct* fu);
    double LG7(Funct* fu);
    double G1Point(Funct* fu);
    double Z1Point(Funct* fu);

};


// Functionoids - a neccesary evil;
double FaceIntegrator::G1Point(Funct* fu)
{
    // Symmetric Gauss-Legendre 1st order.
    // Source https://pdfs.semanticscholar.org/4c92/2ee4effb71a78d9680a8646056e129d22cf1.pdf
    double rv=0;
    /*
    rv+=fu->doit(0.33333333333333,0.33333333333333)*1.000000000000;
    */

    // nagata uses triangle (0,0), (1,0), (1,1) not (0,0) (1,0), (0,1)
    rv+=fu->doit(0.66666666666667,0.33333333333333)*1.000000000000;

    return 0.5*rv;
}


double FaceIntegrator::Z1Point(Funct* fu)
{
    // Symmetric Gauss-Legendre 1st order.
    // Source https://pdfs.semanticscholar.org/4c92/2ee4effb71a78d9680a8646056e129d22cf1.pdf
    double rv=0;
    /*
    rv+=fu->doit(0.33333333333333,0.33333333333333)*1.000000000000;
    */

    // nagata uses triangle (0,0), (1,0), (1,1) not (0,0) (1,0), (0,1)
    rv+=fu->doit(0.0,0.0)*1.000000000000;

    return rv;
}


double FaceIntegrator::Gauss7Point(Funct* fu)
{
    // Symmetric Gauss-Legendre 5th order.
    // Source https://pdfs.semanticscholar.org/4c92/2ee4effb71a78d9680a8646056e129d22cf1.pdf
    double rv=0;
    /*
    rv+=fu->doit(0.33333333333333,0.33333333333333)*0.22500000000000;
    rv+=fu->doit(0.47014206410511,0.47014206410511)*0.13239415278851;
    rv+=fu->doit(0.47014206410511,0.05971587178977)*0.13239415278851;
    rv+=fu->doit(0.05971587178977,0.47014206410511)*0.13239415278851;
    rv+=fu->doit(0.10128650732346,0.10128650732346)*0.12593918054483;
    rv+=fu->doit(0.10128650732346,0.79742698535309)*0.12593918054483;
    rv+=fu->doit(0.79742698535309,0.10128650732346)*0.12593918054483;
    */

    // nagata uses triangle (0,0), (1,0), (1,1) not (0,0) (1,0), (0,1)
    rv+=fu->doit(0.66666666666667,0.33333333333333)*0.22500000000000;
    rv+=fu->doit(0.52985793589489,0.47014206410511)*0.13239415278851;
    rv+=fu->doit(0.52985793589489,0.05971587178977)*0.13239415278851;
    rv+=fu->doit(0.94028412821023,0.47014206410511)*0.13239415278851;
    rv+=fu->doit(0.89871349267654,0.10128650732346)*0.12593918054483;
    rv+=fu->doit(0.89871349267654,0.79742698535309)*0.12593918054483;
    rv+=fu->doit(0.20257301464691,0.10128650732346)*0.12593918054483;

    return 0.5*rv;
}


double FaceIntegrator::Gauss7PointSub(Funct* fu)
{
    // Symmetric Gauss-Legendre 5th order.
    // Source https://pdfs.semanticscholar.org/4c92/2ee4effb71a78d9680a8646056e129d22cf1.pdf
    double rv=0;

    rv += fu->doit(0.33333333333333, 0.16666666666667) * 0.22500000000000;
    rv += fu->doit(0.26492896794745, 0.23507103205256) * 0.13239415278851;
    rv += fu->doit(0.26492896794745, 0.02985793589489) * 0.13239415278851;
    rv += fu->doit(0.47014206410512, 0.23507103205256) * 0.13239415278851;
    rv += fu->doit(0.44935674633827, 0.05064325366173) * 0.12593918054483;
    rv += fu->doit(0.44935674633827, 0.39871349267655) * 0.12593918054483;
    rv += fu->doit(0.10128650732346, 0.05064325366173) * 0.12593918054483;

    rv += fu->doit(0.83333333333333, 0.16666666666667) * 0.22500000000000;
    rv += fu->doit(0.76492896794745, 0.23507103205256) * 0.13239415278851;
    rv += fu->doit(0.76492896794745, 0.02985793589489) * 0.13239415278851;
    rv += fu->doit(0.97014206410512, 0.23507103205256) * 0.13239415278851;
    rv += fu->doit(0.94935674633827, 0.05064325366173) * 0.12593918054483;
    rv += fu->doit(0.94935674633827, 0.39871349267655) * 0.12593918054483;
    rv += fu->doit(0.60128650732346, 0.05064325366173) * 0.12593918054483;

    rv += fu->doit(0.83333333333333, 0.66666666666667) * 0.22500000000000;
    rv += fu->doit(0.76492896794745, 0.73507103205256) * 0.13239415278851;
    rv += fu->doit(0.76492896794745, 0.52985793589489) * 0.13239415278851;
    rv += fu->doit(0.97014206410512, 0.73507103205256) * 0.13239415278851;
    rv += fu->doit(0.94935674633827, 0.55064325366173) * 0.12593918054483;
    rv += fu->doit(0.94935674633827, 0.89871349267655) * 0.12593918054483;
    rv += fu->doit(0.60128650732346, 0.55064325366173) * 0.12593918054483;

    rv += fu->doit(0.66666666666667, 0.33333333333333) * 0.22500000000000;
    rv += fu->doit(0.73507103205255, 0.26492896794744) * 0.13239415278851;
    rv += fu->doit(0.73507103205255, 0.47014206410511) * 0.13239415278851;
    rv += fu->doit(0.52985793589488, 0.26492896794744) * 0.13239415278851;
    rv += fu->doit(0.55064325366173, 0.44935674633827) * 0.12593918054483;
    rv += fu->doit(0.55064325366173, 0.10128650732345) * 0.12593918054483;
    rv += fu->doit(0.89871349267654, 0.44935674633827) * 0.12593918054483;

    return 0.125*rv;
}


double FaceIntegrator::LG7(Funct*fu)
{
    double rv=0; // M.E. Laursen, M. Gellert
                 // International Journal for Numerical Methods in Engineering, vol. 12, no. 1, pp. 67–76, 1978
                 // https://doi.org/10.1002/nme.1620120107
    rv+=fu->doit(0.935069486840835,0.870138973681670)*0.053077801790232;
    rv+=fu->doit(0.935069486840835,0.064930513159165)*0.053077801790232;
    rv+=fu->doit(0.129861026318330,0.064930513159165)*0.053077801790232;
    rv+=fu->doit(0.715424415750830,0.198384476681507)*0.070853083692134;
    rv+=fu->doit(0.686440815615069,0.642577343822696)*0.069274682079417;
    rv+=fu->doit(0.801615523318493,0.517039939069323)*0.070853083692134;
    rv+=fu->doit(0.357422656177304,0.043863471792372)*0.069274682079417;
    rv+=fu->doit(0.482960060930677,0.284575584249170)*0.070853083692134;
    rv+=fu->doit(0.956136528207628,0.313559184384932)*0.069274682079417;
    rv+=fu->doit(0.482960060930677,0.198384476681507)*0.070853083692134;
    rv+=fu->doit(0.956136528207628,0.642577343822696)*0.069274682079417;
    rv+=fu->doit(0.801615523318493,0.284575584249170)*0.070853083692134;
    rv+=fu->doit(0.357422656177304,0.313559184384932)*0.069274682079417;
    rv+=fu->doit(0.715424415750830,0.517039939069323)*0.070853083692134;
    rv+=fu->doit(0.686440815615069,0.043863471792372)*0.069274682079417;

    return 0.5*rv;
}


class FaceProperty
{
    public:
        // Constructors
        FaceProperty(){}
        FaceProperty(string& pnam){
            name=pnam;
        }

        // Methods
        inline void setProperty(double prop) {ppt=prop;}
        inline double getProperty(void) {return ppt;}
        string getName(void) {return name;}

    private:
        double ppt;
        string name;
};


class PropertiedFace
{
    public:
        vector<triple> c; // Coeffs Linked To Edges
        vector<triple*> v; // Vertices Linked to Mesh
        vector<triple*> n; // Normals;
        // Methods
        bool AddProperty(FaceProperty &fp) {
            props[fp.getName()]=fp.getProperty();
            return true;
        }

        double GetProperty(const string &s){
            try{
            return props.at(s);
            }
            catch (std::exception &e)
            {
                return 0;
            }
            
        }
        
        map<string, double >* GetAll(void)
        {return &props;}
        
        bool ExportProperties(vector< pair <string, double > > *vds)
        {
        vds->clear();
        std::map<std::string, double >::iterator it = props.begin();
        while (it!=props.end())
        {
            vds->push_back(make_pair(it->first,it->second));
            
            it++;
        }
            
        return true;
        }
        bool SetProperty(const string &s, double prop)
        {
            props[s]=prop;
            return true;
        }

        virtual triple eval(double eta,double xi)=0;
        virtual triple x_eta(double eta, double xi)=0;
        virtual triple x_xi(double eta, double xi)=0;
        virtual triple normal(double eta, double xi)=0;
        virtual triple eval_p(double eta,double xi, int p)=0;
        virtual triple x_eta_p(double eta, double xi, int p)=0;
        virtual triple x_xi_p(double eta, double xi, int p)=0;
        virtual triple normal_p(double eta, double xi, int p)=0;
        virtual triple cross_p(double eta, double xi, int p)=0;
        virtual triple c_p(int cnum,int p)=0;
        virtual bool updatec()=0;
        virtual void giveC(double eta, double xi, vector<triple> *cOut) = 0;
        virtual bool prepgrad()=0;

    private:
        map<string, double> props;

//    virtual double Area() = 0;
//    virtual double HInt() =0;
//    virtual double H2Int() =0;
//    virtual double KInt() =0;
//    virtual double VIncr() =0;
//    indexlist exposedv;
//    indexlist exposedn;
//    int faceID;
//    virtual double Energy() =0;
//  // etc.
};
