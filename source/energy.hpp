//

#define MEMBP_LAZY -1
#define LAMBDA_DEFAULT -1

struct FrozenCoeff
{
    int vertex;
    bool isnormal;
    int dim; // 0,1,2 x,y,z; or phi/theta respectively;
};

struct Coefficient
{
    double *val;
    int vertex; // Bundle
    int type;
};



/**
 *
 * the principal curvatures of sphere are 1/R^2
 * the Helfrich Energy of sphere is 4\pi\kappa
 *
 */





class Energy : public GeoControl
{
    public:
        bool init;

        Energy () {E0=INVALID_BULK;init=false; lambdas.clear(); defaults.clear();}

        double E0;
        double E;
        vector <double> lambdas;
        vector <double> defaults;
        virtual void initializeValues(PropertyMap &tpm) = 0;
        
        virtual double eval (PropertiedFace *f)=0;
        
        virtual double deep_eval(PropertiedFace *f, bool trye=false)=0;
        
        virtual void deep_props(vector<string> *obs, vector<int> *ind) = 0;

        virtual double evalGradient (PropertiedFace *f, int p)=0;

        virtual double computeTotalEnergy()=0;

 //       virtual vector<Coefficient> giveGradSpecificCoeffs(void){vector<Coefficient> mon;return mon;}

};


class HelfrichEnergy : public Energy
{
    public:
        void initializeValues(PropertyMap &tpm)
        {
        }

        double eval(PropertiedFace *f)
        {
            HelfE e(f);
            double kappa= f->GetProperty("kappa");
            return fi.Gauss7PointSub(&e)*kappa/2.0;
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
           double energy=eval(f); 
           if(trye)
           f->SetProperty("Try_Energy",energy);
           else 
           f->SetProperty("Energy",energy);
           return energy;
        }
        
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
            if(!obs->empty()) obs->clear();
            if(!ind->empty()) ind->clear();
        }

        double evalGradient(PropertiedFace *f, int p)
        {
            dInt h(f,p,f->GetProperty("kappa"));
            return fi.Gauss7PointSub(&h);
        };  
        
      
        class HelfE : 
            public Funct
            {
                public: 
                    HelfE(PropertiedFace *nga)
                {
                    ngt = nga;
                    c0 = ngt->GetProperty("c0");

                }

                double doit(double eta,double xi)
                {
//                      cout << ngt->GetProperty("kappa");
 
                    double H=meanCurvature(eta, xi, *ngt)*2.;
                    H -= c0;
                    return  H*H * (ngt->x_eta(eta, xi) % ngt->x_xi(eta, xi)).abs();//ngt->GetProperty("kappa")*H*H* (ngt->x_xi(eta, xi) % ngt->x_xi(eta, xi)).abs();
                }
                double c0;
                PropertiedFace *ngt;
            };

        class dInt : public Funct
        {
        public:
            dInt(PropertiedFace *nga,int p,double k)
            {
                ngt = nga;
                param = p;
                kappa = k;
            }

            double doit(double eta,double xi)
            {
                double tol = 1e-14;
                triple n   = ngt->normal(eta, xi);
                triple xe  = ngt->x_eta(eta, xi);
                triple xx  = ngt->x_xi(eta, xi);
                triple dxe = ngt->x_eta_p(eta,xi,param);
                triple dxx = ngt->x_xi_p(eta,xi,param);
                triple dn  = ngt->normal_p(eta,xi,param);
                triple crss = xe%xx;
                double xexe = xe*xe;
                double dxexe = 2*(dxe*xe);
                double xxxx = xx*xx;
                double dxxxx= 2*(dxx*xx);
                double xexx = xe*xx;
                double dxexx = (dxe*xx)+(xe*dxx);
                double nc0 = n*(ngt->c[0]);
                double dnc0 = dn*(ngt->c[0])+n*ngt->c_p(0,param);
                double nc1 = n*(ngt->c[1]);
                double dnc1 = dn*(ngt->c[1])+n*ngt->c_p(1,param);
                double detI  = xexe * xxxx - xexx * xexx;
                double ddetI = dxexe*xxxx + xexe*dxxxx - dxexx*xexx - xexx*dxexx;

                // det of 1st fundamental form
                // ddetI = dxexe*xxxx+xexe+dxxxx-dxexxx*xexx-xexx-dxexx;   
                if (detI < tol and detI > - tol)
                {
                    cout << "Unable to calculate Mean curvature, determinant of the 1st fundamental form is zero!";
                    return 0;
                }

                triple innerc  = ngt->c[2] - ngt->c[1] - ngt->c[0];
                triple dinnerc = ngt->c_p(2,param) - ngt->c_p(1,param) - ngt->c_p(0,param);
                double dinner  = dn*innerc + n*dinnerc;
                double H   = 2*(nc0)*xxxx + 2*(nc1)*xexe - 2*(n*(innerc))*xexx;
                double Hm  = H / (2*detI);
                double dH  = 2*dnc0*xxxx + 2*nc0*dxxxx + 2*dnc1*xexe + 2*nc1*dxexe - 2*( dinner*xexx + (n*innerc)*dxexx );
                double u   = 2*detI;
                double dHm = dH/(u) - H/(u*u)*(2*ddetI);

                double dS = crss.abs();
                double ddS = 1/crss.abs()*(ngt->cross_p(eta,xi,param)*crss);

                return  kappa*((2*Hm)*dHm*dS + Hm*Hm*ddS);
            }
            PropertiedFace *ngt;
            int param;
            double kappa;
        };
};


class ConstraintHelfrichEnergy : public HelfrichEnergy
{
    public:
        ConstraintHelfrichEnergy () {E0=INVALID_BULK; 
            init=false;lambdas.resize(2); defaults.resize(2);
            defaults[0]=INVALID_BULK;defaults[1]=INVALID_BULK; 
            lambdas[0]=1e-10; lambdas[1]=1e-18;
            
        }
        double *V;
        double *A;
   
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            {
            Property &p=tpm.GetProperty("MembG",0);
            if(p.state==MEMBP_LAZY) p.update();
            A=p.content;
            cout << "Volume ?" << endl;
            Property &l=tpm.GetProperty("MembG",1);
            if(l.state==MEMBP_LAZY) l.update();
            V=l.content;
            if(defaults[1]==INVALID_BULK) defaults[1]=*V;
            if(defaults[0]==INVALID_BULK) defaults[0]=*A;
            }
            init=true;
        }

        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
    
            double kappa= f->GetProperty("kappa");
            
            double Hint = fi.Gauss7PointSub(&e);
            //double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);

            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa/2.0;
        }
        
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
            if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
            obs->push_back("Area");
            obs->push_back("Volume");
            ind->push_back(0);
            ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
            double Volu=VIncr(f); // TODO Put all increments into one
            double Are=Area(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
           f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
 
           f->SetProperty("Try_Area",Are);
           }
            else
            {
           f->SetProperty("Energy",energy);
           f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
           f->SetProperty("Area",Are);
            }
          return energy;
        }
        
        double evalGradient(PropertiedFace *f, int p)
        {
            dInt h(f,p,f->GetProperty("kappa") , lambdas[1], lambdas[0], *V, defaults[1], *A, defaults[0]);
            return fi.Gauss7PointSub(&h);
        };

        double computeTotalEnergy()
        {
            return 0.5*lambdas[1]*(*V - defaults[1])*(*V - defaults[1]) + 0.5*lambdas[0]*(*A - defaults[0])*(*A - defaults[0]);
        }
      
        class dInt : public Funct
        {
        public:
            dInt(PropertiedFace *nga, int p, double k, double l1, double l2,
                    double Vol, double Vol0, double Ar, double Ar0)
            {
                ngt = nga;
                param = p;
                kappa = k;
                lVol = l1;
                lSur = l2;
                volume = Vol;
                initialVolume = Vol0;
                area = Ar;
                initialArea = Ar0;
            }

            double doit(double eta,double xi)
            {
                double tol = 1e-14;
                triple n   = ngt->normal(eta, xi);
                triple xe  = ngt->x_eta(eta, xi);
                triple xx  = ngt->x_xi(eta, xi);
                triple dxe = ngt->x_eta_p(eta,xi,param);
                triple dxx = ngt->x_xi_p(eta,xi,param);
                triple dn  = ngt->normal_p(eta,xi,param);
                triple crss = xe%xx;
                double xexe = xe*xe;
                double dxexe = 2*(dxe*xe);
                double xxxx = xx*xx;
                double dxxxx= 2*(dxx*xx);
                double xexx = xe*xx;
                double dxexx = (dxe*xx)+(xe*dxx);
                double nc0 = n*(ngt->c[0]);
                double dnc0 = dn*(ngt->c[0])+n*ngt->c_p(0,param);
                double nc1 = n*(ngt->c[1]);
                double dnc1 = dn*(ngt->c[1])+n*ngt->c_p(1,param);
                double detI  = xexe * xxxx - xexx * xexx;
                double ddetI = dxexe*xxxx + xexe*dxxxx - dxexx*xexx - xexx*dxexx;

                // det of 1st fundamental form
                // ddetI = dxexe*xxxx+xexe+dxxxx-dxexxx*xexx-xexx-dxexx;   
                if (detI < tol and detI > - tol)
                {
                    cerr << "Unable to calculate Mean curvature, determinant of the 1st fundamental form is zero!";
                    return 0;
                }

                triple innerc  = ngt->c[2] - ngt->c[1] - ngt->c[0];
                triple dinnerc = ngt->c_p(2,param) - ngt->c_p(1,param) - ngt->c_p(0,param);
                double dinner  = dn*innerc + n*dinnerc;
                double H   = 2*(nc0)*xxxx + 2*(nc1)*xexe - 2*(n*(innerc))*xexx;
                double Hm  = H / (2*detI);
                double dH  = 2*dnc0*xxxx + 2*nc0*dxxxx + 2*dnc1*xexe + 2*nc1*dxexe - 2*( dinner*xexx + (n*innerc)*dxexx );
                double u   = 2*detI;
                double dHm = dH/(u) - H/(u*u)*(2*ddetI);

                double dS = crss.abs();
                double ddS = 1/crss.abs()*(ngt->cross_p(eta,xi,param)*crss);

                double dHelfrich = kappa*((2*Hm)*dHm*dS + Hm*Hm*ddS);

                double dVolume = 2.0/3.0*(volume - initialVolume)*(ngt->eval_p(eta, xi, param)*n
                        + ngt->eval(eta,xi)*dn)*dS
                    + 1.0/3.0*(volume - initialVolume)*(volume - initialVolume)*(ngt->eval(eta,xi)*n)*ddS;

                double dSurface = 2*(area - initialArea)*dS + (area - initialArea)*(area - initialArea)*ddS;

                double result = dHelfrich + 0.5*lVol*dVolume + 0.5*lSur*dSurface;

                return result;
            }

            PropertiedFace *ngt;
            int param;
            double kappa;
            double lVol;
            double lSur;
            double volume;
            double initialVolume;
            double area;
            double initialArea;
        };

};
        class HelfrichElasticAreaPressure : public HelfrichEnergy
{
     public:
        HelfrichElasticAreaPressure  () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(1); defaults.resize(0);
            //defaults[0]=INVALID_BULK; 
            lambdas[0]=1;
            
        }
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            {
            init=true;}
        }
        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
            double Are=Area(f);
            double V=VIncr(f);

            //dA a(f);
            //dV v(f);
   
            double kappa= f->GetProperty("kappa");
            double A0=f->GetProperty("A0");
            double Ka=f->GetProperty("Ka");
            if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
           // cout << "A0: " << A0 << "Area " << Are << "Elastic E " << 0.5*Ka*(A0-Are)*(A0-Are) << endl;
            double Hint = fi.Gauss7PointSub(&e);
            //double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);

            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa/2.0+0.5*Ka*(A0-Are)*(A0-Are)-lambdas[0]*V;
        }
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
                if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
  //          obs->push_back("Area");
       //     obs->push_back("Volume");
//            ind->push_back(0);
         //   ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
         //   double Volu=VIncr(f); // TODO Put all increments into one
        //    double Are=Area(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
         //  f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
            }
            else
            {
           f->SetProperty("Energy",energy);
          // f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
            }
          return energy;
        }
        
   double evalGradient(PropertiedFace *f, int p)
        {
            cerr << "Gradient not implemented yet" << endl;
            return 0.0;
        };

        double computeTotalEnergy()
        {
            return 0;
        }
}; 

            class ElasticAreaHelfrichEnergy : public HelfrichEnergy
{
     public:
        ElasticAreaHelfrichEnergy () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(0); defaults.resize(0);
           // defaults[0]=INVALID_BULK; 
           // lambdas[0]=1e-18;
            
        }
      //  double *V;
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            {
           // cout << "Volume ?" << endl;
           // Property &l=tpm.GetProperty("MembG",1);
            //if(l.state==MEMBP_LAZY) l.update();
            //V=l.content;
            //if(defaults[0]==INVALID_BULK) defaults[0]=*V;
                         init=true;   
            }
        }
        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
            double Are=Area(f);
            //dA a(f);
            //dV v(f);
   
            double kappa= f->GetProperty("kappa");
            double A0=f->GetProperty("A0");
            double Ka=f->GetProperty("Ka");
            if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
           // cout << "A0: " << A0 << "Area " << Are << "Elastic E " << 0.5*Ka*(A0-Are)*(A0-Are) << endl;
            double Hint = fi.Gauss7PointSub(&e);
            //double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);

            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa*0.5+0.5*Ka*(A0-Are)*(A0-Are);
        }
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
                if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
  //          obs->push_back("Area");
          //  obs->push_back("Volume");
//            ind->push_back(0);
        //    ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
           //f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
            }
            else
            {
           f->SetProperty("Energy",energy);
          // f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
            }
          return energy;
        }
        
   double evalGradient(PropertiedFace *f, int p)
        {
            cerr << "Gradient not implemented yet" << endl;
            return 0.0;
        };

        double computeTotalEnergy()
        {
            return 0.0; //0.5*lambdas[0]*(*V - defaults[0])*(*V - defaults[0]) ;
        }
         
        

      /*  vector<Coefficient> giveGradSpecificCoeffs(void)
        {
            vector<Coefficient> coeffs;
            Coefficient lagVol;
            lagVol.vertex = -1;
            lagVol.type = 6;
            lagVol.val = &lambdaVol;
            Coefficient lagSur;
            lagSur.vertex = -1;
            lagSur.type = 7;
            lagSur.val = &lambdaSur;

            coeffs.clear();
            coeffs.push_back(lagVol);
            coeffs.push_back(lagSur);
        unordered_map<string, EnergyProperty> props;
        double Ecur;
        double Vcur;
        double Acur;
        
        
        bool Energy::ReadProperties(std::string file)
{
    // Format: * name value
    // lines that not begin with * are ignored
    ifstream hfile;
    hfile.open(file.c_str());
    if(hfile.fail()) return false;
    string current;
    while ( getline(hfile, current))
    {
        vector<string> line;
        int htokens=Tokenize(current,&line, " \t");
        if(htokens<3)
            continue;
        //cout << line[1] << " "<< line[2] << endl;

        if(line[0]=="*")
        {
            EnergyProperty gprop(line[1]);
            gprop.setProperty(atof(line[2].c_str()));
            AddProperty(gprop);
        }
    }

    hfile.close();
    return true;
}
            return coeffs;
            
                    bool AddProperty(EnergyProperty &ep)
        {
            props.insert(make_pair(ep.getName(), ep));
            return true;
        }

        double GetProperty(const string &s)
        {
            return (props[s]).getProperty();
        }

        bool HaveProperty(const string & s)
        {
            return props.find(s) != props.end();
        }

        bool ReadProperties(string file);

        }*/
    };

  
              class TensionPressureHelfrichEnergy : public HelfrichEnergy
{
     public:
        TensionPressureHelfrichEnergy () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(2); defaults.resize(0);
           // defaults[0]=INVALID_BULK; 
            lambdas[0]=1e-18;
            lambdas[1]=1e-18;
        }
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            
         //   cout << "Volume ?" << endl;
         //   Property &l=tpm.GetProperty("MembG",1);
         //   if(l.state==MEMBP_LAZY) l.update();
          //  V=l.content;
            init=true;
        }
        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
            double Are=Area(f);
            //dA a(f);
            double V=VIncr(f);
   
            double kappa= f->GetProperty("kappa");
         //   if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
           // cout << "A0: " << A0 << "Area " << Are << "Elastic E " << 0.5*Ka*(A0-Are)*(A0-Are) << endl;
            double Hint = fi.Gauss7PointSub(&e);
          //  double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);
            f->SetProperty("Try_Area",Are); // For correct penalty.
            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa*0.5+lambdas[1]*Are-lambdas[0]*V;
        }
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
                if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
  //          obs->push_back("Area");
  //          obs->push_back("Volume");
//            ind->push_back(0);
    //        ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
  //          double Volu=VIncr(f); // TODO Put all increments into one
        //    double Are=Area(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
         //  f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
            }
            else
            {
           f->SetProperty("Energy",energy);
          // f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
            }
          return energy;
  
        }
        
   double evalGradient(PropertiedFace *f, int p)
        {
            cerr << "Gradient not implemented yet" << endl;
            return 0.0;
        };

        double computeTotalEnergy()
        {
            return 0; //-lambdas[0]*(*V);
        }
};

    
            class AreaPressureHelfrichEnergy : public HelfrichEnergy
{
     public:
        AreaPressureHelfrichEnergy () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(1); defaults.resize(0);
           // defaults[0]=INVALID_BULK; 
            lambdas[0]=1e-18;
            
        }
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            
         //   cout << "Volume ?" << endl;
         //   Property &l=tpm.GetProperty("MembG",1);
         //   if(l.state==MEMBP_LAZY) l.update();
          //  V=l.content;
            init=true;
        }
        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
            double Are=Area(f);
            //dA a(f);
            double V=VIncr(f);
   
            double kappa= f->GetProperty("kappa");
             double A0=f->GetProperty("A0");
            double Ka=f->GetProperty("Ka");
         //   if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
           // cout << "A0: " << A0 << "Area " << Are << "Elastic E " << 0.5*Ka*(A0-Are)*(A0-Are) << endl;
            double Hint = fi.Gauss7PointSub(&e);
          //  double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);

            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa*0.5+Ka*abs(A0-Are)-lambdas[0]*V;
        }
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
                if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
  //          obs->push_back("Area");
  //          obs->push_back("Volume");
//            ind->push_back(0);
    //        ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
  //          double Volu=VIncr(f); // TODO Put all increments into one
        //    double Are=Area(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
         //  f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
            }
            else
            {
           f->SetProperty("Energy",energy);
          // f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
            }
          return energy;
  
        }
        
   double evalGradient(PropertiedFace *f, int p)
        {
            cerr << "Gradient not implemented yet" << endl;
            return 0.0;
        };

        double computeTotalEnergy()
        {
            return 0; //-lambdas[0]*(*V);
        }
};


class DCADEHelfrichEnergy : public HelfrichEnergy
{
    public:
        DCADEHelfrichEnergy () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(3); defaults.resize(3);
            defaults[0]=INVALID_BULK;defaults[1]=INVALID_BULK;defaults[2]=INVALID_BULK; 
            lambdas[0]=1e-10; lambdas[1]=1e-18; lambdas[2]=1e-18;
            
        }
        double *V;
        double *A;
        double *M;
//        double V0;
 //       double A0;
     
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            {
            Property &p=tpm.GetProperty("MembG",0);
            if(p.state==MEMBP_LAZY) p.update();
            A=p.content;
            cout << "Volume ?" << endl;
            Property &l=tpm.GetProperty("MembG",1);
            if(l.state==MEMBP_LAZY) l.update();
            Property &m=tpm.GetProperty("MembG",2);
            if(m.state==MEMBP_LAZY) m.update();
            M=m.content;
            V=l.content;
            if(defaults[1]==INVALID_BULK) defaults[1]=*V;
            if(defaults[0]==INVALID_BULK) defaults[0]=*A;
            if(defaults[2]==INVALID_BULK) defaults[2]=*M;
            }
            init=true;
        }

        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
            //dA a(f);
            //dV v(f);
            double kappa= f->GetProperty("kappa");
            double Are=Area(f);
            double Hint = fi.LG7(&e);
            //double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);
            double A0=f->GetProperty("A0");
            double Ka=f->GetProperty("Ka");
            if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa*0.5+Ka*(A0-Are)*(A0-Are);
        }
        
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
            if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
            obs->push_back("Area");
            obs->push_back("Volume");
            obs->push_back("Curvature");
            ind->push_back(0);
            ind->push_back(1);
            ind->push_back(2);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
            double Volu=VIncr(f); // TODO Put all increments into one
            double Are=Area(f);
            double Curv=HInt(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
           f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
           f->SetProperty("Try_Curvature",Curv);
           f->SetProperty("Try_Area",Are);
           }
            else
            {
           f->SetProperty("Energy",energy);
           f->SetProperty("Volume",Volu);
           f->SetProperty("Curvature",Curv);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
           f->SetProperty("Area",Are);
            }
          return energy;
        }
        
        double evalGradient(PropertiedFace *f, int p)
        {
          //  dInt h(f,p,f->GetProperty("kappa") , lambdas[1], lambdas[0], *V, defaults[1], *A, defaults[0]);
            return 0; //fi.Gauss7PointSub(&h);
        };

        double computeTotalEnergy()
        {
            return 0.5*lambdas[1]*(*V - defaults[1])*(*V - defaults[1]) + 0.5*lambdas[0]*(*A - defaults[0])*(*A - defaults[0])+0.5*lambdas[2]*(*M - defaults[2])*(*M - defaults[2]);
        }
      
     
            PropertiedFace *ngt;
            int param;
            double kappa;
            double lVol;
            double lSur;
            double volume;
            double initialVolume;
            double area;
            double initialArea;
        };
 

class DoubleConstraintHelfrichEnergy : public HelfrichEnergy
{
    public:
        DoubleConstraintHelfrichEnergy () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(2); defaults.resize(2);
            defaults[0]=INVALID_BULK;defaults[1]=INVALID_BULK; 
            lambdas[0]=1e-10; lambdas[1]=1e-18;
            
        }
        double *V;
        double *A;
//        double V0;
 //       double A0;
     
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            {
            Property &p=tpm.GetProperty("MembG",0);
            if(p.state==MEMBP_LAZY) p.update();
            A=p.content;
            cout << "Volume ?" << endl;
            Property &l=tpm.GetProperty("MembG",1);
            if(l.state==MEMBP_LAZY) l.update();
            V=l.content;
            if(defaults[1]==INVALID_BULK) defaults[1]=*V;
            if(defaults[0]==INVALID_BULK) defaults[0]=*A;
            }
            init=true;
        }

        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
            //dA a(f);
            //dV v(f);
   
            double kappa= f->GetProperty("kappa");
            double Are=Area(f);
            double Hint = fi.LG7(&e);
            //double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);
            double A0=f->GetProperty("A0");
            double Ka=f->GetProperty("Ka");
            if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa*0.5+Ka*(A0-Are)*(A0-Are);
        }
        
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
            if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
            obs->push_back("Area");
            obs->push_back("Volume");
            ind->push_back(0);
            ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
            double Volu=VIncr(f); // TODO Put all increments into one
            double Are=Area(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
           f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
 
           f->SetProperty("Try_Area",Are);
           }
            else
            {
           f->SetProperty("Energy",energy);
           f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
           f->SetProperty("Area",Are);
            }
          return energy;
        }
        
        double evalGradient(PropertiedFace *f, int p)
        {
          //  dInt h(f,p,f->GetProperty("kappa") , lambdas[1], lambdas[0], *V, defaults[1], *A, defaults[0]);
            return 0; //fi.Gauss7PointSub(&h);
        };

        double computeTotalEnergy()
        {
            return 0.5*lambdas[1]*(*V - defaults[1])*(*V - defaults[1]) + 0.5*lambdas[0]*(*A - defaults[0])*(*A - defaults[0]);
        }
      
     
            PropertiedFace *ngt;
            int param;
            double kappa;
            double lVol;
            double lSur;
            double volume;
            double initialVolume;
            double area;
            double initialArea;
        };
        
class DoubleConstraintGHelfrichEnergy : public HelfrichEnergy
{
    public:
        DoubleConstraintGHelfrichEnergy () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(2); defaults.resize(2);
            defaults[0]=INVALID_BULK;defaults[1]=INVALID_BULK; 
            lambdas[0]=1e-10; lambdas[1]=1e-18;
            
        }
        double *V;
        double *A;
//        double V0;
 //       double A0;
     class GHelfE : 
            public Funct
            {
                public: 
                    GHelfE(PropertiedFace *nga)
                {
                    ngt = nga;
                    c0 = ngt->GetProperty("c0");
                    kappa=ngt->GetProperty("kappa");
                    kg=ngt->GetProperty("kappag");
                }

                double doit(double eta,double xi)
                {
//                      cout << ngt->GetProperty("kappa");
  vector <triple> c(3);
        double tol = 1e-12;
        ngt->giveC(eta, xi, &c);
        triple n   = ngt->normal(eta, xi);
        triple xe  = ngt->x_eta(eta, xi);
        triple xx  = ngt->x_xi(eta, xi);
        double nc0 = n*c[0];
        double nc1 ;//= n*c[1];
        double xxxx=xx*xx;
        double xexe=xe*xe;
        double xexx=xe*xx;
        double detI  = (xexe) * (xxxx) - (xexx) * (xexx);                                       // det of 1st fundamental form
         double detII;
                double cross;
        if(c.size()==10)
        {
            cross = (c[6]+c[8]*2*eta+c[9]*2*xi)*n;
            nc0=nc0+(c[8]*xi+c[3]*3*eta)*n; //
            nc1=(c[7]+c[9]*eta+c[4]*3*xi)*n;
            detII=(4*nc0*nc1 - cross*cross);
        } else 
        {
               double nc2 = n*c[2];
            nc1 = n*c[1];
            cross=n*(c[2] - c[1] - c[0]); //2nd fundamental cross term
            detII = 2*(nc2)*(nc0 + nc1) + 2*nc0*nc1 - nc0*nc0 - nc1*nc1 - nc2*nc2;  // det of 2nd fundamental form
        }

        if (detI < tol and detI > - tol)
        {
            cout << "Unable to calculate Gauss/Mean curvature, determinant of the 1st fundamental form is zero!";
            return 0;
        }

        double H = 2*(nc0)*xxxx + 2*(nc1)*xexe - 2*(cross)*xexx;
        H/=(detI);
        double K=detII/detI;
        H -= c0;
                    return  (kappa*H*H*0.5+kg*K) * (ngt->x_eta(eta, xi) % ngt->x_xi(eta, xi)).abs();//ngt->GetProperty("kappa")*H*H* (ngt->x_xi(eta, xi) % ngt->x_xi(eta, xi)).abs();
                }
                double kappa;
                double kg;
                double c0;
                PropertiedFace *ngt;
            };

        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            {
            Property &p=tpm.GetProperty("MembG",0);
            if(p.state==MEMBP_LAZY) p.update();
            A=p.content;
            cout << "Volume ?" << endl;
            Property &l=tpm.GetProperty("MembG",1);
            if(l.state==MEMBP_LAZY) l.update();
            V=l.content;
            if(defaults[1]==INVALID_BULK) defaults[1]=*V;
            if(defaults[0]==INVALID_BULK) defaults[0]=*A;
            }
            init=true;
        }

        double eval(PropertiedFace *f)
        {
            GHelfE e(f);    
            //dA a(f);
            //dV v(f);
   
            //double kappa= f->GetProperty("kappa");
            double Are=Area(f);
            double Hint = fi.LG7(&e);
            //double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);
            double A0=f->GetProperty("A0");
            double Ka=f->GetProperty("Ka");
            if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint+Ka*(A0-Are)*(A0-Are);
        }
        
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
            if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
            obs->push_back("Area");
            obs->push_back("Volume");
            ind->push_back(0);
            ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
            double Volu=VIncr(f); // TODO Put all increments into one
            double Are=Area(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
           f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
 
           f->SetProperty("Try_Area",Are);
           }
            else
            {
           f->SetProperty("Energy",energy);
           f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
           f->SetProperty("Area",Are);
            }
          return energy;
        }
        
        double evalGradient(PropertiedFace *f, int p)
        {
           // dInt h(f,p,f->GetProperty("kappa") , lambdas[1], lambdas[0], *V, defaults[1], *A, defaults[0]);
            return 0;// fi.Gauss7PointSub(&h);
        };

        double computeTotalEnergy()
        {
            return 0.5*lambdas[1]*(*V - defaults[1])*(*V - defaults[1]) + 0.5*lambdas[0]*(*A - defaults[0])*(*A - defaults[0]);
        }
      
            PropertiedFace *ngt;
            int param;
            double kappa;
            double lVol;
            double lSur;
            double volume;
            double initialVolume;
            double area;
            double initialArea;
//        };


};
class DoubleConstraintLHelfrichEnergy : public HelfrichEnergy
{
    public:
        DoubleConstraintLHelfrichEnergy () {E0=INVALID_BULK; 
            /*lambdaVol=1e-18, lambdaSur=1e-10;*/ init=false;lambdas.resize(2); defaults.resize(2);
            defaults[0]=INVALID_BULK;defaults[1]=INVALID_BULK; 
            lambdas[0]=1e-10; lambdas[1]=1e-18;
            
        }
        double *V;
        double *A;
//        double V0;
 //       double A0;
     
        void initializeValues(PropertyMap &tpm)
        {
            if(!init)
            {
            Property &p=tpm.GetProperty("MembG",0);
            if(p.state==MEMBP_LAZY) p.update();
            A=p.content;
            cout << "Volume ?" << endl;
            Property &l=tpm.GetProperty("MembG",1);
            if(l.state==MEMBP_LAZY) l.update();
            V=l.content;
            if(defaults[1]==INVALID_BULK) defaults[1]=*V;
            if(defaults[0]==INVALID_BULK) defaults[0]=*A;
            }
            init=true;
        }

        double eval(PropertiedFace *f)
        {
            HelfE e(f);    
            //dA a(f);
            //dV v(f);
   
            double kappa= f->GetProperty("kappa");
            double Are=Area(f);
            double Hint = fi.Gauss7PointSub(&e);
            //double Vint = fi.Gauss7Point(&v);
            //double Aint = fi.Gauss7Point(&a);
            double A0=f->GetProperty("A0");
            double Ka=f->GetProperty("Ka");
            if(A0==0) {f->SetProperty("A0",Are);A0=Are;};
            //return Hint*kappa + Vint*lambdaVol*(V - V0)*(V - V0) + Aint*lambdaSur*(A - A0)*(A - A0);
            return Hint*kappa*0.5+Ka*abs(A0-Are);
        }
        
        void deep_props(vector<string> *obs, vector<int> *ind)
        {
            if(!obs->empty()) obs->clear(); if(!ind->empty()) ind->clear();
            obs->push_back("Area");
            obs->push_back("Volume");
            ind->push_back(0);
            ind->push_back(1);
        }
        
        double deep_eval(PropertiedFace *f, bool trye=false)
        {
            double Volu=VIncr(f); // TODO Put all increments into one
            double Are=Area(f);
            double energy= eval(f);

            if(trye)
            {
           f->SetProperty("Try_Energy",energy);
           f->SetProperty("Try_Volume",Volu);
  //         cout << Vol << " "<< f->GetProperty("Try_Volume") << endl;
 
           f->SetProperty("Try_Area",Are);
           }
            else
            {
           f->SetProperty("Energy",energy);
           f->SetProperty("Volume",Volu);
//                  cout << Vol << " "<< f->GetProperty("Volume") << endl;
           f->SetProperty("Area",Are);
            }
          return energy;
        }
        
        double evalGradient(PropertiedFace *f, int p)
        {
            dInt h(f,p,f->GetProperty("kappa") , lambdas[1], lambdas[0], *V, defaults[1], *A, defaults[0]);
            return fi.Gauss7PointSub(&h);
        };

        double computeTotalEnergy()
        {
            return 0.5*lambdas[1]*(*V - defaults[1])*(*V - defaults[1]) + 0.5*lambdas[0]*(*A - defaults[0])*(*A - defaults[0]);
        }
      
        class dInt : public Funct
        {
        public:
            dInt(PropertiedFace *nga, int p, double k, double l1, double l2,
                    double Vol, double Vol0, double Ar, double Ar0)
            {
                ngt = nga;
                param = p;
                kappa = k;
                lVol = l1;
                lSur = l2;
                volume = Vol;
                initialVolume = Vol0;
                area = Ar;
                initialArea = Ar0;
            }

            double doit(double eta,double xi)
            {
                double tol = 1e-14;
                triple n   = ngt->normal(eta, xi);
                triple xe  = ngt->x_eta(eta, xi);
                triple xx  = ngt->x_xi(eta, xi);
                triple dxe = ngt->x_eta_p(eta,xi,param);
                triple dxx = ngt->x_xi_p(eta,xi,param);
                triple dn  = ngt->normal_p(eta,xi,param);
                triple crss = xe%xx;
                double xexe = xe*xe;
                double dxexe = 2*(dxe*xe);
                double xxxx = xx*xx;
                double dxxxx= 2*(dxx*xx);
                double xexx = xe*xx;
                double dxexx = (dxe*xx)+(xe*dxx);
                double nc0 = n*(ngt->c[0]);
                double dnc0 = dn*(ngt->c[0])+n*ngt->c_p(0,param);
                double nc1 = n*(ngt->c[1]);
                double dnc1 = dn*(ngt->c[1])+n*ngt->c_p(1,param);
                double detI  = xexe * xxxx - xexx * xexx;
                double ddetI = dxexe*xxxx + xexe*dxxxx - dxexx*xexx - xexx*dxexx;

                // det of 1st fundamental form
                // ddetI = dxexe*xxxx+xexe+dxxxx-dxexxx*xexx-xexx-dxexx;   
                if (detI < tol and detI > - tol)
                {
                    cout << "Unable to calculate Mean curvature, determinant of the 1st fundamental form is zero!";
                    return 0;
                }

                triple innerc  = ngt->c[2] - ngt->c[1] - ngt->c[0];
                triple dinnerc = ngt->c_p(2,param) - ngt->c_p(1,param) - ngt->c_p(0,param);
                double dinner  = dn*innerc + n*dinnerc;
                double H   = 2*(nc0)*xxxx + 2*(nc1)*xexe - 2*(n*(innerc))*xexx;
                double Hm  = H / (2*detI);
                double dH  = 2*dnc0*xxxx + 2*nc0*dxxxx + 2*dnc1*xexe + 2*nc1*dxexe - 2*( dinner*xexx + (n*innerc)*dxexx );
                double u   = 2*detI;
                double dHm = dH/(u) - H/(u*u)*(2*ddetI);

                double dS = crss.abs();
                double ddS = 1/crss.abs()*(ngt->cross_p(eta,xi,param)*crss);

                double dHelfrich = kappa*((2*Hm)*dHm*dS + Hm*Hm*ddS);

                double dVolume = 2.0/3.0*(volume - initialVolume)*(ngt->eval_p(eta, xi, param)*n
                        + ngt->eval(eta,xi)*dn)*dS
                    + 1.0/3.0*(volume - initialVolume)*(volume - initialVolume)*(ngt->eval(eta,xi)*n)*ddS;

                double dSurface = 2*(area - initialArea)*dS + (area - initialArea)*(area - initialArea)*ddS;

                double result = dHelfrich + 0.5*lVol*dVolume + 0.5*lSur*dSurface;

                return result;
            }

            PropertiedFace *ngt;
            int param;
            double kappa;
            double lVol;
            double lSur;
            double volume;
            double initialVolume;
            double area;
            double initialArea;
        };

};
