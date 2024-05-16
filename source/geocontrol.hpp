#pragma once


#define INVALID_BULK -1

struct GeoControl {
    GeoControl() {}
    FaceIntegrator fi;

    static double gaussCurvature(double eta, double xi, PropertiedFace &p) {
        vector <triple> c(3);
        double tol = 1e-12;
        p.giveC(eta, xi, &c);
        triple n   = p.normal(eta, xi);
        triple xe  = p.x_eta(eta, xi);
        triple xx  = p.x_xi(eta, xi);
        double nc0 = n*c[0];
        double nc1 ;//= n*c[1];
        double detI = (xe*xe) * (xx*xx) - (xe*xx) * (xe*xx);                        // det of 1st fundamental form
        double detII;
        double cross;

        if (c.size()==10) {
            cross = (c[6]+c[8]*2*eta+c[9]*2*xi)*n;
            nc0=nc0+(c[8]*xi+c[3]*3*eta)*n; //
            nc1=(c[7]+c[9]*eta+c[4]*3*xi)*n;
            detII=(4*nc0*nc1 - cross*cross);
        }
        else {
            double nc2 = n * c[2];
            nc1 = n*c[1];
            cross=n*(c[2] - c[1] - c[0]); //2nd fundamental cross term
            detII = 2*(nc2)*(nc0 + nc1) + 2*nc0*nc1 - nc0*nc0 - nc1*nc1 - nc2*nc2;  // det of 2nd fundamental form
        }

        if (detI < tol and detI > - tol) {
            cout << "Unable to calculate Gauss curvature, determinant of the 1st fundamental form is zero!";
            return 0;
        }

        return detII / detI;
    }


    static double meanCurvature(double eta, double xi, PropertiedFace &p) {
        vector <triple> c(3);
        p.giveC(eta, xi, &c);
        double tol = 1e-12;
        triple n   = p.normal(eta, xi);
        triple xe  = p.x_eta(eta, xi);
        triple xx  = p.x_xi(eta, xi);
        double xexe = xe*xe;
        double xxxx = xx*xx;
        double xexx = xe*xx;
        double nc0 = n*c[0];
        double nc1;
        double cross;

//      for (int i=0;i!=c.size();i++) c[i].print(); exit(0);

        if (c.size()==10) {
            cross = (c[6]+c[8]*2*eta+c[9]*(2*xi))*n; //  +cbase[6]  +cbase[8]*2*eta+cbase[9]*xi*2
            nc0=nc0+(c[8]*xi+c[3]*3*eta)*n; //II(0,0)
            nc1=(c[7]+c[9]*eta+c[4]*3*xi)*n; // I (1,1)
        }
        else {
            nc1 = n * c[1];
            cross = n * (c[2] - c[1] - c[0]); //2nd fundamental cross term
        }
                // eta eta
//
//               cbase[0]*2
//               +cbase[8]*2*xi+
//               +cbase[3]*6*eta;

        double detI  = xexe * xxxx - xexx * xexx;   // det of 1st fundamental form

        if (detI < tol and detI > - tol) {
            cerr << "Unable to calculate Mean curvature, determinant of the 1st fundamental form is zero!";
            return 0;
        }

        double H = 2*(nc0)*xxxx + 2*(nc1)*xexe - 2*(cross)*xexx;
        // double dH=2*dnc0*xxxx+2*nc0*dxxxx+2*dnc1*xexe+2*nc1*dxexe-(2*((dinner)*xexx+inner*dxexx)
        //dinner=dn*(pc2-pc1-pc0)+n*(dpc2-dpc1-dpc0)//
        // dhelf=H^2/dp-->2*H*H'
        //cout << "Mean: " << H / (2*detI) << endl;

        return H / (2*detI);
    }


    double Area(PropertiedFace *p) {
        dA a(p);
        //return fi.Gauss7Point(&a);
        return fi.LG7(&a);
    }


    double VIncr(PropertiedFace *p) {
        dV v(p);
        // return fi.Gauss7Point(&v);
        return fi.LG7(&v);
    }


    double DebugI(PropertiedFace *p) {
       Debug v(p);
       return fi.Z1Point(&v);
    }


    double KInt(PropertiedFace *p) {
         dG g(p);
         //return fi.LG7(&g);
         return fi.Gauss7PointSub(&g);
    }


    double HInt(PropertiedFace *p)
    {
         dH h(p);
         return fi.LG7(&h);
         // return fi.Gauss7Point(&h);
    }


   //private:
    // Horrible I know, got any better ideas ?
    class dA : public Funct
    {
    public:
    dA(PropertiedFace *nga)
    {ngt=nga;}
    double doit(double eta,double xi)
    {
        return (ngt->x_eta(eta, xi) % ngt->x_xi(eta, xi)).abs();
    }
   PropertiedFace *ngt;
    };

    class dV : public Funct
    {
    public:
    dV(PropertiedFace *nga)
    {ngt=nga;}
    double doit(double eta,double xi)
    {
        return  ngt->eval(eta,xi)*ngt->normal(eta,xi)*(ngt->x_eta(eta, xi) % ngt->x_xi(eta, xi)).abs()/3.;
    }
    PropertiedFace *ngt;
    };

    class dG :
        public Funct
    {
        public: dG(PropertiedFace *nga)
                {
                    ngt = nga;
                }

                double doit(double eta,double xi)
                {
                    return gaussCurvature(eta, xi,*ngt) * (ngt->x_eta(eta, xi) % ngt->x_xi(eta, xi)).abs();
                }

                PropertiedFace *ngt;
    };

    class dH :
        public Funct
    {
        public: dH(PropertiedFace *nga)
                {
                    ngt = nga;
                }

                double doit(double eta,double xi)
                {
                    return meanCurvature(eta, xi, *ngt) * (ngt->x_eta(eta, xi) % ngt->x_xi(eta, xi)).abs();
                }

                PropertiedFace *ngt;
    };

    class Debug :
    public Funct
    {
        public: Debug(PropertiedFace *nga)
                {
                    ngt = nga;
                }

                double doit(double eta,double xi)
                {
                    return  ngt->x_eta(eta, xi).y;// % ngt->x_xi(eta, xi)).abs();
                }

                PropertiedFace *ngt;
    };


};

