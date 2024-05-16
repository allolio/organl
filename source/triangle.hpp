

bool DostalikHalf(triple *vertex1, triple *vertex2, triple *n1, triple *n2, triple *c1, triple *c2)
{
    double c=(*n1)*(*n2);
    triple d=(*vertex2)-(*vertex1);
    double s=0.5;        double nad=*n1*d;
    //cout << " C: "<< c << endl; 
    if(c<NAGATA_PARALLEL)
    {
        double nbd=*n2*d; double cfac=1./(1.-c*c);
        double snad=nad*s; double snbd=(2.-3.*s)*nbd;
        *c1=(*n1*(snad-c*snbd) + *n2*(-c*snad+snbd))*cfac;
        *c2=d*(1-2*s)-*c1;
    }
    else
    {
        *c1 = (*n1)*(0.5*nad);
        *c2 = *c1*-1.;
    }
    return true;  
}

bool EdgeGradient(triple *d,triple *c0, double t, triple *grad)
{
    *grad=*d+(* c0 * (2*t-1));
    return true; 
}


#define PARALLEL_GRAD 0.001
bool BrokenCurveNormal(triple *d, triple *c0, triple *c1, triple *bcn)
{
    triple g0; triple g1;
    triple dhalf=*d *0.5;
    EdgeGradient(&dhalf,c0,1,&g0);
    EdgeGradient(&dhalf,c1,0,&g1);
    *bcn=g0%g1;
    double bcna=bcn->abs();
    *bcn=*bcn/bcna;
   // cout << bcna << endl;
    if(bcna < PARALLEL_GRAD) return false;
    return true;
}

struct NagataTriangle : public PropertiedFace, public Neighbor
{
    /**
     *
     * Nagata interpolation replaces the linear edge with a quadratic edge.  To do this it uses the normals
     * in the edge vertices, the resulting new edge is then perpendicular to both normals. Further, the interpolation
     * makes the curvatures additional term c the smallest possible. Equivalently, the c vector is the only linear
     * combination of normal vectors that satisfies the requirement of orthogonality.
     * https://www.sciencedirect.com/science/article/pii/S016783960500021X
     *
     * If the normals are parallel, there is nothing to do and the edge remains straight.
     * The nagata algorithm fails when (n1*dp)(n2*dp) > 0 which corresponds to the case of existence of inflection point
     * on the curved edge. This is however impossible to cover with a quadratic mesh. Therefore, the linear edge needs to be used
     * again.
     * Another case of failure is when one of the normals is nearly perpendicular and the second is far from being perpendicular.
     * This will result in a strong localized change in curvature. We will use the linear edge here as well.
     * https://www.sciencedirect.com/science/article/pii/S0010448512002655
     *
     * */
 // Attributes
           double q[3];
    NagataTriangle()
    {
        n.resize(3);
        c.resize(3);
        v.resize(3);
        firstprep=true;
    }
    
    bool updatepair(int i,int j)
    {
            triple n1;
            triple n2;
            triple p1;
            triple p2;

            if (i == 2)
            {
                n1 = *n[j];
                n2 = *n[i];
                p1 = *v[j];
                p2 = *v[i];
            }
            else
            {
                n1 = *n[i];
                n2 = *n[j];
                p1 = *v[i];
                p2 = *v[j];
            }

            triple dp  = p2 - p1;
            //triple dp  = p1 - p2;
            triple nu  = (n1 + n2)/2;
            triple dnu = (n1 - n2)/2;
            double d=dp*nu;
            double dd=dp*dnu;
            double dc=n1*dnu;
            /*
            if (dc >= 0.5)
            {
                cout << "Normals have opposite directions!!\n";
            }
            */
            double c=1-2*dc;                 // = n1*n2
            triple b   = dp / dp.abs();
            double nb1 = n1*b;
            double nb2 = n2*b;
            triple co; co=co*0;
            if ( (abs(c) < NAGATA_PARALLEL)
                    and (nb1*nb2 <= 0)          // inflection point
                    //and (not ((abs(nb1) < NAGATA_EPS1 or abs(nb2) < NAGATA_EPS1) and (abs(nb1 + nb2) > NAGATA_EPS2)))
               )
            { co=nu*(dd/(1-dc))+dnu*(d/dc); this->c[i]=co; return true; }
            this->c[i]=co;    
            return false;
    }
    bool updatec() // Update Coefficents : To be called after Mesh update.
    {
        for(int i=0;i<3;i++)
        {
            int j=i+1; if(j>2) j=0;
            updatepair(i, j);
            // int ci = (i == 0 ? (j == 1 ? 0 : 2) : 1);
        }
       // if(firstprep) firstprep=false;
      //  mid=eval(0.667, 0.333);
        return true;
    }

    inline triple eval(double eta,double xi)
    {
            return (*v[0]*(1-eta))+(*v[1]*(eta-xi))+*v[2]*xi-c[0]*(1-eta)
               *(eta-xi)-c[1]*((eta-xi)*xi)-c[2]*((1-eta)*xi);
    }
    inline triple x_eta(double eta, double xi)
    {
        return ((*v[1] - *v[0]) + c[0]*(2*eta - xi - 1)) + (c[2] - c[1])*xi;
    }

    inline triple x_xi(double eta, double xi)
    {
        return ((*v[2] - *v[1]) + c[1]*(2*xi - eta)) + (c[0] - c[2])*(1 - eta);
    }

    inline triple normal(double eta, double xi)
    {
        triple direction=NagataTriangle::x_eta(eta, xi) % NagataTriangle::x_xi(eta, xi); // LOL CRITICALLY IMPORTANT!
        return direction/direction.abs();
    }
    duple normang(const triple &n)
    {
        duple rd;
        rd.phi=acos(n.z);
        rd.theta=atan2(n.y,n.x);
        //if (n.x<0) rd.theta=2*M_PI-rd.theta;
        return rd;
    }
    triple angnorm(duple &ang)
    {
        triple rn;
        rn.x=cos(ang.theta)*sin(ang.phi);
        rn.y=sin(ang.theta)*sin(ang.phi);
        rn.z=cos(ang.phi);
        return rn;
    }
    // Structure of Parameters for face v1 0-2 , v2 3-5, v3 6-8,n1,9-10,n2, 11-12, n3 13-14;
    // Structure of Grad Matrix dc col: dc1.x, dc1.y, dc1.z [0-2] etc. dp = row
    inline triple colvec(int cols,int row)
    {
        triple rv;
        rv.x=dcdp[cols][row];
        rv.y=dcdp[cols+1][row];
        rv.z=dcdp[cols+2][row];
        return rv;
    }
    inline triple c_p(int cnum,int p)
    {
        return colvec(cnum*3,p);
    }
    inline triple eval_p(double eta,double xi, int p)
    {

       // dev/d_v1 component j  for component k
       // -dc_1k/dvi_j*(1-eta)*(eta-xi)-dc_2k/dvi_j*((eta-xi)*xi)-dc_3k/dvi_j*((1-eta)*xi)
       //  Same for all vectors and even dn_s
       triple c1p=colvec(0,p);
       triple c2p=colvec(3,p);
       triple c3p=colvec(6,p);
       triple rv=c1p*(1-eta)*(eta-xi)-c2p*((eta-xi)*xi)-c3p*((1-eta)*xi);
       double specific=0;
       if(p < 3) specific=(1-eta);
       if(p > 2 && p < 6) specific=(eta-xi);
       if(p > 5 && p < 9) specific=xi;
       if(p > 8 ) return rv;
       rv.p[p%3]+=specific;

       return rv;
    }
    inline triple x_eta_p(double eta, double xi, int p)
    {
        // *v[1] - *v[0] + c[0]*(2*eta - xi - 1) + (c[2] - c[1])*xi;
       triple c1p=colvec(0,p);
       triple c2p=colvec(3,p);
       triple c3p=colvec(6,p);
       //  deveta_k/dv_j = dc_2k/dv_j*(2*eta - xi - 1) + (dc_3k/dv_j - dc_1k_dv_j)*xi
       triple rv= c1p*(2*eta-xi-1)+(c3p-c2p)*xi;
       double specific;
       if(p < 3) specific=-1;
       if(p > 2 && p < 6) specific=1;
       if(p > 5) return rv;
       rv.p[p%3]+=specific;

       return rv;
    }

    inline triple x_xi_p(double eta, double xi, int p)
    {
        //  *v[2] - *v[1] + c[1]*(2*xi - eta) + (c[0] - c[2])*(1 - eta);
       triple c1p=colvec(0,p);
       triple c2p=colvec(3,p);
       triple c3p=colvec(6,p);
       //  devxi_k/dv_j = dc_2k/dv_j (2*xi-eta)+(dc_1k/dv_j-dc_3k/dv_j)(1-eta)
       triple rv= c2p*(2*xi-eta)+(c1p-c3p)*(1-eta);
       double specific=0;
 //      if(p < 3) specific=0;
       if(p > 2 && p < 6) specific=-1;
       if(p > 5 && p < 9) specific=1;
       if(p > 8 ) return rv;
       rv.p[p%3]+=specific;
       //rv.x+=specific;
       //rv.y+=specific;
       //rv.z+=specific;

           //return rv;
       return rv;
    }

    inline triple cross_p(double eta, double xi, int p)
    {
       triple direction = x_eta_p(eta, xi,p) % x_xi(eta, xi)+x_eta(eta, xi) % x_xi_p(eta, xi,p) ;
       return direction;
    }

    inline triple normal_p(double eta, double xi, int p)
    {
        triple crss = x_eta(eta,xi) % x_xi(eta,xi);
        triple dcross= cross_p(eta,xi,p);
        //return dcross*((crss*dcross)/crss.abs());//...*(cross(eta,xi) ...
        //maybe this one:
        return crss*(-(crss*dcross)/pow(crss.abs(),3)) + dcross/crss.abs();
    }
   bool firstprep;

    bool gradcol(int col);
    bool numgradcol(int col);
    inline triple phidiff(duple &n1, duple &n2, triple &l,int ent,int nb, double step=1e-7);

    bool prepgrad()
    {
        if(firstprep)
        {
            dcdp.resize(9);
            for(int i=0;i!=9;i++) dcdp[i].resize(15);
        }
        firstprep=false;
        for(int i=0;i!=3;i++) numgradcol(i);
    // Debug;
/*
    cout << std::fixed;
    cout << std::setprecision(4);
    cout << "------------" <<endl;
    c[0].print();
    c[1].print();
    c[2].print();
    cout << "------------" <<endl;

    for (int i=0;i!=dcdp[0].size();i++)
    {
        for(int j=0;j!=dcdp.size();j++)
            cout << " j: " << j <<" "<< dcdp[j][i] << " " ;
        cout << endl;
    }
    cout << endl;
    cout << "------------" <<endl;
*/

        return true;
    }

    vector< vector < double > > dcdp; // Matrix for gradient

    inline void giveC(double eta, double xi, vector<triple> *cOut)
    {
        *cOut = this->c;
    }
   bool getbound(vector<triple*> *bounds) 
   {
      if(bounds->size()!=4) bounds->resize(4);
     for(int i=0;i!=3;i++) bounds->at(i)=v[i];
       bounds->at(3)=&mid;
       return true;
   }
   triple mid;
};


void Thirdorder(triple *vertex1, triple *vertex2, triple *n1, triple *n2, triple *coeff1, triple *coeff2)
{
    triple dist=*vertex2-*vertex1; triple c1=(dist % *n1 );triple c2=(dist % *n2 );
    double a=*n1 * *n2;
    double b=(dist % *n2 )* *n1;
    double c=c1*c1;
    double d=c2*c1;
    double nad1=*n1 *dist;
    double nad2=- (*n2 *dist);
    double e=c2*c2; 
    double k=81*b*b*b*b-18*(5*c-9*a*d+5*e)*b*b+(9*a*a-10)*(9*d*d-10*c*e);
    double inter=27*d*(b*b+a*d)-30*a*c*e;
    double r[2][4];
    double med[4];
    
    r[0][0]=-5*(9*e*b*b +9*d*d-10*c*e); 
    r[0][1]=inter;
    r[0][2]=45*b*(d-a*e);
    r[0][3]=3*b*(9*b*b-10*c+9*a*d);
    r[1][0]=inter;
    r[1][1]=20*c*e-18*(c*b*b+d*d);
    r[1][2]=-(3*b*(9*b*b+9*a*d-10*e));
    r[1][3]=18*b*(a*c-d);
    for(int i=0;i!=4;i++)
    {   
        med[i]=r[0][i]*nad1+r[1][i]*nad2;
    }
    *coeff1= (*n1)*med[0]+(*n2)*med[1]+c1*med[2]+c2*med[3];
    *coeff2= *coeff1+((*n2)*med[1])+c2*med[3]; *coeff2= *coeff2/k; //(*n1)*med[0]+(*n2)*med[1]*2.+c1*med[2]
    *coeff1= *coeff1/k;
}

struct GNagataTriangle : public NagataTriangle
{
    /**
     *
     * GNagataTriangle has the same properties and functions as the standarn NagataTriangle.
     * The difference is, that GNagataTriangle will generate third order edge interpolants for those triangles
     * that fail the conditions outlined in the manuscript by Neto et al. 
     *
     */
    bool healthyTriangle;
    bool healthyEdge[3];
    vector<triple> cbase;
    GNagataTriangle ()
    {
        healthyTriangle = true;
     //   for (int i = 0; i < 4; i ++)
        cbase.resize(10);
       }
    
    bool updatepair(int i,int j)
    {
            if ((*n[i]*(*v[i] - *v[j])) * (*n[j]*(*v[i] - *v[j])) > 0)
            {    healthyEdge[i]=false;
                 healthyTriangle=false;
                 c[i].x=INVALID_EDGE;
                 return updatec();
            }
            else return NagataTriangle::updatepair(i,j);
    }

    bool updatec()
    {
        // assume nagata representibility in the beginning
        healthyTriangle = true;
        for (int i = 0; i < 3; i ++)
        {
            healthyEdge[i]=true;
            int j = (i + 1) % 3;
            if ((*n[i]*(*v[i] - *v[j])) * (*n[j]*(*v[i] - *v[j])) > 0)
            {
                healthyTriangle = false;
                healthyEdge[i] = false;
            }
        }
        if(healthyTriangle)
        {
//          cbase.clear();
          bool rv= NagataTriangle::updatec();
          mid=eval(0.667, 0.333);
          return rv;
        }
        NagataTriangle::updatec();
//        if(firstprep) {}

         
        if (not healthyTriangle)    // division into subtriangles needed
        {
          //  cout << "HAPPENING" << endl;
            for(int i=0;i<3;i++)
            {   
            int j=i+1; if(j>2) j=0;
            if(!healthyEdge[i]) 
            {
             Thirdorder(v[i], v[j],n[i],n[j],&cbase[i],&cbase[i+3]);   
            }
            else {cbase[i]=NagataTriangle::c[i];cbase[i+3]=cbase[i+3]*0;}
            // int ci = (i == 0 ? (j == 1 ? 0 : 2) : 1);
            }        
            // Now reconstruct the surface polynomial
        //    cbase[6]=(cbase[0]*-3 - cbase[3] + cbase[1]*-2   -cbase[4]     + cbase[2]*3  + cbase[5])*0.25;
          //  cbase[7]=(cbase[0]*-1 + cbase[3] + cbase[1]*2  +cbase[4]     + cbase[2]    - cbase[5])*0.25;
          // ***  cbase[8]=(cbase[0]*-1 - cbase[3]*3+ cbase[1]*-1  +cbase[4]*-3  + cbase[2]    + cbase[5]*3)*.25;
         //   cbase[8]=(cbase[0]*-1 - cbase[3]*3+ cbase[1]*-2  +cbase[4]*-3  + cbase[2]    + cbase[5]*3)*.25;
        //    cbase[9]=(cbase[0]    - cbase[3] + cbase[1]*2    - cbase[4]    - cbase[2]    + cbase[5])*0.25;
            
            cbase[6]=cbase[2]-cbase[0]-cbase[1];
            cbase[7]=cbase[1];
            cbase[8]=cbase[5]-cbase[3]-cbase[4];
            cbase[9]=cbase[9]*0;
 //           cbase[7]=(cbase[0]*-1 + cbase[3] + cbase[1]*2  +cbase[4]     + cbase[2]    - cbase[5])*0.25;
          // ***  cbase[8]=(cbase[0]*-1 - cbase[3]*3+ cbase[1]*-1  +cbase[4]*-3  + cbase[2]    + cbase[5]*3)*.25;
   //         cbase[8]=(cbase[0]*-1 - cbase[3]*3+ cbase[1]*-2  +cbase[4]*-3  + cbase[2]    + cbase[5]*3)*.25;
     //       cbase[9]=(cbase[0]    - cbase[3] + cbase[1]*2    - cbase[4]    - cbase[2]    + cbase[5])*0.25;
        }  
        mid=eval(0.667, 0.333);
        /*
        cout << "CBASE" << endl;
        for(int i=0;i!=10;i++) cbase[i].print();
        
        for(int i=0;i!=100;i++)
        {
            double t=1./100*i;
            triple f=(*v[0]+(*v[1]-*v[0]-cbase[0]-cbase[3])*t  + cbase[0]*t*t + cbase[3]*t*t*t);
            cout << t << " " << f.x << " "<< f.y << " " << f.z << " " << eval(t,0).x << " " << eval(t,0).y << " " << eval(t,0).z << " " << endl;
        }
           for(int i=0;i!=100;i++)
        {
            double t=1./100*i;
            triple f=(*v[1]+(*v[2]-*v[1]-cbase[1]-cbase[4])*t  + cbase[1]*t*t + cbase[4]*t*t*t);
            cout << t << " " << f.x << " "<< f.y << " " << f.z << " " << eval(1,t).x << " " << eval(1,t).y << " " << eval(1,t).z << " " << endl;
        }
        
              for(int i=0;i!=100;i++)
        {
            double t=1./100*i;
            triple f=(*v[2]+(*v[0]-*v[2]-cbase[2]-cbase[5])*t  + cbase[2]*t*t + cbase[5]*t*t*t);
            cout << t << " " << f.x << " "<< f.y << " " << f.z << " " << eval(t,t).x << " " << eval(t,t).y << " " << eval(t,t).z << " " << endl;
        }
     */
        return true;
    }
    
    //inline(

  
    // number of functions rewriting the standard NagataTriangle counterparts
    inline triple eval(double eta, double xi)
    {
        if (healthyTriangle)
            return NagataTriangle::eval(eta, xi);

      // c00 + c10eta + c01x + c11etaxi +c20eta**2+c02*xi**2
      // c12*eta*xi**2 + c30*eta**3 + c03*xi**3
        
        return *v[0]+(*v[1]-*v[0]-cbase[0]-cbase[3])*eta
               +(*v[2]-*v[1]+cbase[0]-cbase[2]+cbase[3]-cbase[5])*xi
               +cbase[6]*eta*xi+cbase[0]*eta*eta+cbase[7]*xi*xi
               +cbase[8]*eta*eta*xi+cbase[9]*eta*xi*xi
               +cbase[3]*eta*eta*eta+cbase[4]*xi*xi*xi;
    }

    inline triple x_eta(double eta, double xi)
    {
        if (healthyTriangle)
            return NagataTriangle::x_eta(eta, xi);

     //   return *v[0]+(*v[1]-*v[0]-cbase[0]-cbase[3])*eta
     //          +(*v[2]-*v[1]+cbase[1]-cbase[2]+cbase[3]-cbase[5])*xi
     //          +cbase[6]*eta*xi+cbase[0]*eta*eta+cbase[7]*xi*xi
     //          +cbase[8]*eta*eta*xi+cbase[9]*eta*xi*xi
     //          +cbase[3]*eta*eta*eta+cbase[4]*xi*xi*xi;
     return (*v[1]-*v[0]-cbase[0]-cbase[3])
               +cbase[6]*xi+cbase[0]*2*eta
               +cbase[8]*2*eta*xi+cbase[9]*xi*xi
               +cbase[3]*3*eta*eta;
    }
// eta xi
//
//               +cbase[6]
//               +cbase[8]*2*eta+cbase[9]*2*xi
//
    // eta eta
//
//               cbase[0]*2
//               +cbase[8]*2*xi+
//               +cbase[3]*6*eta;
//
    inline triple x_xi(double eta, double xi)
    {
        if (healthyTriangle)
            return NagataTriangle::x_xi(eta, xi);
         return      (*v[2]-*v[1]+cbase[0]-cbase[2]+cbase[3]-cbase[5])
               +cbase[6]*eta+cbase[7]*2*xi
               +cbase[8]*eta*eta+cbase[9]*eta*xi*2
               +cbase[4]*3*xi*xi;
   
    }
//  xi eta
//  
//               +cbase[6]  +cbase[8]*2*eta+cbase[9]*xi*2
// 
//
        // xi xi
//
//               +cbase[7]*2
//               +cbase[9]*eta*2
//               +cbase[4]*6*xi;
//
//
//


    inline triple normal(double eta, double xi)
    {
        if (healthyTriangle)
            return NagataTriangle::normal(eta, xi);
        triple direction=x_eta(eta, xi) % x_xi(eta, xi); 
        return direction/direction.abs();
       
    }

    inline void giveC(double eta, double xi, vector<triple> *cOut)
    {
        if (healthyTriangle)
            return NagataTriangle::giveC(eta, xi, cOut);
        else
        {cOut->resize(10);
            for(int i=0;i!=10;i++)
                cOut->at(i)=cbase[i];
        }
    }
};


#ifndef HAVE_GRADIENT
bool NagataTriangle::gradcol(int col)
{return true;}
bool NagataTriangle::numgradcol(int col)
{return true;}
    inline triple phidiff(duple &n1, duple &n2, triple &l,int ent,int nb, double step=1e-5){ triple rv; return rv;}
#endif
