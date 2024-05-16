struct triple
{ 
    triple()
    {}

    triple(double *in)
    {
        x=*in;
        y=*(in+1);
        z=*(in+2);
    }

    triple(double a, double b, double c)
    {
    x=a;y=b;z=c;    
    }

    union
    {
        double p[3];
        struct
        {
            double x;
            double y;
            double z;
        };
    };   

    triple operator-(const triple &b)
    {
        triple t;
        t.x=this->x-b.x;
        t.y=this->y-b.y;
        t.z=this->z-b.z;
        return t;
    }

    triple operator%(const triple &b)
    {
       triple t;
        t.x=this->y*b.z-this->z*b.y;
        t.y=this->z*b.x-this->x*b.z;
        t.z=this->x*b.y-this->y*b.x;
        return t;
    }

    triple operator+(const triple &b)
    {
        triple t;
        t.x=this->x+b.x;
        t.y=this->y+b.y;
        t.z=this->z+b.z;
        return t;
    }

    triple operator*(double b)
    {
        triple t;
        t.x=this->x*b;
        t.y=this->y*b;
        t.z=this->z*b;
        return t;
    }

    triple operator/(double b)
    {
        triple t;
        t.x=this->x/b;
        t.y=this->y/b;
        t.z=this->z/b;
        return t;
    }
  
     // Scalar Product
     double operator*(const triple &b)
     {
        double t=0;
        t+=x*b.x;
        t+=y*b.y;
        t+=z*b.z;
        return t;
     }

    double abs()
    {
        return sqrt(x*x+y*y+z*z);
    }

    triple positive()
    {
        triple r;
        r.x=fabs(this->x);
        r.y=fabs(this->y);
        r.z=fabs(this->z);
        return r;
    }

    double distanceSquared(triple &t)
    {
        double res = 0;
        res += (this->x - t.x)*(this->x - t.x);
        res += (this->y - t.y)*(this->y - t.y);
        res += (this->z - t.z)*(this->z - t.z);
        return res;
    }

    double distance(triple &t)
    {
        return sqrt(distanceSquared(t));
    }

    void print()
    {
        cout << x << " " << y << " " << z << endl;
    }
};


struct Matrix3
{
    triple row[3];

    double p(int i,int j)
    {
        return this->row[i].p[j];
    }

    void setp(int i,int j,double val)
    {
        this->row[i].p[j]=val;
    }

    Matrix3 operator*(double l)
    { 
        Matrix3 mat;
        mat.row[0]=this->row[0]*l;
        mat.row[1]=this->row[1]*l;
        mat.row[2]=this->row[2]*l;
        return mat;    
    }

    Matrix3 operator/(double l)
    { 
        Matrix3 mat;
        mat.row[0]=row[0]/l;
        mat.row[1]=row[1]/l;
        mat.row[2]=row[2]/l;
        return mat;    
    }

    triple operator*(triple x)
    { 
        triple a;
        a.x = row[0]*x;
        a.y = row[1]*x;
        a.z = row[2]*x;
        return a;
    }

    double det()
    {
        Matrix3 *m=this;
        double dval= m->p(0, 0) * (m->p(1, 1) * m->p(2, 2) - m->p(2, 1) * m->p(1, 2)) -
             m->p(0, 1) * (m->p(1, 0) * m->p(2, 2) - m->p(1, 2) * m->p(2, 0)) +
             m->p(0, 2) * (m->p(1, 0) * m->p(2, 1) - m->p(1, 1) * m->p(2, 0));
	    return dval;
    }

    Matrix3 inv()
    {
        // computes the inverse of a matrix m
        double invdet = 1 / det();
        Matrix3 *m=this;
        Matrix3 minv; // inverse of matrix m
        minv.setp(0, 0, (m->p(1, 1) * m->p(2, 2) - m->p(2, 1) * m->p(1, 2)) * invdet);
        minv.setp(0, 1, (m->p(0, 2) * m->p(2, 1) - m->p(0, 1) * m->p(2, 2)) * invdet);
        minv.setp(0, 2,(m->p(0, 1) * m->p(1, 2) - m->p(0, 2) * m->p(1, 1)) * invdet);
        minv.setp(1, 0, (m->p(1, 2) * m->p(2, 0) - m->p(1, 0) * m->p(2, 2)) * invdet);
        minv.setp(1, 1, (m->p(0, 0) * m->p(2, 2) - m->p(0, 2) * m->p(2, 0)) * invdet);
        minv.setp(1, 2, (m->p(1, 0) * m->p(0, 2) - m->p(0, 0) * m->p(1, 2)) * invdet);
        minv.setp(2, 0, (m->p(1, 0) * m->p(2, 1) - m->p(2, 0) * m->p(1, 1)) * invdet);
        minv.setp(2, 1, (m->p(2, 0) * m->p(0, 1) - m->p(0, 0) * m->p(2, 1)) * invdet);
        minv.setp(2, 2, (m->p(0, 0) * m->p(1, 1) - m->p(1, 0) * m->p(0, 1)) * invdet);

        return minv;
    }

    void print(void)
    {
        row[0].print();
        row[1].print();
        row[2].print();
    }

    Matrix3 transpose(void)
    {
        Matrix3 m=(*this);
        m.setp(0,1 , this->p(1,0));
        m.setp(1,0 , this->p(0,1));
        m.setp(2,0 , this->p(0,2));
        m.setp(2,1 , this->p(1,2));
        m.setp(1,2 , this->p(2,1));
        m.setp(0,2 , this->p(2,0));
        return m;  
    }
};  

union duple{
  double p[2];
  struct {
  double phi;
  double theta;
  };
};

/*  void T(void)
  {
    Matrix3 temp;
    for(int i=0;i!=3;i++)
    {
      row[i];
    }
  }*/
  
