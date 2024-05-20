class CustomDepEnergy {
  public:
    virtual double eval(bool tay) = 0;
    virtual bool store() = 0;
    virtual bool restore() = 0;
    void *next; // Interface for linking Energies
    double tray;
    double stored;
    map<string, double>
        props; // This to store, e.g. Geodesic curvature if necessary.
};

class TrivialEdge : public CustomDepEnergy {
  public:
    EdgeProperty *EP;
    double eval(bool tay) override {
      if (tay)
        tray = 1.0;
      else
        stored = 1.0;
      return 1.0;
    }
    bool store() override {
      stored = tray;
      return true;
    }
};

class CustomDepEnergyFactory {
  public:
    virtual CustomDepEnergy *GetNew(void *parent) = 0;
};

class BonnetEdge : public CustomDepEnergy {
  public:
    void setTriangle(NTriangle *n) {
      nt = n;
      v1 = nt->v[edgeindex];
    }
    void setParent(EdgeEnerProperty *eep) { ep = eep; }
    double eval(bool tay) override {
      double kg = nt->GetProperty("kappag");
      double tau = nt->GetProperty("tau");
      double l;
      // 1, Line
      double rv = kg * Fourpoint(&l);
      // 2. Corner Turning for Discontinous tangent.
      Property *e = (Property *)next;
      rv += crosstan(e) * kg;
      // DEBUG

      //       cout << "DBG100" << endl;
      //

      rv += tau * l;
      if (tay)
        tray = rv;
      else
        stored = rv;
      return rv;
    }
    bool store() override {
      stored = tray;
      return true;
    }
    bool restore() override {
      Property *edgep = (Property *)ep->parent;
      NTriangle *ntnew = NULL;
      // Find and Set new triangles
      bool found = false;
      for (int i = 0; i != edgep->deps.size(); i++) {
        if (edgep->deps[i]->ptype->collidable) {
          ntnew = (NTriangle *)edgep->deps[i]->content;
          if (edgep->deps[i]->content == (void *)nt) {
            found = true;
            break;
          }
        }
      }
      // If triangle has changed fix deplist
      if (found != true) {
        ep->deps.clear();
        for (int i = 0; i != edgep->deps.size(); i++) {
          if (edgep->deps[i]->ptype->collidable)
            ep->deps.push_back(edgep->deps[i]);
        }
      }
      triple *v1old = v1;
      if (ntnew != NULL) {
        setTriangle(ntnew);
        // Rebuild Edge index:
        int faceind = -1;
        for (int l = 0; l != 3; l++) {
          if (ntnew->v[l] == v1old) {
            edgeindex = l;
            break;
          }
        }
        return true;
      } else
        return false;
    }
    int edgeindex;

  private:
    double crosstan(Property *e) {
      triple dist;
      triple t1, t2;
      triple norm;
      // Eval my own tan at x=1;
      if (edgeindex == 2) {
        dist = *nt->v[2] - *nt->v[0];
        norm = *nt->n[2];
      } else {
        dist = *nt->v[edgeindex + 1] - *nt->v[edgeindex];
        norm = *nt->n[edgeindex + 1];
      }
      if (nt->healthyEdge[edgeindex])
        t1 = dist + (nt->c[edgeindex]);
      else {
        t1 = dist + nt->c[edgeindex] + (nt->c[edgeindex + 3]) * 2.0;
      }

      // Eval the other tan at x=0;
      NTriangle *nt2 = (NTriangle *)e->deps[0]->content; // Not safe, dude
      int edgeind2 = e->state;
      if (edgeind2 == 2)
        dist = *nt2->v[2] - *nt2->v[0];
      else
        dist = *nt2->v[edgeind2 + 1] - *nt2->v[edgeind2];

      if (nt->healthyEdge[edgeind2]) {
        t2 = dist - (nt2->c[edgeind2]);
      } else {
        t2 = dist - nt2->c[edgeind2] - (nt2->c[edgeind2 + 3]);
      }
      t1 = t1 / t1.abs();
      t2 = t2 / t2.abs();
      // Find orientation.;
      int sgn = 1;
      if (norm * (t1 % t2) < 0)
        sgn = -1;
      return acos(t1 * t2) * sgn;
    }
    double evol(const double x, double *l) {
      // We'll do this the dumb way:
      // d^T/ds = d T^/dt * dt/ds.
      triple t; // Tangent vector
      triple dt;
      triple norm;
      triple c; // curvature vector
      triple dist;
      double dtds;
      if (edgeindex == 2)
        dist = *nt->v[2] - *nt->v[0];
      else
        dist = *nt->v[edgeindex + 1] - *nt->v[edgeindex];

      if (nt->healthyEdge[edgeindex]) {
        t = dist + (nt->c[edgeindex]) * (2 * x - 1);
        dt = nt->c[edgeindex] * (2.0); // a constant!
      } else {
        t = dist + nt->c[edgeindex] * (2 * x - 1) +
            (nt->c[edgeindex + 3]) * (3 * x * x - 1);
        dt = nt->c[edgeindex] * (2.0) + (nt->c[edgeindex + 3]) * (6 * x);
      }
      dtds = 1. / t.abs(); // CHeck this;
      if (edgeindex == 0)
        norm = nt->normal(x, 0); // x 0
      else if (edgeindex == 1)
        norm = nt->normal(0, x); // 0 x
      else
        norm = nt->normal(x, x); // x x
      // https://people.math.wisc.edu/~angenent/561/kgsolutions.pdf
      *l = 1. / dtds; //
      return dt * (norm % t) *
            (dtds * dtds); //*dtds); // Geodesic c// kg Last one removed due to
                            //line element. --> kg dl
      //  else
    }
    double Fourpoint(double *l) {
      // Four-Point Gauss-Legendre Shifted to [0:1] interval
      double rv = 0;
      double lincr = 0;
      rv += 0.173927422568727 * evol(0.0694318442029737, &lincr);
      (*l) += 0.173927422568727 * lincr;
      rv += 0.173927422568727 * evol(0.930568155797026, &lincr);
      (*l) += 0.173927422568727 * lincr;
      rv += 0.326072577431273 * evol(0.330009478207572, &lincr);
      (*l) += 0.326072577431273 * lincr;
      rv += 0.326072577431273 * evol(0.330009478207572, &lincr);
      (*l) += 0.326072577431273 * lincr;
      return rv;
    }
    NTriangle *nt;
    EdgeEnerProperty *ep;
    triple *v1;
    triple c1, c2, di; // dummy edge coefficients;
    bool ishealthy;
};

class BonnetEdgeFactory : public CustomDepEnergyFactory {
  public:
    CustomDepEnergy *GetNew(void *parent) {
    EdgeEnerProperty *ep = (EdgeEnerProperty *)parent;
    BonnetEdge *bne = new BonnetEdge;
    bne->setParent(ep);
    bne->edgeindex = ep->state;
    for (int i = 0; i != ep->deps.size(); i++) {
      if (ep->deps[i]->ptype->collidable) {
        bne->setTriangle((NTriangle *)ep->deps[i]->content);
        break;
      }
    }
    bne->next = NULL;
    return bne;
  }
};

class TauEdge : public CustomDepEnergy {
  public:
    void setTriangle(NTriangle *n) {
      nt = n;
      edgeindex = edgeIndexDetect();
    }
    void setParent(EdgeEnerProperty *eep) { ep = eep; }
    double eval(bool tay) override {
      double tau = nt->GetProperty("tau");
      // 1, Line
      double rv = tau * Fourpoint();
      // 2. Corner Turning for Discontinous tangent.
      if (tay)
        tray = rv;
      else
        stored = rv;
      return rv;
    }
    bool store() override {
      stored = tray;
      return true;
    }
    int edgeindex;

  private:
    int edgeIndexDetect() {
      // Peel Vertex Addresses from mesh.
      FaceVertexMesh *fvm =
          (FaceVertexMesh *)((Property *)(ep->parent))->ptype->parent;
      Edge *edge = (Edge *)((Property *)(ep->parent))->content;
      v1 = &fvm->vertices[edge->vertA].pos;
      v2 = &fvm->vertices[edge->vertB].pos;
      int redgeindex = -1;
      if ((nt->v[0] == v1) or (nt->v[0] == v2))
        if ((nt->v[1] == v2) or (nt->v[1] == v1))
          redgeindex = 0;
        else
          redgeindex = 2;
      else if ((nt->v[1] == v1) or (nt->v[1] == v2))
        if ((nt->v[2] == v2) or (nt->v[2] == v1))
          redgeindex = 1;
        else
          redgeindex = 0;
      return redgeindex;
    }
    double evol(const double x) {
      // We'll do this the dumb way:
      // d^T/ds = d T^/dt * dt/ds.
      triple t; // Tangent vector
      triple dt;
      triple norm;
      triple c; // curvature vector
      triple dist;
      double dtds;
      //  nt->SetProperty("DBG",100);
      if (edgeindex == 2)
        dist = *nt->v[2] - *nt->v[0];
      else
        dist = *nt->v[edgeindex + 1] - *nt->v[edgeindex];

      if (nt->healthyEdge[edgeindex]) {
        t = dist + (nt->c[edgeindex]) * (2 * x - 1);
      } else {
        t = dist + nt->c[edgeindex] * (2 * x - 1) +
            (nt->c[edgeindex + 3]) * (3 * x * x - 1);
      }
      dtds = 1. / t.abs(); // CHeck this;

      // https://people.math.wisc.edu/~angenent/561/kgsolutions.pdf
      //
      return 1. / dtds; //*dtds); // Geodesic c// kg Last one removed due to line
                        //element. --> kg dl
                        //  else
    }
    double Fourpoint() {
      // Four-Point Gauss-Legendre Shifted to [0:1] interval
      double rv = 0;

      rv += 0.173927422568727 * evol(0.0694318442029737);
      rv += 0.173927422568727 * evol(0.930568155797026);
      rv += 0.326072577431273 * evol(0.330009478207572);
      rv += 0.326072577431273 * evol(0.330009478207572);
      return rv;
    }
    bool restore() override {
      Property *edgep = (Property *)ep->parent;
      NTriangle *ntnew = NULL;
      // Find and Set new triangles
      bool found = false;
      for (int i = 0; i != edgep->deps.size(); i++) {
        if (edgep->deps[i]->ptype->collidable) {
          ntnew = (NTriangle *)edgep->deps[i]->content;
          if (edgep->deps[i]->content == (void *)nt) {
            found = true;
            break;
          }
        }
      }
      // If triangle has changed fix deplist
      if (found != true) {
        ep->deps.clear();
        for (int i = 0; i != edgep->deps.size(); i++) {
          if (edgep->deps[i]->ptype->collidable)
            ep->deps.push_back(edgep->deps[i]);
        }
      }
      triple *v1old = v1;
      if (ntnew != NULL) {
        setTriangle(ntnew);
        // Rebuild Edge index:
        int faceind = -1;
        edgeindex = edgeIndexDetect();
      }
      return true;
    }
    triple *v1;
    triple *v2;
    Property *ep;
    NTriangle *nt;
    triple c1, c2, di; // dummy edge coefficients;
    bool ishealthy;
};

class TauEdgeFactory : public CustomDepEnergyFactory {
  public:
    CustomDepEnergy *GetNew(void *parent) {
    EdgeEnerProperty *ep = (EdgeEnerProperty *)parent;
    TauEdge *bne = new TauEdge;
    bne->setParent(ep);
    bne->edgeindex = ep->state;
    for (int i = 0; i != ep->deps.size(); i++) {
      if (ep->deps[i]->ptype->collidable) {
        bne->setTriangle((NTriangle *)ep->deps[i]->content);
        break;
      }
    }

    bne->next = NULL;
    return bne;
  }
};
