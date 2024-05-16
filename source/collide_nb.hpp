#pragma once

#include "collision.hpp"

//#include "config_ng.hpp"
//#include"nagata.hpp"

class Neighbor
{
public:
    Neighbor(){szBounds=4;node_id=-1;}
    virtual bool getbound(vector<triple*> *bounds) =0;
    int szBounds;
    int node_id;
};

class NeighborCollection
{
public:
    virtual bool exposeNeighbors(vector<Neighbor*> &col)=0;
};


class Neighborlist
{
public:
        Neighborlist() {col.clear();dsq.clear();neighbors.clear();}
        void createNearestNeighborList(double radius);       
        vector<indexlist> neighbors;        // list of neighbor objects
        indexlist parent_index;        // number inside parent entity
        indexlist parent_entity;        // number inside parent entity
        vector< vector <double> > dsq; // Distance squared to neighbors for Energy cutoffs
        deque<Neighbor *> col; // Imagine Face/Ion or Face RandomObject Neighbors;
        //deque<Collidable *> colB;
    
};

void Neighborlist::createNearestNeighborList(double radius)
{
    /**
     * Given a radius, this function will creates for each Collidable (face) a list of faces within the neighborhood of this radius.
     *
     * To evaluate the distance we use the coordinates of face vertices and one point inside the curved face
     * (in case the face is highly curved this point can have a significant distace from the vertices)
     *
     * This function is used for collision detection. The collision is than checked only between the neighbors and
     * does not need to be checked for all pairs of faces.
     */

    double r = radius;
 /*   if (radius <= 0)
    {
        r = 5*sqrt(TotalArea()/fv->faces.size());
    }
*/
    int lsize=col.size();
    neighbors.resize(lsize);
    dsq.resize(lsize);
    double radiusSquared = r * r;
  
    // compute mieghbors for the rest of the entities
    int ind1=0,ind2=0;
    for(int i=0;i!=neighbors.size();i++)
    {
    neighbors[i].clear(); dsq[i].clear();
    }
    //  return;
///#pragma omp parallel for
    for (auto i = col.begin(); i!=col.end(); i++)
     {
         col[ind1]->node_id=ind1; // Make Neighbors accesible via properties.
//#pragma omp atomic
         ind2=ind1+1;
//         neighbors[ind2].clear();
         // Vectorize here!
#pragma omp parallel for
      for (auto j = i + 1; j< col.end(); j++)
      {
            vector<triple*> box1(4);
            (*i)->getbound(&box1);
            vector<triple*> box2(4);
            (*j)->getbound(&box2);
            bool nghbrs = false;
            for (int k = 0;k!=box1.size();k++)
            {
                for (int l = 0; l!=box2.size();l++)
                {
                    double dsquared=box1[k]->distanceSquared(*box2[l]);
                   // cout << dsquared << endl;
                    if ( dsquared < radiusSquared)
                    {                                                   
#pragma omp critical
                        {
                     //   cout << dsquared << endl;
                        if(!nghbrs) 
                        {
                            nghbrs = true;

                         neighbors[ind1].push_back(ind2);
                         neighbors[ind2].push_back(ind1);
                         dsq[ind1].push_back(dsquared);
                         dsq[ind2].push_back(dsquared);
                            
                        }
                        else if(nghbrs and dsquared==0)
                        {
                         *(dsq[ind1].end()--)=0;
                         *(dsq[ind2].end()--)=0;
                        }
                        //if(dsquared==0) break;
                        }
                  }
                  if (nghbrs)
                    break;
               }
                  if (nghbrs)
                    break;
            }
#pragma omp atomic
           ind2++;
        }
        ind1++;
        /*if (nearestNeighbors[i].empty())
        {
            cout << "No neighbors found for face with index " << i << ".\n";
        }*/
    }
}

class CollisionDetection
{
 public:
    CollisionDetection() {n=NULL;}
    CollisionDetection(Neighborlist *ni) {n=ni;}
    Neighborlist *n;
    bool CheckCollision(int index)
    {
        bool collided=false;
        vector<triple*> f1(4,0);
        vector<triple*> f2(4,0);
        Neighbor *center=n->col[index];
        center->getbound(&f1); 
    //    cout << "INDEX " << index << endl; 
        for(int i=0;i!=n->neighbors[index].size();i++)
        {
   //         cout << "NB" << n->neighbors[index][i] << endl;
            Neighbor *off=n->col[n->neighbors[index][i]];
        off->getbound(&f2); 

//            off->getbound(f2);
            collided=TetraCollision(f1, f2); //TODO ADD More options for collisions
            if(collided) return collided;
        }
        return false;
    }
    
    bool TetraCollision(vector<triple*> &f1, vector<triple*> &f2)
    {
    
    bool res;
    
 //    mid2=f2->eval(0.667, 0.333);

    for (int i = 0; i < 4; i ++)
    {
        for (int j = 0; j < 4; j ++)
        {
            int i1 = (i + 1) % 4;
            int i2 = (i + 2) % 4;
            int i3 = (i + 3) % 4;
            int j1 = (j + 1) % 4;
            int j2 = (j + 2) % 4;
            int j3 = (j + 3) % 4;

            res = moeller::TriangleIntersects< double [3]>::triangle(
                            f1[i1]->p, f1[i2]->p, f1[i3]->p, f2[j1]->p, f2[j2]->p,f2[j3]->p);
            if (res)
            {
                         return true;
            }
        }
    }

    return false;
   }
    
};
