/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_EMBED_SURFACE_H
#define MA_EMBED_SURFACE_H

#include <apf.h>
#include <pcu_util.h>
#include "maMesh.h"
#include "maTables.h"
#include "maRefine.h"

namespace ma {

class Adapt;

class Surface
{
  public:
    virtual ~Surface() {};
    virtual int getEdgeIntersectionXi(Entity* e, Vector& xi, Tag* t) = 0;
    virtual apf::Field* getField() = 0;
    /* virtual bool isIntersecting(Entity* e) = 0; */
};

class LevelSetSurface : public Surface
{
  public:
    LevelSetSurface(Mesh* m, apf::Field* p, double t) :
      mesh(m), phi(p), tol(t)
    {}
    ~LevelSetSurface()
    {}
    virtual int getEdgeIntersectionXi(Entity* e, Vector& xi, Tag* t)
    {
      Entity* dv[2];
      mesh->getDownward(e, 0, dv);
      apf::MeshElement* me = apf::createMeshElement(mesh, e);
      apf::Element* el = apf::createElement(phi, me);

      Vector xi_low(-1.,0.,0.);
      Vector xi_high(1.,0.,0.);

      // Note: this function is called with t = 0 during vertex moves. Hence the
      // reason for checking t.
      // The purpose of this function is to avoid finding the edge intersection if
      // one of the verts is already moved to the surface.
      if (t)
	if (mesh->hasTag(dv[0], t) || mesh->hasTag(dv[1], t))
	  return 0;



      xi = (xi_low + xi_high) * 0.5;

      double f_low, f_high, f;

      f_low = apf::getScalar(el, xi_low);
      f_high = apf::getScalar(el, xi_high);
      f = apf::getScalar(el, xi);

      if (f_low * f_high > 0.) return 0;

      while (std::abs(f) > tol) {
      	if (f*f_low < 0) {
      	  xi_high = xi;
      	  f_high = apf::getScalar(el, xi_high);
	}
	else {
	  xi_low = xi;
	  f_low = apf::getScalar(el, xi_low);
	}
	xi = (xi_low + xi_high) * 0.5;
        f = apf::getScalar(el, xi);
      }

      apf::destroyElement(el);
      apf::destroyMeshElement(me);

      // The following check is to make sure not too many short edges are created.
      // Using a hard-coded threshold value here is OK since we are comparing
      // parametric coordinates.
      // TODO: Can we do better?
      // Updated Note: having this logic makes keeping track of entities on the
      // surface much more difficult. So for now let's have this commented and try to
      // removed short edges with collapse operations, instead
      /* Entity* dv[2]; */
      /* mesh->getDownward(e, 0, dv); */
      /* int onSruface = 1; */
      /* if (std::abs(xi[0] - 1.) < 0.05) { */
      /* 	if (t) */
      /* 	  mesh->setIntTag(dv[1], t, &onSruface); */
      /* 	return 0; */
      /* } */
      /* if (std::abs(xi[0] + 1.) < 0.05) { */
      /* 	if (t) */
	  /* mesh->setIntTag(dv[0], t, &onSruface); */
      /* 	return 0; */
      /* } */
      return 1;
    }
    virtual apf::Field* getField()
    {
      return phi;
    }
    /* virtual bool isIntersecting(Entity* e) */
    /* { */
    /*   Vector xi; */
    /*   int type = mesh->getType(e); */
    /*   if (type == apf::Mesh::VERTEX) */
    /*   	return false; */
    /*   else if (type == apf::Mesh::EDGE) { */
    /*   	int n = getEdgeIntersectionXi(e, xi, 0); */
    /*   	return n > 0; */
    /*   } */
    /*   else if (type == apf::Mesh::TRIANGLE || */ 
    /*   	       type == apf::Mesh::TET) { */
    /*   	Entity* ds[12]; */
    /*   	int nd = mesh->getDownward(e, 1, ds); */
    /*   	for(int i = 0; i < nd; i++) { */
    /*   	  int n = getEdgeIntersectionXi(ds[i], xi, 0); */
    /*   	  if (n > 0) return true; */
	/* } */
    /*   } */
    /*   else */
    /*   	PCU_ALWAYS_ASSERT(0); */
    /*   return false; */
    /* } */
    apf::Field* phi;
    double tol;
    Mesh* mesh;
};

long markIntersectionEdgesToSplit(Adapt* a, Surface* s, Tag* ei, Tag* vc);

}

#endif
