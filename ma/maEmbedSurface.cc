/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maRefine.h"
#include "maCoarsen.h"
#include "maEmbedSurface.h"
#include "maTemplates.h"
#include "maAdapt.h"
#include "maMesh.h"
#include "maOperator.h"

#include "maTables.h"
#include "maMatch.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"
#include "maSnap.h"
#include "maLayer.h"
#include <apf.h>
#include <apfShape.h>
#include <pcu_util.h>

namespace ma {

static void updateOnSurfaceTag(Mesh* m, Tag* t)
{
  Entity* e;
  Iterator* it;

  int onSurface = 1;

  it = m->begin(1);
  while ( (e = m->iterate(it)) ) {
    Entity* dv[2];
    m->getDownward(e, 0, dv);
    if (m->hasTag(dv[0], t) && m->hasTag(dv[1], t))
      m->setIntTag(e, t, &onSurface);
  }
  m->end(it);

  it = m->begin(2);
  while ( (e = m->iterate(it)) ) {
    Entity* dv[3];
    m->getDownward(e, 1, dv);
    if (m->hasTag(dv[0], t) && m->hasTag(dv[1], t) && m->hasTag(dv[2], t))
      m->setIntTag(e, t, &onSurface);
  }
  m->end(it);

  it = m->begin(3);
  while ( ( e = m->iterate(it)) ) {
    Entity* dv[4];
    m->getDownward(e, 2, dv);
    if (m->hasTag(dv[0], t) && m->hasTag(dv[1], t) && m->hasTag(dv[2], t) &&
          m->hasTag(dv[3], t))
      printf("WARNING \n");
      /* PCU_ALWAYS_ASSERT(0); */
  }
  m->end(it);
}

long markVertsToMove(Adapt* a, Surface* s, Tag* t, double xi_min = 0.05)
{
  int trueFlag = SPLIT;
  int FalseFlag = DONT_SPLIT;

  Entity* e;
  long count = 0;
  Mesh* m = a->mesh;
  int dim = m->getDimension();
  Iterator* it = m->begin(1);
  while ((e = m->iterate(it)))
  {
    Vector xi;
    int n = s->getEdgeIntersectionXi(e, xi, 0);
    if (n > 0)
    {
      int edgeModelType = m->getModelType(m->toModel(e));
      apf::MeshElement* me = apf::createMeshElement(m, e);
      Vector p;
      double pd[3];
      Entity* dv[2];
      m->getDownward(e, 0, dv);
      if (std::abs(xi[0] - 1.) < xi_min &&
      	  !getFlag(a, dv[1], trueFlag) &&
      	  edgeModelType == m->getModelType(m->toModel(dv[1]))) {
      	setFlag(a, dv[1], trueFlag);
	apf::mapLocalToGlobal(me, xi, p);
	p.toArray(pd);
	m->setDoubleTag(dv[1], t, pd);
	if (a->mesh->isOwned(dv[1]))
	  ++count;
      }
      if (std::abs(xi[0] + 1.) < xi_min &&
      	  !getFlag(a, dv[0], trueFlag) &&
      	  edgeModelType == m->getModelType(m->toModel(dv[0]))) {
      	setFlag(a, dv[0], trueFlag);
	apf::mapLocalToGlobal(me, xi, p);
	p.toArray(pd);
	m->setDoubleTag(dv[0], t, pd);
	if (a->mesh->isOwned(dv[0]))
	  ++count;
      }
      apf::destroyMeshElement(me);
    } 
  }
  m->end(it);

  /* updateOnSurfaceTag(m, vc); */

  return PCU_Add_Long(count);
}

class VertexMover : public Operator
{
  public:
    VertexMover(Adapt* a, Surface* s, Tag* t, Tag* on):
      adapter(a), surface(s), tag(t), onSurface(on)
    {
      successCount = 0;
      cavity.init(a);
    }
    int getTargetDimension()
    {
      return 0;
    }
    bool shouldApply(Entity* e)
    {
      if (!getFlag(adapter, e, SPLIT))
      	return false;
      vert = e;
      return true;
    }
    bool requestLocality(apf::CavityOp* o)
    {
      if (!o->requestLocality(&vert, 1))
      	return false;
      else
      	return true;
    }
    void apply()
    {
      Mesh* m = adapter->mesh;
      Model* c = m->toModel(vert);
      int dim = m->getDimension();
      double x[3];
      m->getDoubleTag(vert, tag, x);
      Vector point;
      point.fromArray(x);
      Vector param(0., 0., 0.); // have to figure this out.
      Entity* newVert = buildVertex(adapter, c, point, param);
      // also solution transfer and size field transfer onVertex has to be taken care of.
      // easy solution: level_set value is expected to be 0 so set that for the vert
      apf::setScalar(surface->getField(), newVert, 0, 0.);
      
      // get old elements
      EntityArray oldElements;
      EntityArray newElements;
      apf::Adjacent adj;
      m->getAdjacent(vert, dim, adj);
      oldElements.setSize(adj.getSize());
      newElements.setSize(adj.getSize());
      for (size_t i = 0; i < adj.getSize(); ++i)
      	oldElements[i] = adj[i];

      cavity.beforeBuilding();
      for (size_t i = 0; i < adj.getSize(); ++i)
      	newElements[i] = rebuildElement(m, oldElements[i], vert, newVert,
      	    adapter->buildCallback);
      cavity.afterBuilding();

      if (cavity.shouldFit)
      	cavity.fit(oldElements);
      // check validity
      bool isCavityValid = true;

      if (isCavityValid) {
      	cavity.transfer(oldElements);
      	// 1. destroy old elements
      	// 2. destroy vert
      	// 3. set on surface flag for the newVert
      	// 4. *clear the split flag if need be
      	for (size_t i = 0; i < oldElements.getSize(); ++i)
      	  destroyElement(adapter, oldElements[i]);
      	/* m->destroy(vert); */
        int onSurf = 1;
        m->setIntTag(newVert, onSurface, &onSurf);
      	successCount++;
      }
      else {
      	// 1. destroy new elements
      	// 2. destroy new vertex
      	// 3. clear the split flag on the old vert
      }
    }
    int successCount;
  private:
    Adapt* adapter;
    Surface* surface;
    Tag* tag;
    Tag* onSurface;
    Entity* vert;
    Cavity cavity;
};

static void convertTagToField(Mesh* m, Tag* t, const char* name)
{
  apf::Field* f = apf::createField(m, name, apf::VECTOR, apf::getLagrange(1));
  Entity* v;
  Iterator* it = m->begin(0);

  while ( (v = m->iterate(it)) ) {
    Vector point;
    m->getPoint(v, 0, point);
    Vector newPoint;
    if (m->hasTag(v, t)) {
      double p[3];
      m->getDoubleTag(v, t, p);
      newPoint.fromArray(p);
    }
    else {
      newPoint = point;
    }

    apf::setVector(f, v, 0, newPoint-point);
  }
  m->end(it);
}

static void moveToSurface(Adapt* a, Surface* s, Tag* onSurfaceTag)
{
  double t0 = PCU_Time();
  // Note: for now let's use the SPLIT flag on verts that need to be moved.
  // TODO: we may need to come up with a better flag. In that case we have to
  // update markVertsToMove
  Mesh* m = a->mesh;
  Tag* moveToTag = m->createDoubleTag("embed_vertex_move", 3);
  long count = markVertsToMove(a, s, moveToTag, 0.1);
  if ( ! count) {
    return;
  }
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a,0,SPLIT));

  /* convertTagToField(m, moveToTag, "move_field"); */
  /* apf::writeVtkFiles("mesh_with_move_tag", m); */

  /* printf("checkpoint 01\n"); */
  VertexMover mover(a, s, moveToTag, onSurfaceTag);
  applyOperator(a, &mover);

  Entity* e;
  Iterator* it;

  for (int d = 0; d <= 3; d++) {
    printf("checking entities of dim %d\n", d);
    it = m->begin(d);
    int count = 0;
    while ( (e = m->iterate(it)) ) {
      if (getFlag(a, e, SPLIT))
      	printf("%d's entity has the split flag \n", count);
      count++;
    }
    m->end(it);
  }

  updateOnSurfaceTag(m, onSurfaceTag);
  // remove the double tag
  m->destroyTag(moveToTag);
  double t1 = PCU_Time();
  print("embedSurface moved %d/%li verts in %f seconds", mover.successCount, count,t1-t0);
}



long markIntersectingEdgesToSplit(Adapt* a, Surface* s, Tag* ei, Tag* vc)
{
  int dimension = 1;
  int trueFlag = SPLIT;
  int setFalseFlag = NEED_NOT_SPLIT;
  int allFalseFlags = DONT_SPLIT | NEED_NOT_SPLIT;

  if (!allFalseFlags) allFalseFlags = setFalseFlag;
  Entity* e;
  long count = 0;
  int count2 = 0;
  Mesh* m = a->mesh;
  Iterator* it = m->begin(dimension);
  while ((e = m->iterate(it)))
  {
    count2++;
    PCU_ALWAYS_ASSERT( ! getFlag(a,e,trueFlag));
    /* this skip conditional is powerful: it affords us a
       3X speedup of the entire adaptation in some cases */
    if (allFalseFlags & getFlags(a,e))
      continue;
    Vector xi;
    printf("count2 is %d\n", count2);
    if (m->hasTag(e, vc))
      printf("edge has the on surface tag\n");
    else
      printf("edge does not have the on surface tag\n");
    int n = s->getEdgeIntersectionXi(e, xi, vc);
    printf("after count2 is %d\n", count2);
    if (n > 0)
    {   
      setFlag(a,e,trueFlag);
      m->setDoubleTag(e, ei, &xi[0]);
      if (a->mesh->isOwned(e))
        ++count;
    }   
    else
      setFlag(a,e,setFalseFlag);
  }
  m->end(it);

  /* updateOnSurfaceTag(m, vc); */

  return PCU_Add_Long(count);
}

long markSurfaceEdgesToCollapse(Adapt* a, Tag* t, double l)
{
  int dimension = 1;
  int trueFlag = COLLAPSE;
  int setFalseFlag = NEED_NOT_COLLAPSE;
  int allFalseFlags = DONT_COLLAPSE | NEED_NOT_COLLAPSE;

  if (!allFalseFlags) allFalseFlags = setFalseFlag;
  Entity* e;
  long count = 0;
  Mesh* m = a->mesh;
  Iterator* it = m->begin(dimension);
  while ((e = m->iterate(it)))
  {
    PCU_ALWAYS_ASSERT( ! getFlag(a,e,trueFlag));
    /* this skip conditional is powerful: it affords us a
       3X speedup of the entire adaptation in some cases */
    if (allFalseFlags & getFlags(a,e))
      continue;
    if (!m->hasTag(e, t))
      continue;
    Entity* dv[2];
    m->getDownward(e, 0, dv);
    Vector p0;
    Vector p1;
    m->getPoint(dv[0], 0, p0);
    m->getPoint(dv[1], 0, p1);

    double length = (p1 - p0).getLength();

    if (length < l) {
      setFlag(a,e,trueFlag);
      if (m->isOwned(e))
      	count++;
    }
    else
      setFlag(a,e,setFalseFlag);
  }
  m->end(it);

  /* updateOnSurfaceTag(m, vc); */

  return PCU_Add_Long(count);
}

static Tag* transferOnSurfaceTag(Mesh* m, const char* name)
{
  Tag* inTag = m->findTag("ma_on_embed_surface");
  PCU_ALWAYS_ASSERT(inTag);

  Tag* outTag = m->createIntTag(name, 1);

  int onSurf = 1;
  for (int d = 0; d < m->getDimension(); d++)
  {
    Entity* e;
    Iterator* it = m->begin(d);
    while ( (e = m->iterate(it)) ) {
      if (m->hasTag(e, inTag))
        m->setIntTag(e, outTag, &onSurf);
    }
    m->end(it);
  }

  return outTag;
}




static void refineOnSurface(Adapt* a, Surface* s, Refine* r)
{
  double t0 = PCU_Time();
  printf("checkpoint refine 01\n");
  long count = markIntersectingEdgesToSplit(a, s, r->intersectionTag, r->onSurfaceTag);
  printf("checkpoint refine 02\n");
  if ( ! count) {
    return;
  }
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a,1,SPLIT));

  resetCollection(r);
  collectForTransfer(r);
  collectForMatching(r);
  addAllMarkedEdges(r);
  splitElements(r);
  processNewElements(r);
  destroySplitElements(r);
  forgetNewEntities(r);
  updateOnSurfaceTag(a->mesh, r->onSurfaceTag);

  double t1 = PCU_Time();
  print("embedSurface refined %li edges in %f seconds",count,t1-t0);
}


static void coarsenOnSurface(Adapt* a, Tag* t, double length, int iter = 3)
{
  for (int i = 0; i < iter; i++) {
    double t0 = PCU_Time();
    long markedCount = 0;
    long successCount = 0;
    markedCount = markSurfaceEdgesToCollapse(a, t, length);
    for (int d = 1; d <= a->mesh->getDimension(); d++) {
      ma::checkAllEdgeCollapses(a, d);
      findIndependentSet(a);
      successCount += collapseAllEdges(a, d);
    }
    updateOnSurfaceTag(a->mesh, t);
    double t1 = PCU_Time();
    markedCount = PCU_Add_Long(markedCount);
    successCount = PCU_Add_Long(successCount);
    print("embedSurface coarsened %li/%li edges in %f seconds", markedCount, successCount, t1-t0);
  }
}

void embedSurface(Input* in, apf::Field* phi)
{
  //// general setups
  bool shouldMove = true;
  bool shouldRefine = true;
  bool shouldCoarsen = true;

  validateInput(in);
  Adapt* a = new Adapt(in);
  Refine* r  = a->refine;
  Surface* s = new LevelSetSurface(a->mesh, phi, 1.e-12);

  //// moving close verts to the surface
  if (shouldMove) {
    moveToSurface(a, s, r->onSurfaceTag);
  }
  // for debugging
  /* Tag* t1 = transferOnSurfaceTag(a->mesh, "on_surface_after_moves"); */

  //// surface refinements
  if (shouldRefine) {
    refineOnSurface(a, s, r);
  }
  // for debugging
  /* Tag* t2 = transferOnSurfaceTag(a->mesh, "on_surface_after_splits"); */

  //// surface collapses
  if (shouldCoarsen) {
    double length = 0.4; // user specified
    double halfLength = length/2.;
    /* in->shouldForceAdaptation = true; */
    coarsenOnSurface(a, r->onSurfaceTag, 0.05*halfLength);
    coarsenOnSurface(a, r->onSurfaceTag,     halfLength);
  }
  printf("checkpoint hh 01\n");
  // for debugging
  /* Tag* t3 = transferOnSurfaceTag(a->mesh, "on_surface_after_collapses"); */
  /* printf("checkpoint 02\n"); */

  //// cleanups
  delete s;
  delete a;
  delete in;
  printf("checkpoint hh 02\n");

  /* a->mesh->destroyTag(t1); */
  /* printf("checkpoint 03a\n"); */
  /* a->mesh->destroyTag(t2); */
  /* printf("checkpoint 03b\n"); */
  /* a->mesh->destroyTag(t3); */
  /* printf("checkpoint 04\n"); */

  a->mesh->removeField(phi);
  apf::destroyField(phi);
  printf("checkpoint hh 03\n");
}

}
