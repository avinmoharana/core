/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maRefine.h"
#include "maEmbedSurface.h"
#include "maTemplates.h"
#include "maAdapt.h"
#include "maMesh.h"
#include "maTables.h"
#include "maMatch.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"
#include "maSnap.h"
#include "maLayer.h"
#include <apf.h>
#include <pcu_util.h>

namespace ma {

long markIntersectingEdgesToSplit(Adapt* a, Surface* s, Tag* ei)
{
  int dimension = 1;
  int trueFlag = SPLIT;
  int setFalseFlag = NEED_NOT_SPLIT;
  int allFalseFlags = DONT_SPLIT | NEED_NOT_SPLIT;

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
    Vector xi;
    int n = s->getEdgeIntersectionXi(e, xi);
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
  return PCU_Add_Long(count);
}

void embedSurface(Input* in, apf::Field* phi)
{
  double t0 = PCU_Time();
  validateInput(in);
  Adapt* a = new Adapt(in);
  Refine* r  = a->refine;
  Surface* s = new LevelSetSurface(a->mesh, phi, 1.e-12);
  long count = markIntersectingEdgesToSplit(a, s, r->intersectionTag);
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
  delete s;
  delete a;
  delete in;
  double t1 = PCU_Time();
  print("embedSurface refined %li edges in %f seconds",count,t1-t0);
}

}
