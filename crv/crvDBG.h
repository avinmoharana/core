/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef CRV_DBG_H
#define CRV_DBG_H

#include <apfMesh2.h>
#include <apfShape.h>
#include <apfNumbering.h>
#include <maStats.h>
#include <ma.h>

#include <vector>
#include <assert.h>

namespace crv_dbg {


void createCavityMeshIndividualTets(ma::Adapt* a,
    ma::EntityArray& tets,
    const char* prefix);

void createCavityMesh(ma::Adapt* a,
    ma::EntityArray& tets,
    const char* prefix);

void createCavityMesh(ma::Adapt* a,
    ma::EntitySet& tets,
    const char* prefix);

}
#endif
