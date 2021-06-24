
/*******************************************************************************

Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#ifndef MA_SOLUTIONTRANSFERHELPER_H
#define MA_SOLUTIONTRANSFERHELPER_H

#include <apf.h>
#include "maMesh.h"

#include "maAffine.h"
#include "maMap.h"

namespace ma {

int getBestElement(
    apf::Mesh* mesh,
    int n,
    apf::Element** elems,
    Affine* elemInvMaps,
    Vector const& point,
    Vector& bestXi);
int getMinimumDimension(apf::FieldShape* s);
void transfer(
    apf::Field* field,
    double *value,
    int n, // size of the cavity
    Entity** cavity,
    EntityArray& newEntities);
}

#endif
