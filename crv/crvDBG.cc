/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "crvDBG.h"
#include "crvAdapt.h"
#include "crv.h"
#include <gmi.h>
#include <gmi_null.h>
#include <PCU.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfConvert.h>
#include <apfDynamicArray.h>
#include <pcu_util.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cfloat>
#include <cassert>
#include <stdarg.h>
// ==== includes for safe_dir =====
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h>  /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h>     /*for checking the error from mkdir */
// ================================

namespace crv_dbg {

void safe_mkdir(const char* path);

void createCavityMeshIndividualTets(
    ma::Adapt* a,
    ma::EntityArray& ents,
    const char* prefix)
{
  /* std::stringstream ss; */
  /* ss<<prefix<<"_AllTets"; */
  /* const char* outPrefix = ss.str().c_str(); */
  /* safe_mkdir(outPrefix); */

  apf::Mesh2* m = a->mesh;
  int nent = ents.getSize();
  int dim = m->getDimension();

  int count = 0;
  for (int i = 0; i < nent; i++) {
    if (m->getType(ents[i]) != apf::Mesh::TET) continue;
    count++;
  }
  ma::EntityArray tets;
  tets.setSize(count);
  count = 0;
  for (int i = 0; i < nent; i++) {
    if (m->getType(ents[i]) != apf::Mesh::TET) continue;
    tets[count] = ents[i];
    count++;
  }

  for (int ii = 0; ii < count; ii++) {
    gmi_model* g = gmi_load(".null");
    apf::Mesh2* outMesh = apf::makeEmptyMdsMesh(g, dim, false);

    // Verts
    apf::MeshEntity* vs[12];
    apf::MeshEntity* newVs[12];
    int nv = m->getDownward(tets[ii], 0, vs);
    for(int i = 0; i < nv; i++) {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[i], 0, p);
      newVs[i] = outMesh->createVertex(0, p, param);
    }

    // Edges
    apf::MeshEntity* es[12];
    apf::MeshEntity* newEs[12];
    int ne = m->getDownward(tets[ii], 1, es);
    for(int i = 0; i < ne; ++i) {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* new_evs[2];
      m->getDownward(es[i], 0, evs);
      for (int j = 0; j < 2; j++) {
	new_evs[j] = newVs[apf::findIn(vs, nv, evs[j])];
      }

      newEs[i] = outMesh->createEntity(apf::Mesh::EDGE, 0, new_evs);
    }
    // Faces
    apf::MeshEntity* fs[12];
    apf::MeshEntity* newFs[12];
    int nf = m->getDownward(tets[ii], 2, fs); 
    for(int i = 0; i < nf; ++i) 
    {
      apf::MeshEntity* fes[3];
      apf::MeshEntity* new_fes[3];
      m->getDownward(fs[i], 1, fes);
      for (int j = 0; j < 3; j++) {
	new_fes[j] = newEs[apf::findIn(es, ne, fes[j])];
      }
      newFs[i] = outMesh->createEntity(apf::Mesh::TRIANGLE, 0, new_fes);
    }

    // Regions
    apf::MeshEntity* tet;
    tet = outMesh->createEntity(apf::Mesh::TET, 0, newFs);

    apf::MeshEntity* newEnt;
    newEnt = tet;

    PCU_ALWAYS_ASSERT(m->getType(tets[ii]) == outMesh->getType(newEnt));
    outMesh->acceptChanges();

    apf::changeMeshShape(outMesh,
      crv::getBezier(m->getShape()->getOrder()), true);
    outMesh->acceptChanges();
    for (int d = 1; d <= dim; d++)
    {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* eds[12];
      int counter = m->getDownward(tets[ii], d, eds);
      apf::MeshEntity* new_eds[12];
      outMesh->getDownward(newEnt, d, new_eds);
      int non = outMesh->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for(int n = 0; n < counter; ++n) {
	for(int i = 0; i < non; ++i) {
	  apf::Vector3 p;
	  m->getPoint(eds[n], i, p);
	  outMesh->setPoint(new_eds[n], i, p);
	}
      }
    }
    outMesh->acceptChanges();

    std::stringstream ss;
    ss << prefix<< "TET"<<ii;
    /* ss << outPrefix<<"/"<<"TET"<<ii; */
    crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TET, 20, ss.str().c_str());
    crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TRIANGLE, 20, ss.str().c_str());
    crv::writeCurvedWireFrame(outMesh, 50, ss.str().c_str());
    ss.str("");

    outMesh->destroyNative();
    apf::destroyMesh(outMesh);
  }
}

void createCavityMesh(ma::Adapt* a,
    ma::EntityArray& ents,
    const char* prefix)
{
  apf::Mesh2* m = a->mesh;
  int nent = ents.getSize();
  int dim = m->getDimension();

  gmi_model* g = gmi_load(".null");
  apf::Mesh2* outMesh = apf::makeEmptyMdsMesh(g, dim, false);

  int count = 0;
  for (int i = 0; i < nent; i++) {
    if (m->getType(ents[i]) != apf::Mesh::TET) continue;
    count++;
  }
  ma::EntityArray tets;
  ma::EntityArray newEnt;
  tets.setSize(count);
  newEnt.setSize(count);
  count = 0;
  for (int i = 0; i < nent; i++) {
    if (m->getType(ents[i]) != apf::Mesh::TET) continue;
    tets[count] = ents[i];
    count++;
  }

  for (int ii = 0; ii < count; ii++) {
    // Verts
    apf::MeshEntity* vs[12];
    apf::MeshEntity* newVs[12];
    int nv = m->getDownward(tets[ii], 0, vs);
    for(int i = 0; i < nv; ++i) {
      apf::Vector3 p;
      apf::Vector3 param(0., 0., 0.);
      m->getPoint(vs[i], 0, p);
      newVs[i] = outMesh->createVertex(0, p, param);
    }

    // Edges
    apf::MeshEntity* es[12];
    apf::MeshEntity* newEs[12];
    int ne = m->getDownward(tets[ii], 1, es);
    for(int i = 0; i < ne; ++i) {
      apf::MeshEntity* evs[2];
      apf::MeshEntity* new_evs[2];
      m->getDownward(es[i], 0, evs);
      for (int j = 0; j < 2; j++) {
        new_evs[j] = newVs[apf::findIn(vs, nv, evs[j])];
      }
      newEs[i] = outMesh->createEntity(apf::Mesh::EDGE, 0, new_evs);
    }

    // Faces
    apf::MeshEntity* fs[12];
    apf::MeshEntity* newFs[12];
    int nf = m->getDownward(tets[ii], 2, fs);
    for(int i = 0; i < nf; ++i) {
      apf::MeshEntity* fes[3];
      apf::MeshEntity* new_fes[3];
      m->getDownward(fs[i], 1, fes);
      for (int j = 0; j < 3; j++) {
        new_fes[j] = newEs[apf::findIn(es, ne, fes[j])];
      }
      newFs[i] = outMesh->createEntity(apf::Mesh::TRIANGLE, 0, new_fes);
    }

    // Regions
    apf::MeshEntity* tet;
    tet = outMesh->createEntity(apf::Mesh::TET, 0, newFs);
    newEnt[ii] = tet;

    PCU_ALWAYS_ASSERT(m->getType(tets[ii]) == outMesh->getType(newEnt[ii]));
    outMesh->acceptChanges();
  }

  apf::changeMeshShape(outMesh,
      crv::getBezier(m->getShape()->getOrder()), true);
  outMesh->acceptChanges();

  for (int ii = 0; ii < count; ii++) {
    for (int d = 1; d <= dim; d++) {
      if (!m->getShape()->hasNodesIn(d)) continue;
      apf::MeshEntity* eds[12];
      int counter = m->getDownward(tets[ii], d, eds);
      apf::MeshEntity* new_eds[12];
      outMesh->getDownward(newEnt[ii], d, new_eds);
      int non = outMesh->getShape()->countNodesOn(apf::Mesh::simplexTypes[d]);
      for(int n = 0; n < counter; ++n) {
        for(int i = 0; i < non; ++i) {
        apf::Vector3 p;
        m->getPoint(eds[n], i, p);
        outMesh->setPoint(new_eds[n], i, p);
        }
      }
    }
    outMesh->acceptChanges();
  }

  std::stringstream ss;
  ss << prefix;
  crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TET, 8, ss.str().c_str());
  crv::writeCurvedVtuFiles(outMesh, apf::Mesh::TRIANGLE, 16, ss.str().c_str());
  crv::writeCurvedWireFrame(outMesh, 16, ss.str().c_str());

  outMesh->destroyNative();
  apf::destroyMesh(outMesh);
}

void createCavityMesh(ma::Adapt* a,
    ma::EntitySet& ents,
    const char* prefix)
{
  ma::EntityArray entArray;
  entArray.setSize(ents.size());
  int count = 0;
  APF_ITERATE(ma::EntitySet,ents,it)
    entArray[count++] = *it;
  createCavityMesh(a, entArray, prefix);
}

}
