/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>
#include <lionPrint.h>
#include <maCoarsen.h>
#include <maEdgeSwap.h>
#include <maOperator.h>
#include <maShapeHandler.h>
#include <maShape.h>
#include <pcu_util.h>
#include <iostream>
#include "crv.h"
#include "crvAdapt.h"
#include "crvQuality.h"
#include "crvShape.h"
#include "crvTables.h"
#include "crvShapeFixer.h"
#include "crvOptimizations.h"

/* This is similar to maShape.cc, conceptually, but different enough
 * that some duplicate code makes sense */

namespace crv {

bool isBoundaryEntity(apf::Mesh* m, apf::MeshEntity* e)
{
  return m->getModelType(m->toModel(e)) < m->getDimension();
}
/** \brief checks if any entity has two entities of
 * dimension on the boundary
 * \details this is useful for some shape correction assessments,
 * and in general, curved elements with multiple entities on the boundary
 * are at risk for poor quality, since this strongly constrains
 * their shape */
static bool hasTwoEntitiesOnBoundary(apf::Mesh* m, apf::MeshEntity* e, int dimension)
{
  apf::Downward down;
  int count = 0;
  int nd = m->getDownward(e,dimension,down);
  for (int i = 0; i < nd; ++i){
    if(isBoundaryEntity(m,down[i]))
      ++count;
    if(count == 2)
      return true;
  }
  return false;
}

/* Mark Edges based on the invalidity code the element has been
 * tagged with.
 */

static std::vector<int> getEdgeSequenceFromInvalidVertex(ma::Mesh* mesh, ma::Entity* e, int index)
{
  apf::MeshEntity* edges[6];
  mesh->getDownward(e, 1, edges);

  apf::MeshEntity* f[4];
  int nf = mesh->getDownward(e, 2, f);
  apf::MeshEntity* vf[3];
  apf::MeshEntity* vt[4];
  mesh->getDownward(e, 0, vt);

  std::vector<int> a;
  std::vector<int> b;

  // get the vertex local number for face from the invalid vertex
  for (int i = 0; i < nf; i++) {
    mesh->getDownward(f[i], 0, vf);
    int j = apf::findIn(vf, 3, vt[index]);

    if (j != -1) {
      a.push_back(i);
      b.push_back(j);
    }
  }

  //getJacobianDeterminant on each face at
  double c[3]; //contains jacobian value
  int cc[3];   //face index
  apf::Vector3 xi;
  apf::Matrix3x3 Jac;
  for(size_t k = 0; k < a.size(); k++) {
    if (b[k] == 0) xi = {0, 0, 0};
    else if (b[k] == 1) xi = {1, 0, 0};
    else xi = {0, 1, 0};

    cc[k] = k;

    apf::MeshElement* me = apf::createMeshElement(mesh, f[a[k]]);
    apf::getJacobian(me, xi, Jac);
    c[k] = apf::getJacobianDeterminant(Jac, 2);

    apf::destroyMeshElement(me);
  }

  // sort min to max jacobian determinant
  for (int i = 0; i < 3; i++) {
    for (int j = i-1; j >= 0; --j) {
      double k = c[j+1];
      int kk = cc[j+1];
      if ( k < c[j]) {
      	c[j+1] = c[j];
      	c[j] = k;
      	cc[j+1] = cc[j];
      	cc[j] = kk;
      }
    }
  }

  apf::MeshEntity* ef[3];
  apf::MeshEntity* ve[2];
  std::vector<int> aa;
  //std::vector<int> bb;

  // need to flag all adjacent edges
  // hence only 2 faces would do
  for (int i = 0; i < 2; i++) {
    int nef = mesh->getDownward(f[a[cc[i]]], 1, ef);
    for (int j = 0; j < nef; j++) {
      mesh->getDownward(ef[j], 0, ve);
      if (apf::findIn(ve, 2, vt[index]) != -1) {
      	int k = apf::findIn(edges, 6, ef[j]);
      	if (k != -1) {
      	  if (aa.size() >= 2) {
      	    bool contains = false;
      	    for (std::size_t ii = 0; ii < aa.size(); ii++)
      	      contains = contains || (k == aa[ii]);
      	    if (contains == false)
      	      aa.push_back(k);
	  }
	  else {
	    aa.push_back(k);
      	    //bb.push_back(mesh->getModelType(mesh->toModel(ef[j])));
	  }
	}
      }
    }
  }

  // aa has the index of the edges ordered
  // min to max of adj face jacobian

  return aa;

}

static void sortDesWrtFreqVector(std::vector<int> &f, std::vector<int> &a)
{
  for (std::size_t i = 0; i < f.size(); i++) {
    for (int j = i-1; j >= 0; j--) {
      int ko = f[j+1];
      int ke = a[j+1];
      if (ko > f[j]) {
      	f[j+1] = f[j];
      	a[j+1] = a[j];
      	f[j] = ko;
      	a[j] = ke;
      }
    }
  }
}

static std::vector<int> sortEdgeIndexByFrequency(std::vector<int> &all)
{
  std::vector<int> freq, unqEntries;
  freq.push_back(1);
  unqEntries.push_back(all[0]);

  for (std::size_t i = 1; i < all.size(); i++) {
    bool milila = false;

    for (std::size_t j = 0; j < unqEntries.size(); j++) {
      milila = milila || (all[i] == unqEntries[j]);

      if (all[i] == unqEntries[j]) {
      	freq[j] = freq[j] + 1;
      }
    }
    if (milila == false) {
      unqEntries.push_back(all[i]);
      freq.push_back(1);
    }
  }

  sortDesWrtFreqVector(freq, unqEntries);

  return unqEntries;
}

static std::vector<int> numOfUniqueFaces(std::vector<int> &ai)
{
  std::vector<int> aiEdgeNVtx;
  std::vector<int> aiOnlyFace;
  std::vector<int> aiOnlyEdge;
  std::vector<int> aiOnlyVtx;
  std::vector<int> aiOnlyTet;

  std::vector<int> aiUniqueFace;
  std::vector<int> FFE;

  for (size_t i = 0; i < ai.size(); i++) {
    if (ai[i] < 8) aiOnlyVtx.push_back(ai[i]);
    else if (ai[i] > 7 && ai[i] < 14) aiOnlyEdge.push_back(ai[i]);
    else if (ai[i] > 13 && ai[i] < 20) aiOnlyFace.push_back(ai[i]);
    else aiOnlyTet.push_back(ai[i]);

  }

  int faceToEdgeInd[4][4] = {{-1, 0, 1, 2},
			     {0, -1, 4, 3},
			     {1, 4, -1, 5},
			     {2, 3, 5, -1}};

  int edgeToFaceInd[6][2] = {{0, 1},
			     {0, 2},
			     {0, 3},
			     {1, 3},
			     {1, 2},
			     {2, 3}};

  int edgeToVtx[6][2] = {{0, 1},
			 {1, 2},
			 {0, 2},
			 {0, 3},
			 {1, 3},
			 {2, 3}};

  for (std::size_t i = 0; i < aiOnlyFace.size(); i++) {
    for (std::size_t j = 0; j < aiOnlyFace.size(); j++) {
      if ((i != j) && (j > i)) {
      	for (std::size_t k = 0; k < aiOnlyEdge.size(); k++) {
      	  if (faceToEdgeInd[aiOnlyFace[i]-14][aiOnlyFace[j]-14] + 8 == aiOnlyEdge[k]) {
      	    FFE.push_back(aiOnlyEdge[k]);
      	    FFE.push_back(aiOnlyFace[j]);
      	    FFE.push_back(aiOnlyFace[i]);
      	    break;
	  }
	}
	break;
      }
    }
  }

  bool hasDownVtx, alreadyIn;

  for (size_t j = 0; j < aiOnlyEdge.size(); j++) {
    hasDownVtx = false;
    for (int jj = 0; jj < 2; jj++) {
      //hasDownVtx = false;
      for (size_t ii = 0; ii < aiOnlyVtx.size(); ii++) {
      	hasDownVtx = hasDownVtx || (edgeToVtx[aiOnlyEdge[j]-8][jj] == aiOnlyVtx[ii] - 2);
      }
    }

    if (hasDownVtx == false) {
      for (int i = 0; i < 2; i++) {
      	alreadyIn = false;
      	for (size_t k = 0; k < aiOnlyFace.size(); k++) {
      	  alreadyIn = alreadyIn || (edgeToFaceInd[aiOnlyEdge[j]-8][i] == aiOnlyFace[k]-14);
	}
	if (alreadyIn == false) {
	  aiOnlyFace.push_back(edgeToFaceInd[aiOnlyEdge[j]-8][i] + 14);
	}
      }
    }
    else {
      for (int i = 0; i < 2; i++) {
      	for (size_t k = 0; k < aiOnlyFace.size(); k++) {
      	  if (edgeToFaceInd[aiOnlyEdge[j]-8][i] == aiOnlyFace[k]-14) {
      	    //aiOnlyFace.erase(aiOnlyFace.begin()+k);
	  }
	}
      }
    }
  }

  for (size_t j = 0; j < aiOnlyFace.size(); j++) {
    aiUniqueFace.push_back(aiOnlyFace[j]);
  }

  ai.clear();

  std::vector<int> forFaceMarking;

  // TODO: to be cleaned up
  //[number of unique faces, {unique face index}, {pairs of face face coomon edge invalids}]

  if (aiUniqueFace.size() > 0) {
    forFaceMarking.push_back(aiUniqueFace.size());
    for (std::size_t i = 0; i < aiOnlyVtx.size(); i++)
      ai.push_back(aiOnlyVtx[i]);
    for (std::size_t i = 0; i < aiOnlyEdge.size(); i++)
      ai.push_back(aiOnlyEdge[i]);
    for (std::size_t i = 0; i < aiOnlyFace.size(); i++)
      ai.push_back(aiOnlyFace[i]);
    for (std::size_t i = 0; i < aiOnlyTet.size(); i++)
      ai.push_back(aiOnlyTet[i]);
    for (std::size_t i = 0; i < aiUniqueFace.size(); i++) {
      forFaceMarking.push_back(aiUniqueFace[i]);
    }
    for (std::size_t i = 0; i < FFE.size(); i++)
      forFaceMarking.push_back(FFE[i]);

  }
  else if (aiUniqueFace.size() == 0) {
    forFaceMarking.push_back(0);
    if (FFE.size() > 0) {
      for (std::size_t i = 0; i < FFE.size(); i++)
      	forFaceMarking.push_back(FFE[i]);
    }
    for (std::size_t i = 0; i < aiOnlyVtx.size(); i++)
      ai.push_back(aiOnlyVtx[i]);
    for (std::size_t i = 0; i < aiOnlyEdge.size(); i++)
      ai.push_back(aiOnlyEdge[i]);
    for (std::size_t i = 0; i < aiOnlyFace.size(); i++)
      ai.push_back(aiOnlyFace[i]);
  }
  return forFaceMarking;
}

static int markAllEdges(ma::Mesh* m, ma::Entity* e,
    std::vector<int> ai, ma::Entity* edges[6])
{
  std::vector<int> bb;
  apf::MeshEntity* edf[3];
  apf::MeshEntity* faces[4];
  ma::Downward ed;
  m->getDownward(e,1,ed);

 //std::vector<int> ai = inspectInvalidies(aiall);
  std::vector<int> faceInvalid = numOfUniqueFaces(ai);

  if (ai.size() == 0) return 0;
  if (ai.size() == (std::size_t) faceInvalid[0]) return 0;

  for (size_t ii = 0; ii < ai.size(); ii++) {

    int dim = (ai[ii]-2)/6;
    int index = (ai[ii]-2) % 6;

    switch (dim) {
      case 0:
      {
      	//ma::Downward ed;
      	//m->getDownward(e,1,ed);

      	std::vector<int> aa = getEdgeSequenceFromInvalidVertex(m, e, index);
      	PCU_ALWAYS_ASSERT(index < 4);
      	for (size_t i = 0; i < aa.size(); i++)
      	  bb.push_back(aa[i]);

	break;
      }

      case 1:
      {
      	//ma::Downward ed;
        //m->getDownward(e,1,ed);
        bb.push_back(index);
        bb.push_back(index);
        break;
      }
      //break;
      case 2:
      {
        // if we have an invalid face, operate on its edges
        //ma::Downward edf, faces;
        m->getDownward(e,2,faces);
        m->getDownward(faces[index],1,edf);

	// TODO : This loop seems to be doing nothing
        for (int i = 0; i < 3; i++) {
          int j = apf::findIn(ed, 6, edf[i]);
          j++;
          //if (j != -1)
            //bb.push_back(j);
	}
        break;
      }
      //break;
      case 3:
      {
        m->getDownward(e,1,edges);
        return 6;
      }
      default:
        fail("invalid quality tag in markEdges\n");
        break;
      }
  }

  int n = 0;
  //std::vector<int> allinvEdges = sortEdgeIndexByType(m, e, bb);
  std::vector<int> allinvEdges = sortEdgeIndexByFrequency(bb);

  for (size_t i = 0; i < allinvEdges.size(); i++) {
    if (m->getModelType(m->toModel(ed[allinvEdges[i]])) != 1) {
      n++;
      edges[n-1] = ed[allinvEdges[i]];
    }
  }

  return n;
}

static int markEdges(ma::Mesh* m, ma::Entity* e, int tag,
    ma::Entity* edges[6])
{
  if ( tag <= 1 ) // if its valid, or not checked, don't worry about it
    return 0;
  int dim = (tag-2)/6;
  int index = (tag-2) % 6;
  int n = 0;
  int md = m->getDimension();
  switch (dim) {
    case 0:
    {
      // if we have an invalid vertex, operate on its edges
      ma::Downward ed;
      m->getDownward(e,1,ed);
      n = md;

      if(md == 2){
        edges[0] = ed[index];
        edges[1] = ed[(index+2) % 3];
      } else {
        PCU_ALWAYS_ASSERT(index < 4);
        edges[0] = ed[tetVertEdges[index][0]];
        edges[1] = ed[tetVertEdges[index][1]];
        edges[2] = ed[tetVertEdges[index][2]];
      }
    }
    break;
    case 1:
    {
      // if we have a single invalid edge, operate on it
      ma::Downward ed;
      m->getDownward(e,1,ed);
      edges[0] = ed[index];
      n = 1;
      //break;
    }
    break;
    case 2:
    {
      // if we have an invalid face, operate on its edges
      ma::Downward ed, faces;
      m->getDownward(e,2,faces);
      m->getDownward(faces[index],1,ed);
      n = 3;
      edges[0] = ed[0];
      edges[1] = ed[1];
      edges[2] = ed[2];
      //break;
    }
    break;
    case 3:
    {
      m->getDownward(e,1,edges);
      n = 6;
      //break;
    }
    break;
    default:
      fail("invalid quality tag in markEdges\n");
      break;
  }

  return n;
}

static int markUniqueFaces(ma::Mesh* m, ma::Entity* e, std::vector<int> ai,
    ma::Entity* faces[4])
{

  if (ai.size() == 0) return 0;
  std::vector<int> faceInvalid = numOfUniqueFaces(ai);
  int n = 0;

  ma::Downward fc;
  m->getDownward(e, 2, fc);
  ma::Downward ed;
  m->getDownward(e, 1, ed);
  if (faceInvalid[0] > 0) {
    for (int i = 1; i < faceInvalid[0]+1; i++) {
      int index = (faceInvalid[i] -2) % 6;

      if (m->getModelType(m->toModel(fc[index])) == 3) {
      	faces[n] = fc[index];
      	n++;
      }
    }
  }
  return n;
}

static apf::Vector3 getTriOrientation(apf::Mesh2* m, apf::MeshEntity* tri)
{
  PCU_ALWAYS_ASSERT(m->getType(tri) == apf::Mesh::TRIANGLE);

  apf::MeshEntity* vs[3];
  m->getDownward(tri, 0, vs);
  apf::Vector3 ps[3];
  for (int i = 0; i < 3; i++) {
    m->getPoint(vs[i], 0, ps[i]);
  }
  apf::Vector3 e0 = ps[1] - ps[0];
  apf::Vector3 e1 = ps[2] - ps[1];
  apf::Vector3 n0 = apf::cross(e0, e1);
  return n0/n0.getLength();
}

static int getEdgeIdInSimplex(apf::Mesh2* m, apf::MeshEntity* tri, apf::MeshEntity* edge)
{
  PCU_ALWAYS_ASSERT(m->getType(tri) == apf::Mesh::TRIANGLE);
  PCU_ALWAYS_ASSERT(m->getType(edge) == apf::Mesh::EDGE);
  apf::MeshEntity* es[3];
  m->getDownward(tri, 1, es);
  return apf::findIn(es, 3, edge);
}

class EdgeSwapper : public ma::Operator
{
public:
  EdgeSwapper(Adapt* a)
  {
    adapter = a;
    mesh = a->mesh;
    edges[0] = edges[1] = edges[2] = edges[3] = edges[4] = edges[5] = 0;
    simplex = 0;
    edgeSwap = ma::makeEdgeSwap(a);
    md = mesh->getDimension();
    ne = ns = 0;
  }
  virtual ~EdgeSwapper()
  {
    delete edgeSwap;
  }
  virtual int getTargetDimension() {return md;}
  virtual bool shouldApply(ma::Entity* e)
  {
    /*
    std::vector<int> ai = crv::getAllInvalidities(mesh, e);
    mesh->getDownward(e, 1, edges);
    int niv = ai.size();
    ne = 0;
    if (niv != 0) {
      ne = markAllEdges(mesh, e, ai, edges);
      simplex = e;
    }
    return (ne > 0);
    */
    int tag = crv::getTag(adapter,e);
    ne = markEdges(mesh,e,tag,edges);

    simplex = e;
    return (ne > 0);
  }
  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(edges,ne);
  }
  virtual void apply()
  {
    for (int i = 0; i < ne; ++i){
      if (edgeSwap->run(edges[i])){
        ns++;
        //crv::clearTag(adapter,simplex);
        ma::clearFlag(adapter,edges[i],ma::COLLAPSE | ma::BAD_QUALITY);
        break;
      }
    }
  }
private:
  Adapt* adapter;
  ma::Mesh* mesh;
  ma::Entity* simplex;
  ma::Entity* edges[6];
  ma::EdgeSwap* edgeSwap;
  int md;
  int ne;
public:
  int ns;
};

class EdgeReshaper : public ma::Operator
{
public:
  EdgeReshaper(Adapt* a)
  {
    adapter = a;
    mesh = a->mesh;
    edges[0] = edges[1] = edges[2] = edges[3] = edges[4] = edges[5] = 0;
    simplex = 0;
    md = mesh->getDimension();
    ne = nr = 0;
    qual = makeQuality(mesh,2);
  }
  virtual ~EdgeReshaper()
  {
    delete qual;
  }
  virtual int getTargetDimension() {return md;}
  virtual bool shouldApply(ma::Entity* e)
  {
    int tag = crv::getTag(adapter,e);

    ne = markEdges(mesh,e,tag,edges);
    simplex = e;
    return (ne > 0);
  }
  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(edges,ne);
  }
  virtual void apply()
  {
    for (int i = 0; i < ne; ++i){
      if (!isBoundaryEntity(mesh,edges[i]) &&
          repositionEdge(edges[i])){
        nr++;
        crv::clearTag(adapter,simplex);
        ma::clearFlag(adapter,edges[i],ma::COLLAPSE | ma::BAD_QUALITY);
        break;
      }
    }
  }
private:
  /** \brief reposition second order edge control point based on XJ Luo's
      thesis and bezier.tex in SCOREC/docs repo, only works for second order */
  bool repositionEdge(ma::Entity* edge)
  {
    // lets assume we have an edge we want to fix
    // only support second order for now
    int P = mesh->getShape()->getOrder();
    if (P != 2) return false;

    ma::Entity* verts[4];
    ma::Entity* edges[6];

    ma::Vector pivotPoint;
    ma::Vector edgeVectors[3];
    mesh->getDownward(simplex,0,verts);
    mesh->getDownward(simplex,1,edges);

    // pick a pivotVert, the vertex with the worse jacobian determinant
    ma::Entity* pivotVert;
    int pivotIndex;
    {
      apf::MeshElement* me = apf::createMeshElement(mesh,simplex);

      ma::Entity* edgeVerts[2];
      mesh->getDownward(edge,0,edgeVerts);
      apf::Matrix3x3 J;
      pivotIndex = apf::findIn(verts,4,edgeVerts[0]);
      PCU_ALWAYS_ASSERT(pivotIndex >= 0);

      ma::Vector xi;
      if (md == 3)
      	xi = crv::elem_vert_xi[apf::Mesh::TET][pivotIndex];
      else if (md == 2)
      	xi = crv::elem_vert_xi[apf::Mesh::TRIANGLE][pivotIndex];
      else
      	PCU_ALWAYS_ASSERT(0);

      apf::getJacobian(me,xi,J);

      double j;
      if (md == 3)
      	j = apf::getJacobianDeterminant(J,3);
      else if (md == 2)
      	j = apf::getJacobianDeterminant(J,2);
      else
      	PCU_ALWAYS_ASSERT(0);

      pivotVert = edgeVerts[0];

      int index = apf::findIn(verts,4,edgeVerts[1]);
      PCU_ALWAYS_ASSERT(index >= 0);
      if (md == 3)
      	xi = crv::elem_vert_xi[apf::Mesh::TET][index];
      else if (md == 2)
      	xi = crv::elem_vert_xi[apf::Mesh::TRIANGLE][index];
      else
      	PCU_ALWAYS_ASSERT(0);


      apf::getJacobian(me,xi,J);

      if ((md == 3 && apf::getJacobianDeterminant(J,3) < j) ||
      	  (md == 2 && apf::getJacobianDeterminant(J,2) < j)){
        pivotVert = edgeVerts[1];
        pivotIndex = index;
      }
      apf::destroyMeshElement(me);
    }

    mesh->getPoint(pivotVert,0,pivotPoint);

    // local, of edges around vert, [0,2]
    int edgeIndex = 0;

    for (int i = 0; i < md; ++i){
      // theres only one point, so reuse this...
      if (md == 3)
	edgeVectors[i] = ma::getPosition(mesh,edges[tetVertEdges[pivotIndex][i]])
		      - pivotPoint;
      else if (md == 2)
	edgeVectors[i] = ma::getPosition(mesh,edges[triVertEdges[pivotIndex][i]])
		      - pivotPoint;
      else
      	PCU_ALWAYS_ASSERT(0);

      if ((md == 3 && edges[tetVertEdges[pivotIndex][i]] == edge) ||
	  (md == 2 && edges[triVertEdges[pivotIndex][i]] == edge))
        edgeIndex = i;
    }

    PCU_ALWAYS_ASSERT(edgeIndex >= 0);

    ma::Entity* edge1;
    ma::Entity* edge2;
    // For 3D meshes get the other 2 edges to form a reflection plane
    // Note that the order does not matter
    if (md == 3) {
      edge1 = edges[tetVertEdges[pivotIndex][(edgeIndex+1)%3]];
      edge2 = edges[tetVertEdges[pivotIndex][(edgeIndex+2)%3]];
    }
    // For 2D meshes the invalid edge has to be reflected wrt the
    // other edge connected to the pivot vertex
    // So let:
    // edge1 = edge
    // edge2 = other edge
    else if (md == 2) {
      edge1 = edge;
      edge2 = edges[triVertEdges[pivotIndex][(edgeIndex+1)%2]];
    }


    if (md == 3) {
      // Find the normal to the plane of edge1 and edge2 by first computing the edge
      // tangents at the pivotVert
      ma::Vector t1 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      ma::Vector t2 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      ma::Vector normal = apf::cross(t1, t2);
      double length = normal.getLength();
      double validity = edgeVectors[edgeIndex]*normal;

      if(validity > 1e-10)
	return false;
      ma::Vector oldPoint = ma::getPosition(mesh,edge);
      apf::Adjacent adjacent;
      mesh->getAdjacent(edge,3,adjacent);

      /* mirror the vector edgeVectors[edgeIndex] with respect to the plane
      * perpendicular to the normal. The parameter alpha scales the normal
      * (to the plane) component of the mirrored vector.
      */
      double alpha = 0.5;

      ma::Vector newPoint = pivotPoint + edgeVectors[edgeIndex] -
	normal * (normal * edgeVectors[edgeIndex]) * (1 + alpha) / length / length;

      mesh->setPoint(edge,0,newPoint);

      for (std::size_t i = 0; i < adjacent.getSize(); ++i){
	if (qual->checkValidity(adjacent[i]) < 0){
	  mesh->setPoint(edge,0,oldPoint);
	  return false;
	}
      }
    }
    else if (md == 2) {
      // edge1 which should be the same as edge needs to be modified
      PCU_ALWAYS_ASSERT(edge1 == edge);
      // get the orientation of the simplex using its vertex coordinates
      apf::Vector3 normal = getTriOrientation(mesh, simplex);
      int edge1Idx, edge2Idx;
      edge1Idx = getEdgeIdInSimplex(mesh, simplex, edge1);
      PCU_ALWAYS_ASSERT(edge1Idx > -1);
      edge2Idx = getEdgeIdInSimplex(mesh, simplex, edge2);
      PCU_ALWAYS_ASSERT(edge2Idx > -1);
      ma::Vector t1 = ma::Vector(0.0, 0.0, 0.0);
      ma::Vector t2 = ma::Vector(0.0, 0.0, 0.0);
      switch(pivotIndex) {
      	case 0:
      	  if (edge1Idx == 0 && edge2Idx == 2) {
      	    t1 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      	    t2 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
	  }
	  else if (edge1Idx == 2 && edge2Idx == 0) {
      	    t1 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      	    t2 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
	  }
	  else
	    PCU_ALWAYS_ASSERT(0);
      	  break;
	case 1:
      	  if (edge1Idx == 1 && edge2Idx == 0) {
      	    t1 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      	    t2 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
	  }
	  else if (edge1Idx == 0 && edge2Idx == 1) {
      	    t1 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      	    t2 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
	  }
	  else
	    PCU_ALWAYS_ASSERT(0);
	  break;
	case 2:
      	  if (edge1Idx == 2 && edge2Idx == 1) {
      	    t1 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      	    t2 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
	  }
	  else if (edge1Idx == 1 && edge2Idx == 2) {
      	    t1 = computeEdgeTangentAtVertex(mesh, edge2, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
      	    t2 = computeEdgeTangentAtVertex(mesh, edge1, pivotVert, ma::Matrix(1,0,0,0,1,0,0,0,1));
	  }
	  else
	    PCU_ALWAYS_ASSERT(0);
	  break;
	default:
	  break;
      }

      // At this point if t1xt2 is in the same direction as the normal the triangle is valid
      apf::Vector3 normal2 = apf::cross(t1, t2);
      if (normal * normal2 > 0)
      	return false;

      ma::Vector oldPoint = ma::getPosition(mesh,edge);
      apf::Adjacent adjacent;
      mesh->getAdjacent(edge,md,adjacent);

      /* mirror the vector edgeVectors[edgeIndex] with respect the other edge2 */
      double alpha = 0.5;

      apf::MeshEntity* edge1OppV = 0;
      apf::MeshEntity* edge2OppV = 0;
      apf::MeshEntity* vs[2];

      mesh->getDownward(edge1, 0, vs);
      edge1OppV = (vs[0] == pivotVert) ? vs[1] : vs[0];
      PCU_ALWAYS_ASSERT(edge1OppV);

      mesh->getDownward(edge2, 0, vs);
      edge2OppV = (vs[0] == pivotVert) ? vs[1] : vs[0];
      PCU_ALWAYS_ASSERT(edge2OppV);

      apf::Vector3 L    = ma::getPosition(mesh, edge1OppV) - ma::getPosition(mesh, pivotVert);
      apf::Vector3 refL = ma::getPosition(mesh, edge2OppV) - ma::getPosition(mesh, pivotVert);
      refL = refL/refL.getLength(); // refL is unit

      // get the part of L that is perpendicular to refL
      apf::Vector3 x = L - refL * (refL * L);
      // make sure x is in the plane of edge1 and edge2
      x = x - normal * (normal * x);
      apf::Vector3 deltaX = L - x * alpha * 2;
      ma::Vector newPoint = pivotPoint + deltaX;

      mesh->setPoint(edge,0,newPoint);
      for (std::size_t i = 0; i < adjacent.getSize(); ++i){
	if (qual->checkValidity(adjacent[i]) < 0){
	  mesh->setPoint(edge,0,oldPoint);
	  return false;
	}
      }
    }
    else
      PCU_ALWAYS_ASSERT(0);

    return true;
  }
  Adapt* adapter;
  ma::Mesh* mesh;
  Quality* qual;
  ma::Entity* simplex;
  ma::Entity* edges[6];
  int md;
  int ne;
public:
  int nr;
};

class FaceOptimizer : public ma::Operator
{
public:
  FaceOptimizer(Adapt* a) {
    adapter = a;
    mesh = a->mesh;
    faces[0] = faces[1] = faces[2] = faces[3] = 0;
    simplex = 0;
    md = mesh->getDimension();
    numf = 0;
    ns = 0;
    nf = 0;
  }
  ~FaceOptimizer() {
  }
  virtual int getTargetDimension() {return md;}
  virtual bool shouldApply(ma::Entity* e) {

    std::vector<int> ai = crv::getAllInvalidities(mesh, e);
    int niv = ai.size();
    mesh->getDownward(e, 2, faces);
    numf = 0;
    if (niv != 0) {
      numf = markUniqueFaces(mesh, e, ai, faces);
      simplex = e;
    }
    return (numf > 0);
  }

  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(faces, numf);
  }

  virtual void apply(){
    int invaliditySize = 0;
    //mesh->getDownward(simplex, 2, faces);
    if (numf == 0) return;
    for (int i = 0; i < numf; i++ ) {
      if (mesh->getModelType(mesh->toModel(faces[i])) != 3) continue;
      //CrvFaceOptim *cfo = new CrvFaceOptim(adapter, faces[i], simplex, DETJ);
      CrvFaceOptim *cfo = new CrvFaceOptim(adapter, faces[i], simplex, NIJK);
      cfo->setMaxIter(100);
      cfo->setTol(1e-8);
      if (cfo->run(invaliditySize)) {
	if (invaliditySize < 1) {
	  ns++;
	  delete cfo;
	  break;
	}
      }
      else {
	nf++;
	if (invaliditySize < 1) {
	  delete cfo;
	  break;
	}
      }
      delete cfo;
    }
  }
protected:
  Adapt* adapter;
  ma::Mesh* mesh;
public:
  ma::Entity* faces[4];
  ma::Entity* simplex;
  int numf;
  int md;
  int ns;
  int nf;
};

class EdgeOptimizer : public ma::Operator
{
public:
  EdgeOptimizer(Adapt* a) {
    adapter = a;
    mesh = a->mesh;
    edges[0] = edges[1] = edges[2] = edges[3] = edges[4] = edges[5] = 0;
    simplex = 0;
    md = mesh->getDimension();
    ns = 0;
    nf = 0;
    ne = 0;
  }
  ~EdgeOptimizer() {
  }
  virtual int getTargetDimension() {return md;}
  virtual bool shouldApply(ma::Entity* e) {
    std::vector<int> ai = crv::getAllInvalidities(mesh, e);
    mesh->getDownward(e, 1, edges);
    int niv = ai.size();
    ne = 0;
    if (niv != 0) {
      ne = markAllEdges(mesh, e, ai, edges);
      simplex = e;
    }
    return (ne > 0);
  }

  virtual bool requestLocality(apf::CavityOp* o)
  {
    return o->requestLocality(edges, ne);
  }

  virtual void apply() {
    if (ne == 0) return;
    for (int i = 0; i < ne; i++ ) {
      int invaliditySize = 0;
      CrvEntityOptim* ceo;

      if (mesh->getModelType(mesh->toModel(edges[i])) == 3) {
      	//ceo = new CrvInternalEdgeOptim(adapter, edges[i], simplex, DETJ);
      	ceo = new CrvInternalEdgeOptim(adapter, edges[i], simplex, NIJK);
      }
      else {
      	//ceo = new CrvBoundaryEdgeOptim(adapter, edges[i], simplex, DETJ);
      	ceo = new CrvBoundaryEdgeOptim(adapter, edges[i], simplex, NIJK);
      }

      ceo->setMaxIter(100);
      ceo->setTol(1e-8);

      if (ceo->run(invaliditySize)) {
	ns++;
	if (invaliditySize < 1) {
	  break;
	  delete ceo;
	}
      }
      else {
	nf++;
	if (invaliditySize < 1) {
	  delete ceo;
	  break;
	}
      }
      delete ceo;
    }
  }
protected:
  Adapt* adapter;
  ma::Mesh* mesh;
  ma::Entity* edge;
public:
  int ne;
  int md;
  ma::Entity* edges[6];
  ma::Entity* simplex;
  int ns;
  int nf;
};

static bool isCornerTriAngleLargeMetric(crv::Adapt *a,
    ma::Entity* tri, int index)
{
  ma::Mesh* m = a->mesh;
  // get downward vertexes
  ma::Entity* down[3];
  m->getDownward(tri, 0, down);
  // first get the size field at the center of the tri
  ma::SizeField* sf = a->sizeField;
  apf::Matrix3x3 Q;
  apf::MeshElement* element = apf::createMeshElement(m, tri);
  sf->getTransform(element, apf::Vector3(1./3.,1./3.,1./3.), Q);
  ma::Vector cornerNormal = computeFaceNormalAtVertex(m, tri, down[index], Q);
  apf::destroyMeshElement(element);

  ma::Vector normal = ma::getTriNormal(m,tri);

  // this statement is not exactly a fair comparison, but at least gives
  // some level of control over what is considered an invalid angle
  // an angle is "too large" if the dot product between the corner triangle
  // and the triangle formed by its vertices is negative
  if (cornerNormal*normal < a->input->validQuality)
    return true;

  return false;
}

/* Checks if an angle of a triangle is large (>= 180 degrees)
 * which can be caused by two edges on the boundary curved to it
 *
 */
static ma::Entity* isLargeAngleTriMetric(crv::Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;
  ma::Entity* edges[3];
  m->getDownward(e,1,edges);
  for (int i = 0; i < 3; ++i)
  {
    ma::Entity* e0 = edges[i];
    ma::Entity* e1 = edges[(i+1) % 3];
    if(isBoundaryEntity(m,e0) && isBoundaryEntity(m,e1))
    {
      // TODO: This conditions seems problematic. Think
      // about a better version. My intuition is that the
      // following to edges has to be marked
      // (i+1)%3 and i%3
      if(isCornerTriAngleLargeMetric(a,e,(i+1) % 3)){
        ma::Entity* edge = edges[(i+2) % 3];
        if(!ma::getFlag(a,edge,ma::SPLIT) && !isBoundaryEntity(m,edge)){
          return edge;
        }
      }
    }
  }
  return 0;
}

/* Checks if an angle of a tet is large (>= 180 degrees)
 * which can be caused by two edges on the boundary curved to it
 *
 * An analytic approach, looking at the control point net of points
 * by comparing surface normals of each pair of control net points
 * adjacent to the triangle is an incredibly complex ordering exercise
 * Rather than attempt to do the ordering, sampling the Jacobian at
 * P+1 points is used. A validity check on this edge could also be used
 */

static ma::Entity* isLargeAngleTetMetric(crv::Adapt* a, ma::Entity* e)
{
  ma::Mesh* m = a->mesh;
  // first get the size field at the center of the entity e
  ma::SizeField* sf = a->sizeField;
  apf::Matrix3x3 Q;
  apf::MeshElement* element = apf::createMeshElement(m, e);

  sf->getTransform(element, apf::Vector3(0.25,0.25,0.25), Q);

  apf::destroyMeshElement(element);

  ma::Entity* edges[6];
  ma::Entity* faces[4];
  m->getDownward(e,1,edges);
  m->getDownward(e,2,faces);


  // find edge that matters
  int index = -1;
  for (int i = 0; i < 6; ++i){
    if(isBoundaryEntity(m,faces[edgeFaces[i][0]]) &&
        isBoundaryEntity(m,faces[edgeFaces[i][1]])){
      index = i;
      break;
    }
  }
  if(index < 0) return 0;

  if(!isBoundaryEntity(m,edges[index])) return 0;


  ma::Entity* leftFace  = faces[edgeFaces[index][0]];
  ma::Entity* rightFace = faces[edgeFaces[index][1]];
  double cosAngle = apf::computeCosAngleInTet(m, e, leftFace, rightFace, Q);

  ma::Entity* edge = 0;
  if (cosAngle < -0.9)
    edge = edges[oppEdges[index]];

  return edge;
}


// BAD_QUALITY flag is used on edges to identify them
// as splits for quality, rather than for size refinement
// These two functions handle two seperate situations

// First, triangles are looked at to see if they have an angle > 180 degrees
// There is also a check for the weird situation described above
// Second, tets are looked at to see if they have two faces on the boundary
// where the jacobian determinant is negative along the shared edge,
// indicative of a large angle (curving around a cylinder or sphere).

static int markEdgesOppLargeAnglesTri(Adapt* a)
{
  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  do {
    ma::Iterator* it = m->begin(2);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      if(!hasTwoEntitiesOnBoundary(m,e,1)) continue;
      ma::Entity* edge = isLargeAngleTriMetric(a,e);
      if (edge)
      {
        PCU_ALWAYS_ASSERT(m->getType(edge) == 1);
        ma::setFlag(a,edge,ma::SPLIT);
        ma::setFlag(a,edge,ma::BAD_QUALITY);
        if (a->mesh->isOwned(edge))
          ++count;
      }
    }
    m->end(it);
  } while(count > prev_count);
  return PCU_Add_Long(count);
}

static int markEdgesOppLargeAnglesTet(Adapt* a)
{
  int count = 0;
  int prev_count;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  do {
    ma::Iterator* it = m->begin(3);
    prev_count = count;
    while ((e = m->iterate(it)))
    {
      ma::Entity* edge = isLargeAngleTetMetric(a,e);
      if (edge && !ma::getFlag(a,edge,ma::SPLIT))
      {
        PCU_ALWAYS_ASSERT(m->getType(edge) == 1);
        ma::setFlag(a,edge,ma::SPLIT);
        ma::setFlag(a,edge,ma::BAD_QUALITY);
        if (a->mesh->isOwned(edge))
          ++count;
      }
    }
    m->end(it);
  } while(count > prev_count);
  return PCU_Add_Long(count);
}

/* The whole idea is to do the quality check once,
 * and then use the results to mark edges, etc for
 * fixing
 */

/* static int markEdgesToFix(Adapt* a, int flag) */
/* { */
/*   // do a invalidity check first */
/*   int invalid = markInvalidEntities(a); */
/*   if ( !invalid ) */
/*     return 0; */
/*   int count = 0; */

/*   ma::Mesh* m = a->mesh; */
/*   ma::Entity* e; */

/*   // markEdges could have upto 6 edges marked!!! */
/*   ma::Entity* edges[6]; */
/*   ma::Iterator* it = m->begin(m->getDimension()); */
  
/*   while ((e = m->iterate(it))) */
/*   { */
/*     //int tag = crv::getTag(a,e); */
/*     //int n = markEdges(m,e,tag,edges); */
    
/*     std::vector<int> ai = crv::getAllInvalidities(m, e); */
/*     int niv = ai.size(); */

/*     if (niv != 0) { */
/*       int n = markAllEdges(m, e, ai, edges); */
/*       for (int i = 0; i < n; ++i){ */
/*       	ma::Entity* edge = edges[i]; */
/*       	PCU_ALWAYS_ASSERT(edge); */

/*       	if (edge && !ma::getFlag(a,edge,flag)) { */
/*       	  ma::setFlag(a,edge,flag); */
/*       	  if (a->mesh->isOwned(edge)) */
/*       	    ++count; */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   m->end(it); */

/*   return PCU_Add_Long(count); */
/* } */

/*
static int markFacesToFix(Adapt* a, int flag)
{
  int invalid = markInvalidEntities(a);
  if (!invalid)
    return 0;
  int count = 0;

  ma::Mesh* m = a->mesh;
  ma::Entity* e;

  ma::Entity* faces[4];
  ma::Iterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it)))
  {
    std::vector<int> ai = crv::getAllInvalidities(m, e);
//    int tag = crv::getTag(a, e);
//    int n = markFaces(m, e, tag, faces);
    int n = markUniqueFaces(m, e, ai, faces); 
    for (int i = 0; i < n; i++) {
      ma::Entity* face = faces[i];
      PCU_ALWAYS_ASSERT(face);
      if (face && !ma::getFlag(a, face, flag)) {
      	if (m->getModelType(m->toModel(face)) != 2)
      	  ma::setFlag(a, face, flag);
      	if (a->mesh->isOwned(face))
      	  count++;
      }
    }
  }
  m->end(it);

  return PCU_Add_Long(count);
}
*/
int fixLargeBoundaryAngles(Adapt* a)
{
  double t0 = PCU_Time();
  int count = markEdgesOppLargeAnglesTet(a);
  count += markEdgesOppLargeAnglesTri(a);

  if ( ! count){
    return 0;
  }
  splitEdges(a);
  double t1 = PCU_Time();
  ma::print("split %d boundary edges with "
      "large angles in %f seconds",count,t1-t0);
  return 0;
}

static int collapseInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  ma::Mesh* m = a->mesh;
  int maxDimension = m->getDimension();
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a,1,ma::COLLAPSE));
  int successCount = 0;
  for (int modelDimension=1; modelDimension <= maxDimension; ++modelDimension)
  {
    checkAllEdgeCollapses(a,modelDimension);
    findIndependentSet(a);
    successCount += ma::collapseAllEdges(a, modelDimension);
  }
  successCount = PCU_Add_Long(successCount);
  double t1 = PCU_Time();
  ma::print("Collapsed %d bad edges "
      "in %f seconds",successCount, t1-t0);
  return successCount;
}

static int swapInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  EdgeSwapper es(a);
  ma::applyOperator(a,&es);
  double t1 = PCU_Time();
  ma::print("Swapped %d bad edges "
      "in %f seconds",es.ns, t1-t0);
  return es.ns;
}

static int repositionInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  EdgeReshaper es(a);
  ma::applyOperator(a,&es);
  double t1 = PCU_Time();
  ma::print("Repositioned %d bad edges "
      "in %f seconds",es.nr, t1-t0);
  return es.nr;
}

static int optimizeInvalidEdges(Adapt* a)
{
  double t0 = PCU_Time();
  EdgeOptimizer eo(a);
  ma::applyOperator(a,&eo);
  double t1 = PCU_Time();
  ma::print("Optimized %d bad edges, failed %d edges "
      "in %f seconds",eo.ns, eo.nf, t1-t0);
  return eo.ns;
}

static int optimizeInvalidFaces(Adapt* a)
{
  double t0 = PCU_Time();
  FaceOptimizer fo(a);
  ma::applyOperator(a, &fo);
  double t1 = PCU_Time();
  ma::print("Optimized %d bad faces, failed %d faces "
      "in %f seconds", fo.ns, fo.nf, t1-t0);
  return fo.ns;
}

int fixInvalidFaces(Adapt* a)
{
  return optimizeInvalidFaces(a);
}

int fixInvalidEdges(Adapt* a)
{
  int count = 0;

  if(a->mesh->getShape()->getOrder() == 2)
    count += repositionInvalidEdges(a);
  else {
    if (a->mesh->getDimension() == 3)
      count += optimizeInvalidEdges(a);
  }

  count += collapseInvalidEdges(a);
  count += swapInvalidEdges(a);
  return count;
}

struct IsBadCrvQuality : public ma::Predicate
{
  IsBadCrvQuality(Adapt* a_):a(a_)
  {
    sh = crv::getShapeHandler(a);
  }
  ~IsBadCrvQuality()
  {
    delete sh;
  }
  bool operator()(apf::MeshEntity* e)
  {
    return sh->getQuality(e) < a->input->goodQuality;
  }
  Adapt* a;
  ma::ShapeHandler* sh;
};

int markCrvBadQuality(Adapt* a)
{
  IsBadCrvQuality p(a);
  return ma::markEntities(a, a->mesh->getDimension(), p,
      ma::BAD_QUALITY, ma::OK_QUALITY);
}


int fixLargeAngles(Adapt *a)
{
  if (a->mesh->getDimension() == 3) {
    CrvLargeAngleTetFixer tetFixer(a);
    applyOperator(a, &tetFixer);
    /* return tetFixer.getSuccessCount(); */
  }
  else {
    CrvLargeAngleTriFixer triFixer(a);
    applyOperator(a, &triFixer);
    /* return triFixer.getSuccessCount(); */
  }
  return 0;
}

static int fixShortEdgeElements(Adapt* a)
{
  CrvShortEdgeFixer fixer(a);
  applyOperator(a,&fixer);
  return fixer.nr;
}

void fixCrvElementShapes(Adapt* a)
{
  if ( ! a->input->shouldFixShape)
    return;
  a->input->shouldForceAdaptation = true;
  double t0 = PCU_Time();
  int count = markCrvBadQuality(a);
  int originalCount = count;
  int prev_count;
  int i = 0;
  do {
    if ( ! count)
      break;
    prev_count = count;
    fixLargeAngles(a); // update this
    /* int numOpSuccess = fixLargeAngles(a); // update this */
    /* PCU_Add_Ints(&numOpSuccess,1); */
    /* if (PCU_Comm_Self() == 0) */
    /*   lion_oprint(1,"==> %d large angle fix operations succeeded.\n", numOpSuccess); */
    markCrvBadQuality(a);
    fixShortEdgeElements(a); // update this
    /* int numEdgeRemoved = fixShortEdgeElements(a); // update this */
    /* PCU_Add_Ints(&numEdgeRemoved,1); */
    /* if (PCU_Comm_Self() == 0) */
    /*   lion_oprint(1,"==> %d edges removal operations succeeded.\n", numEdgeRemoved); */
    //fixInvalidFaces(a);
    count = markCrvBadQuality(a);
    ++i;
  } while(count < prev_count && i < 6); // the second conditions is to make sure this does not take long
  double t1 = PCU_Time();
  ma::print("bad shapes down from %d to %d in %f seconds",
        originalCount,count,t1-t0);
  a->input->shouldForceAdaptation = false;
}

}
