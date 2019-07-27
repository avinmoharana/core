#include "crvModelEdgeOptim.h"
#include "LBFGS.h"
#include "crv.h"
#include "crvQuality.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvTables.h"
#include "crvSnap.h"
#include "crvMath.h"
#include <iostream>
#include "apfMatrix.h"
#include <mth_def.h>

namespace crv{

int CrvModelEdgeReshapeObjFunc :: getSpaceDim()
{
  return P*d;
}

void CrvModelEdgeReshapeObjFunc :: getInitEdgeN()
{
  apf::Vector3 intEdgeX;
  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));

  for (int i = 0; i < numENodes; i++) {
    mesh->getPoint(edge, i, intEdgeX);
    ien.push_back(intEdgeX);
  }
}

apf::Vector3 getInterpolatingPointOnFace(apf::Mesh* mesh, apf::MeshEntity* face, int odr, int whichNode) 
{
  std::vector<apf::Vector3> faceIp;
  apf::Vector3 xi;
  apf::NewArray<apf::Vector3> allCntrlP;
  apf::Element* Fa = apf::createElement(mesh->getCoordinateField(), face);
  apf::getVectorNodes(Fa, allCntrlP);

  int nFn = mesh->getShape()->countNodesOn(mesh->getType(face));
  apf::Vector3 blTri;
  //apf::NewArray<apf::Vector3> rhs(n);
  int j = 0;

  getBezierNodeXi(apf::Mesh::TRIANGLE, odr, whichNode, xi);
  blTri.zero();

  for (int ii = 0; ii < odr+1; ii++) {
    for (int jj = 0; jj < odr+1-ii; jj++) {
      double bFactor = trinomial(odr, ii, jj) * Bijk(ii, jj, odr-ii-jj, 1.-xi[0]-xi[1], xi[0], xi[1]);
      blTri = blTri + allCntrlP[getTriNodeIndex(odr, ii, jj)] * bFactor;
    }
  }

  apf::destroyElement(Fa);
  return blTri;
}

void CrvModelEdgeReshapeObjFunc :: getInitFaceN()
{
  apf::Adjacent adjF;
  apf::Vector3 intFaceX;
  apf::Vector3 ipFN;
  mesh->getAdjacent(edge, 2, adjF);
  int numFNodes = mesh->getShape()->countNodesOn(mesh->TRIANGLE);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    for (int j = 0; j < numFNodes; j++) {
      mesh->getPoint(adjF[i], j, intFaceX);
      ifn.push_back(intFaceX);
      //ipFN = getInterpolatingPointOnFace(mesh, adjF[i], P, j);
      //itpfn.push_back(ipFN);
    }
  }
}

void CrvModelEdgeReshapeObjFunc :: getInitTetN()
{
  apf::Adjacent adjT;
  apf::Vector3 intTetX;
  mesh->getAdjacent(edge, 3, adjT);
  int numTNodes = mesh->getShape()->countNodesOn(mesh->TET);
  for (std::size_t i = 0; i < adjT.getSize(); i++) {
    for (int j = 0; j < numTNodes; j++) {
      mesh->getPoint(adjT[i], j, intTetX);
      itn.push_back(intTetX);
    }
  }
}

std::vector<double> CrvModelEdgeReshapeObjFunc :: getParamCoords()
{
  apf::Vector3 xi;
  apf::Vector3 param;
  std::vector<double> xp;

  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));

  for (int i = 0; i < numENodes; i++) {
    getBezierNodeXi(mesh->getType(edge), P, i, xi);
    transferParametricOnEdgeSplit(mesh, edge, 0.5*(xi[0]+1.0), param);
    for (int j = 0; j < 3; j++) 
      xp.push_back(param[j]);
  }

  apf::Vector3 xif;
  apf::Vector3 paramf;
  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
    if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      for (int j = 0; j < numFNodes; j++) {
      	getBezierNodeXi(mesh->getType(adjF[i]), P, j, xif);
      	transferParametricOnTriSplit(mesh, adjF[i], xif, paramf);
      	for (int k = 0; k < 3; k++)
      	  xp.push_back(paramf[k]);
      }
    }
    else {
      apf::Vector3 intFN;
      for (int j = 0; j < numFNodes; j++) {
      	mesh->getPoint(adjF[i], j, intFN);
      	for (int k = 0; k < 3; k++)
      	  xp.push_back(intFN[k]);
      }
    }
  }
  
  return xp;
}

std::vector<apf::Vector3> CrvModelEdgeReshapeObjFunc :: convertParamCoordsToNodeVector(const std::vector<double> &x)
{
  apf::ModelEntity* me = mesh->toModel(edge);
  std::vector<apf::Vector3> edn;
  std::vector<apf::Vector3> vn = convertXtoNodeVector(x);
  int nENodes = mesh->getShape()->countNodesOn(mesh->EDGE);
  
  for (int i = 0; i < nENodes; i++) {
    apf::Vector3 coorde;
    mesh->snapToModel(me, vn[i], coorde);
    edn.push_back(coorde);
  }

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  int numFNTotal = 0;
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    int nFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));

    if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      apf::ModelEntity* mef = mesh->toModel(adjF[i]);
      for (int j = 0; j < nFNodes; j++) {
	apf::Vector3 coordf;
	mesh->snapToModel(mef, vn[nENodes+numFNTotal], coordf);
        edn.push_back(coordf);
        numFNTotal++;
      }
    }
    else {
      for (int j = 0; j < nFNodes; j++) {
      	edn.push_back(vn[nENodes+numFNTotal]);
      	numFNTotal++;
      }
    }
    
  }

  return edn;
}

std::vector<double> CrvModelEdgeReshapeObjFunc :: getInitialGuess()
{
  return getParamCoords();
}

std::vector<apf::Vector3> CrvModelEdgeReshapeObjFunc :: convertXtoNodeVector(const std::vector<double> &x)
{
  std::vector<apf::Vector3> a;
  apf::Vector3 v;

  if (d == 3) {
    std::size_t num = x.size()/d;
    //int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
    for (std::size_t i = 0; i < num; i++) {
      v = {x[d*i], x[d*i + 1], x[d*i + 2]};
      a.push_back(v);
    }
  }
  return a;
}

void CrvModelEdgeReshapeObjFunc :: blendTris(const std::vector<apf::Vector3> &egn, std::vector<apf::Vector3> &faceNodes)
{
  apf::Vector3 xi;
  apf::Adjacent adjF;
  apf::Adjacent adjE;

  mesh->getAdjacent(edge, 2, adjF);

  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      mesh->getAdjacent(adjF[i], 1, adjE);
      int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));

      for (std::size_t j = 0; j < adjE.getSize(); j++) {
      	if (adjE[j] == edge) {
      	  int jj = 1;

      	  if ( j == 0)
      	    jj = 2;
      	  else if ( j == 1)
      	    jj = 0;
      	  else
      	    jj = 1;

      	  for (int k = 0; k < numFNodes; k++) {
      	    getBezierNodeXi(mesh->TRIANGLE, P, k, xi);
      	    for (std::size_t ii = 0; ii < egn.size(); ii++) {
      	      double factor = 0.5 * (xi[j]/(1-xi[jj])) * binomial(P, ii+1) * intpow(1-xi[jj], ii+1) * intpow(xi[jj], P-ii-1)
      	      	+ 0.5 * (xi[jj]/(1-xi[j])) * binomial(P, ii+1) * intpow(xi[j], ii+1) * intpow(1-xi[j], P-ii-1);
      	      faceNodes[numFNodes*i+k] = faceNodes[numFNodes*i+k] + (egn[ii] - ien[ii])*factor;
	    }
	  }
	}
      }
    }
  }
}

std::vector<apf::Vector3> CrvModelEdgeReshapeObjFunc :: getFaceControlPointsFromInterpolatingPoints(apf::MeshEntity* face, const std::vector<apf::Vector3> &faceInterpolatingP)
{
  std::vector<apf::Vector3> faceControlP;
  apf::Vector3 xi;
  apf::NewArray<apf::Vector3> allCntrlP;
  apf::Element* Fa = apf::createElement(mesh->getCoordinateField(), face);
  apf::getVectorNodes(Fa, allCntrlP);

  int n = mesh->getShape()->countNodesOn(mesh->getType(face));
  mth::Matrix<double> A(n, n);
  mth::Matrix<double> Ainv(n, n);
  apf::NewArray<apf::Vector3> rhs(n);
  int j = 0;

  for (int i = 0; i < n; i++) {
    getBezierNodeXi(apf::Mesh::TRIANGLE, P, i, xi);
    rhs[i].zero();
    for (int ii = 0; ii < P+1; ii++) {
      for (int jj = 0; jj < P+1-ii; jj++) {
      	if (ii == 0 || jj == 0 || (ii+jj == P)) {
      	  double bFactor = trinomial(P, ii, jj) * Bijk(ii, jj, P-ii-jj, 1.-xi[0]-xi[1], xi[0], xi[1]);
          rhs[i] += allCntrlP[getTriNodeIndex(P, ii, jj)] * bFactor;
	}
	else {
	  j = getTriNodeIndex(P, ii, jj) - 3*P; //3P is the total number of nodes on all edges
	  A(i, j) = trinomial(P, ii, jj) * Bijk(ii, jj, P-ii-jj, 1.-xi[0]-xi[1], xi[0], xi[1]);
	}
      }
    }
    rhs[i] = faceInterpolatingP[i] - rhs[i];
  }
  apf::destroyElement(Fa);

  if (n > 1)
    invertMatrixWithPLU(n, A, Ainv);
  else
    Ainv(0,0) = 1./A(0,0);

  for (int i = 0; i < n; i++) {
    apf::Vector3 fcp(0., 0., 0.);
    for (int j = 0; j < n; j++)
      fcp += rhs[j]*Ainv(i, j);   
    faceControlP.push_back(fcp);
  }

  return faceControlP;
}

void CrvModelEdgeReshapeObjFunc :: updateNodes(std::vector<apf::Vector3> ed, std::vector<apf::Vector3> fa, std::vector<apf::Vector3> te, bool isInitialX)
{
  if (!isInitialX) {
    apf::NewArray<apf::Vector3> eIntpCords;
    apf::Element* Ed = apf::createElement(mesh->getCoordinateField(), edge);
    apf::getVectorNodes(Ed, eIntpCords);
    apf::NewArray<double> trsCoff;

    int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
    apf::NewArray<apf::Vector3> contP(nEN);

    for (int i = 0; i < nEN; i++) {
      for (int j = 0; j < 3; j++)
      	eIntpCords[2+i][j] = ed[i][j];
    }

    getBezierTransformationCoefficients(P, mesh->getType(edge), trsCoff);
    crv::convertInterpolationPoints(nEN+2, nEN, eIntpCords, trsCoff, contP);

    for (int i = 0; i < nEN; i++)
      mesh->setPoint(edge, i, contP[i]);

    apf::destroyElement(Ed);
  }
  else {
    int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
    for (int i = 0; i < nEN; i++)
      mesh->setPoint(edge, i, ed[i]);
  }

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  std::vector<apf::Vector3> fNd;

  if (!isInitialX) {
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
      fNd.clear();

      for (int j = 0; j < numFNodes; j++)
      	fNd.push_back(fa[numFNodes*i + j]);

      if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      	std::vector<apf::Vector3> fCN = getFaceControlPointsFromInterpolatingPoints(adjF[i], fNd);

      	for (int j = 0; j < numFNodes; j++)
      	  mesh->setPoint(adjF[i], j, fCN[j]);
      }
      else {
      	for (int j = 0; j < numFNodes; j++)
      	  mesh->setPoint(adjF[i], j, fa[numFNodes*i + j]);
      }
    }
  }
  else {
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
      for (int j = 0; j < numFNodes; j++)
      	mesh->setPoint(adjF[i], j, fa[numFNodes*i+j]);
    }
  }      

  if (d == 3 && P > 3) {
    apf::Adjacent adjT;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      int numTNodes = mesh->getShape()->countNodesOn(mesh->getType(adjT[i]));
      for (int j =0; j < numTNodes; j++)
        mesh->setPoint(adjT[i], j, te[numTNodes*i + j]);
    }
  }

  apf::synchronize(mesh->getCoordinateField()); 
}

void CrvModelEdgeReshapeObjFunc :: setNodes(std::vector<double> &x)
{
  std::vector<apf::Vector3> en;
  std::vector<apf::Vector3> fn;
  std::vector<apf::Vector3> tn;
  std::vector<apf::Vector3> nod = convertParamCoordsToNodeVector(x);
  //std::vector<apf::Vector3> fnewN;

  int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
  for (int i = 0; i < nEN; i++)
    en.push_back(nod[i]);

  //for (size_t i = 0; i < itpfn.size(); i++)
  //  fnewN.push_back(itpfn[i]);

  //blendTris(en, fnewN);

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  int kk = 0;
  int nFN = mesh->getShape()->countNodesOn(mesh->TRIANGLE);

  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    for (int j = 0; j < nFN; j++)
      fn.push_back(nod[nEN+i*nFN+j]);
  }

  if (d > 2 && P > 3) {
   //int nTN = mesh->getShape()->countNodesOn(mesh->TET);
   for (std::size_t i = nEN+nFN; i <nod.size(); i++)
     tn.push_back(nod[i]);
  }
  
  updateNodes(en, fn, tn, false);
}

std::vector<double> CrvModelEdgeReshapeObjFunc :: getVolume()
{
  if (d == 3) {
    apf::Adjacent adjT;
    apf::Matrix3x3 m;
    mesh->getAdjacent(edge, 3, adjT);
    apf::Vector3 point0, point;
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      apf::Adjacent adjV;
      mesh->getAdjacent(adjT[i], 0, adjV);
      for (std::size_t j = 0; j < adjV.getSize(); j++) {
	if ( j == 0)
	  mesh->getPoint(adjV[j], 0, point0);
	else {
	  mesh->getPoint(adjV[j], 0, point);
	  for (int k = 0; k < 3; k++)
	    m[j-1][k] = point[k] - point0[k];
	}
      }
      double v = getDeterminant(m)/6.0;
      vol.push_back(std::abs(v));
    }
  }
  if (d == 2) {
    apf::Adjacent adjF;
    apf::Matrix3x3 m;
    mesh->getAdjacent(edge, 2, adjF);
    apf::Vector3 point;
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      apf::Adjacent adjV;
      mesh->getAdjacent(adjF[i], 0, adjV);
      for (std::size_t j = 0; j < adjV.getSize(); j++) {
	mesh->getPoint(adjV[j], 0, point);
	for (int k = 0; k < 3; k++) {
	  if (k < 2) m[j][k] = point[k];
	  else m[j][2] = 1.0;
	}
      }
      double v = getDeterminant(m);
      vol.push_back(v);
    }
  }
  return vol;
}

static double getAr(apf::Vector3 p0, apf::Vector3 p1, apf::Vector3 p2)
{
  p1 = p1 - p0;
  p2 = p2 - p0;
  double area = (apf::cross(p1, p2)).getLength();
  return (area/2.0);
}

static double getArea(apf::Mesh* mesh, apf::MeshEntity* e)
{
  apf::MeshEntity* ver[3];
  mesh->getDownward(e, 0, ver);
  apf::Vector3 point0, point1, point2;

  mesh->getPoint(ver[0], 0, point0);
  mesh->getPoint(ver[1], 0, point1);
  mesh->getPoint(ver[2], 0, point2);

  point1 = point1 - point0;
  point2 = point2 - point0;

  double area = (apf::cross(point1, point2)).getLength();

  return (area/2.0);
}


double CrvModelEdgeReshapeObjFunc :: computeFValOfElement(apf::NewArray<apf::Vector3> &nodes, double volm)
{
  int weight = 1;
  double sumf = 0;
  if (d == 3) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
	for (int K = 0; K <= d*(P-1); K++) {
	  for (int L = 0; L <= d*(P-1); L++) {
	    if ((I == J && J == K && I == 0) || (J == K && K == L && J == 0) || (I == K && K == L && I == 0) || (I == J && J == L && I == 0))
	      weight = 14;
	    else if ((I == J && I == 0) || (I == K && I == 0) || (I == L && I == 0) || (J == K && J == 0) || (J == L && J == 0) || (K == L && K == 0))
	      weight = 4;
	    else
	      weight = 1;
	    if (I + J + K + L == d*(P-1)) {
	      double f = Nijkl(nodes,P,I,J,K)/(6.0*volm) - 1.0;
	      //std::cout<<"["<<I<<","<<J<<","<<K<<","<<L<<"]   "<<f<<std::endl;
	      sumf = sumf + weight*f*f;
	    }
	  }
	}
      }
    }
  }

  if (d == 2) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
	for (int K = 0; K <= d*(P-1); K++) {
	  if ((I == J && I == 0) || (J == K && J == 0) || (I == K && I == 0))
	    weight = 2;
	  else
	    weight = 1;
	  if (I + J + K == d*(P-1)) {
	    double f = Nijk(nodes,P,I,J)/(4.0*volm) - 1.0;
	    sumf = sumf + weight*f*f;
	  }
	}
      }
    }
  }
  return sumf;
}

void CrvModelEdgeReshapeObjFunc :: restoreInitialNodes()
{
  updateNodes(ien, ifn, itn, true);
}

double CrvModelEdgeReshapeObjFunc :: getValue(std::vector<double> &x)
{
  setNodes(x);

  double sum = 0.0;
  if (d == 2) {
    apf::Adjacent adjF;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 2, adjF);
    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      apf::Element* el = apf::createElement(mesh->getCoordinateField(), adjF[i]);
      apf::getVectorNodes(el, allNodes);
      sum = sum + computeFValOfElement(allNodes, vol[i]);
      apf::destroyElement(el);
    }

    restoreInitialNodes();
  }

  if (d == 3) {
    apf::Adjacent adjT;
    apf::NewArray<apf::Vector3> allNodes;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      apf::Element* el = apf::createElement(mesh->getCoordinateField(), adjT[i]);
      apf::getVectorNodes(el, allNodes);
      sum = sum + computeFValOfElement(allNodes, vol[i]);
      apf::destroyElement(el);
    }

    apf::NewArray<apf::Vector3> eCP;
    apf::Element* Edl = apf::createElement(mesh->getCoordinateField(), edge);
    apf::getVectorNodes(Edl, eCP);
    int nEN = mesh->getShape()->countNodesOn(mesh->getType(edge));
    double xr = 1.0;
    double xir = 1.0;
    //double alpha = 1.0;
    double ad = 0.0;
    apf::Vector3 x1, x0;
    apf::Vector3 xi2, xi1, xi0;
    double beta = 0.0;

    apf::MeshElement* mEdl = apf::createMeshElement(mesh, edge);
    for (int i = 0; i <nEN; i++) {
      apf::Vector3 scord;
      getBezierNodeXi(mesh->getType(edge), P, i, xi1);
      apf::mapLocalToGlobal(mEdl, xi1, scord);
      eCP[2+i] = scord;
    }
    apf::destroyMeshElement(mEdl);

    for (int i = 0; i < nEN; i++) {
      getBezierNodeXi(mesh->getType(edge), P, i, xi1);
      if (i > 0 && i < nEN - 1) {
        getBezierNodeXi(mesh->getType(edge), P, i+1, xi2);
        getBezierNodeXi(mesh->getType(edge), P, i-1, xi0);
      	x1 = eCP[2+i+1] - eCP[2+i];
      	x0 = eCP[2+i] - eCP[2+i-1];
      	xir = (xi2[0] - xi1[0])/(xi1[0] - xi0[0]);
      }
      else if ( i == 0) {
        getBezierNodeXi(mesh->getType(edge), P, i+1, xi2);
      	x1 = eCP[2+i+1] - eCP[2+i];
      	x0 = eCP[2+i] - eCP[0];
      	xir = (xi2[0] - xi1[0])/(xi1[0] + 1); // parent coordinate[-1,1]
      }
      else {
        getBezierNodeXi(mesh->getType(edge), P, i-1, xi0);
      	x1 = eCP[1] - eCP[2+i];
      	x0 = eCP[2+i] - eCP[2+i-1];
      	xir = (1 - xi1[0])/(xi1[0] - xi0[0]);
      }

      //if (0.5*(xi1[0]+1.0) < 1.0 - 0.5*(xi1[0]+1.0))
      //	alpha = 0.5*(xi1[0]+1.0);
      //else
      //	alpha = 1.0 - 0.5*(xi1[0]+1.0);

      xr = (x1.getLength()/x0.getLength());
      ad = (xr/xir - 1);   //(alpha*alpha);
      beta = beta + ad*ad;
      //sum = sum + ad*ad;
    }
    //sum = sum*(1 + beta);
    apf::destroyElement(Edl);
/*
    //if (mesh->getModelType(mesh->toModel(edge)) == 2) {
    //  if (mesh->getModelTag(mesh->toModel(edge)) == 26) {
//	double arandomnumber = 2;
 //     }
 //   }
 
    apf::Adjacent adjF;
    mesh->getAdjacent(edge, 2, adjF);
    apf::NewArray<apf::Vector3> fCP;
    apf::Vector3 xif;
    double a[3] = {1.0, 1.0, 1.0};
    double b[3] = {1.0, 1.0, 1.0};
    double aratio = 0.0;
    //double wfactor = 1.0;
    double gamma = 0.0;

    for (std::size_t i = 0; i < adjF.getSize(); i++) {
      if (mesh->getModelType(mesh->toModel(adjF[i])) == 2) {
      	apf::Element* Fal = apf::createElement(mesh->getCoordinateField(), adjF[i]);
      	apf::getVectorNodes(Fal, fCP);
      	int nFN = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));

      	int vN[3] = {getTriNodeIndex(P, P, 0), getTriNodeIndex(P, 0, P), getTriNodeIndex(P, 0, 0)};

	apf::MeshElement* mfel = apf::createMeshElement(mesh, adjF[i]);
      	for (int j = 0; j < nFN; j++) {
      	  getBezierNodeXi(mesh->getType(adjF[i]), P, j, xif);
	  apf::Vector3 sfcord;
	  apf::mapLocalToGlobal(mfel, xif, sfcord);
	  fCP[3*P + j] = sfcord;
	}
	apf::destroyMeshElement(mfel);

      	double triAphys = getAr(fCP[0], fCP[1], fCP[2]);
      	apf::Vector3 prt0 = {0, 0, 0};
      	apf::Vector3 prt1 = {0, 1, 0};
      	apf::Vector3 prt2 = {0, 0, 1};
      	double triAparnt = getAr(prt0, prt1, prt2);

      	for (int j = 0; j < nFN; j++) {
      	  getBezierNodeXi(mesh->getType(adjF[i]), P, j, xif);
      	  apf::Vector3 xifm = {1.0-xif[0]-xif[1], xif[0], xif[1]};
      	  a[0] = getAr(fCP[3*P + j], fCP[vN[0]], fCP[vN[1]]);
      	  b[0] = getAr(xifm, prt0, prt1);
      	  a[1] = getAr(fCP[3*P + j], fCP[vN[1]], fCP[vN[2]]);
      	  b[1] = getAr(xifm, prt1, prt2);
      	  a[2] = getAr(fCP[3*P + j], fCP[vN[2]], fCP[vN[0]]);
      	  b[2] = getAr(xifm, prt2, prt1);

      	  for (int k = 0; k < 3; k++) {
      	    aratio = (a[k]*triAparnt/(b[k]*triAphys) - 1.0);
      	    gamma = gamma + aratio*aratio;
	  }
	}
	apf::destroyElement(Fal);
      }
    }
*/
    sum = sum*(1 + beta );//+ 0.3*gamma);

    restoreInitialNodes();
  }
  return sum;
}

std::vector<double> CrvModelEdgeReshapeObjFunc :: getGrad(std::vector<double> &x)
{
  //double fold = getValue(x);
  double eps = 1.0e-5;
  double h = eps;
  std::vector<double> g;
  /*
  std::vector<apf::Vector3> par;

  apf::MeshEntity* v[2];
  mesh->getDownward(edge, 0, v);

  for (int i = 0; i<2; i++)
    mesh->getParam(v[i], par[i]);

  double delx = std::abs(par[0][0] - par[1][0]);
  double dely = std::abs(par[0][1] - par[1][1]);
  double delz = std::abs(par[0][2] - par[1][2]);
  double delta = 1.0;
*/
  for (std::size_t i = 0; i < x.size(); i++) {
    /*
    if (i % 3 == 0) delta = delx;
    //if (i % 3 == 1) delta = dely;
    //if (i % 3 == 2) delta = delz;
    
    if (delta < eps) 
      h = eps;
    else 
      h = eps * delta;
*/
    if (std::abs(x[i]) > eps)
      h = eps * std::abs(x[i]);
    else
      h = eps;

    x[i] = x[i] + h;
    double ff = getValue(x);
    x[i] = x[i] - h;
    
    x[i] = x[i] - h;
    double fb = getValue(x);
    x[i] = x[i] + h;
/*
    dff = (ff - fold)/h;
    dbb = (fold - fb)/h;
    df = (dff + dbb)/(2.0);
*/
    double df = (ff - fb)/(2.0 * h);
    g.push_back(df);
  }
  return g;
}

void CrvModelEdgeOptim :: setMaxIter(int n)
{
  iter = n;
}

void CrvModelEdgeOptim :: setTol(double t)
{
  tol = t;
}

bool CrvModelEdgeOptim :: run()
{
  CrvModelEdgeReshapeObjFunc *objF = new CrvModelEdgeReshapeObjFunc(mesh, edge);
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  LBFGS *l = new LBFGS(tol, iter, x0, objF);

  if (l->run()) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);
   
    apf::Adjacent adjT;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      if (checkValidity(mesh, adjT[i], 2) > 1) {
      	objF->restoreInitialNodes();
//	std::cout<<"invalid entity after edop with code "<<checkValidity(mesh, adjT[i], 1) <<std::endl;
      	return false;
      }
    }
    
    //std::cout<<"---------------------------------------------"<<std::endl;
    return true;
  }
  else {
    std::cout<<"*****Optim FAILURE"<<std::endl;
    return false;
  }
}

}
