#include "crvEdgeOptim.h"
#include "LBFGS.h"
#include "crv.h"
#include "crvQuality.h"
#include "crvBezier.h"
#include "crvMath.h"
#include <iostream>
#include "apfMatrix.h"

namespace crv{

int CrvEdgeReshapeObjFunc :: getSpaceDim()
{
  return P*d;
}

void CrvEdgeReshapeObjFunc :: getInitEdgeN()
{
  apf::Vector3 intEdgeX;
  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));

  for (int i = 0; i < numENodes; i++) {
    mesh->getPoint(edge, i, intEdgeX);
    ien.push_back(intEdgeX);
  }
}

void CrvEdgeReshapeObjFunc :: getInitFaceN()
{
  apf::Adjacent adjF;
  apf::Vector3 intFaceX;
  mesh->getAdjacent(edge, 2, adjF);
  int numFNodes = mesh->getShape()->countNodesOn(mesh->TRIANGLE);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    for (int j = 0; j < numFNodes; j++) {
      mesh->getPoint(adjF[i], j, intFaceX);
      ifn.push_back(intFaceX);
    }
  }
}

void CrvEdgeReshapeObjFunc :: getInitTetN()
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

std::vector<double> CrvEdgeReshapeObjFunc :: getInitialGuess()
{
  return convertNodeVectorToX(ien);
}

std::vector<double> CrvEdgeReshapeObjFunc :: convertNodeVectorToX(std::vector<apf::Vector3> eftn)
{
  std::vector<double> x0;
  for (int i = 0; i < P-1; i++)
    for (int j = 0; j < 3; j++)
      x0.push_back(eftn[i][j]);

  return x0;
}

std::vector<apf::Vector3> CrvEdgeReshapeObjFunc :: convertXtoNodeVector(const std::vector<double> &x)
{
  std::vector<apf::Vector3> a;
  apf::Vector3 v;
  if (d == 2) {
    int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
    for (int i = 0; i < numENodes; i++) {
      v = {x[d*i], x[d*i + 1], 0.0};
      a.push_back(v);
    }
  }
  if (d == 3) {
    int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
    for (int i = 0; i < numENodes; i++) {
      v = {x[d*i], x[d*i + 1], x[d*i + 2]};
      a.push_back(v);
    }
  }
  return a;
}

void CrvEdgeReshapeObjFunc :: blendTris(const std::vector<apf::Vector3> &egn, std::vector<apf::Vector3> &faceNodes)
{
  apf::Vector3 xi;
  apf::Adjacent adjF;
  apf::Adjacent adjE;
  std::vector<apf::Vector3> cien (ien.begin(),ien.end());

  mesh->getAdjacent(edge, 2, adjF);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
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
	    faceNodes[numFNodes*i+k] = faceNodes[numFNodes*i+k] + (egn[ii] - cien[ii])*factor;
	  }
	}
      }
    }
  }
}

void CrvEdgeReshapeObjFunc :: updateNodes(std::vector<apf::Vector3> ed, std::vector<apf::Vector3> fa, std::vector<apf::Vector3> te)
{
  int numENodes = mesh->getShape()->countNodesOn(mesh->getType(edge));
  for (int i = 0; i < numENodes; i++)
    mesh->setPoint(edge, i, ed[i]);

  apf::Adjacent adjF;
  mesh->getAdjacent(edge, 2, adjF);
  for (std::size_t i = 0; i < adjF.getSize(); i++) {
    int numFNodes = mesh->getShape()->countNodesOn(mesh->getType(adjF[i]));
    for (int j = 0; j < numFNodes; j++)
      mesh->setPoint(adjF[i], j, fa[numFNodes*i + j]);
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
}


void CrvEdgeReshapeObjFunc :: setNodes(std::vector<double> &x)
{
  std::vector<apf::Vector3> en (ien.begin(), ien.end());
  std::vector<apf::Vector3> fn (ifn.begin(), ifn.end());
  std::vector<apf::Vector3> tn (itn.begin(), itn.end());
  en = convertXtoNodeVector(x);
  blendTris(en, fn);

  updateNodes(en, fn, tn);
}

std::vector<double> CrvEdgeReshapeObjFunc :: getVolume()
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
      vol.push_back(v);
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

double CrvEdgeReshapeObjFunc :: computeFValOfElement(apf::NewArray<apf::Vector3> &nodes, double volm)
{
  int weight = 1;
  double sumf = 0;
  if (d == 3) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
	for (int K = 0; K <= d*(P-1); K++) {
	  for (int L = 0; L <= d*(P-1); L++) {
	    if ((I == J && J == K && I == 0) || (J == K && K == L && J == 0) || (I == K && K == L && I == 0) || (I == J && J == L && I == 0))
	      weight = 4;
	    else if ((I == J && I == 0) || (I == K && I == 0) || (I == L && I == 0) || (J == K && J == 0) || (J == L && J == 0) || (K == L && K == 0))
	      weight = 2;
	    else
	      weight = 1;
	    if (I + J + K + L == d*(P-1)) {
	      double f = Nijkl(nodes,P,I,J,K)/(6.0*volm) - 1.0;
	      std::cout<<f<<std::endl;
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

void CrvEdgeReshapeObjFunc :: restoreInitialNodes()
{
  updateNodes(ien, ifn, itn);
}

double CrvEdgeReshapeObjFunc :: getValue(std::vector<double> &x)
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
    restoreInitialNodes();
  }
  return sum;
}

std::vector<double> CrvEdgeReshapeObjFunc :: getGrad(std::vector<double> &x)
{
  double fold = getValue(x);
  double eps = 1.0e-4;
  double h;
  std::vector<double> g;

  for (std::size_t i = 0; i < x.size(); i++) {
    if (std::abs(x[i]) > eps)
      h = eps * std::abs(x[i]);
    else
      h = eps;

    x[i] = x[i] + h;
    double fnew = getValue(x);
    x[i] = x[i] - h;
    g.push_back((fnew-fold)/h);
  }
  return g;
}

void CrvEdgeOptim :: setMaxIter(int n)
{
  iter = n;
}

void CrvEdgeOptim :: setTol(double t)
{
  tol = t;
}

bool CrvEdgeOptim :: run()
{
  CrvEdgeReshapeObjFunc *objF = new CrvEdgeReshapeObjFunc(mesh, edge);
  std::vector<double> x0 = objF->getInitialGuess();
  double f0 = objF->getValue(x0);
  std::cout<< "fval at x0 " << f0<<std::endl;
  LBFGS *l = new LBFGS(tol, iter, x0, objF);

  if (l->run()) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);
    apf::Adjacent adjT;
    mesh->getAdjacent(edge, 3, adjT);
    for (std::size_t i = 0; i < adjT.getSize(); i++) {
      if (checkValidity(mesh, adjT[i], 1) > 1) {
      	objF->restoreInitialNodes();
      	return false;
      }
    }
    return true;
  }
  else
    return false;
}

}
