#include "MeshEdge.hpp"
#include "MeshFace.hpp"
#include "MeshVertex.hpp"

MeshEdge *
MeshEdge::nextAroundEndpoint(int i)
{
  debugAssertM(i == 0 || i == 1, "MeshEdge: Invalid endpoint index");

  if (numFaces() > 2)  // non-manifold
    return NULL;

  // Find which incident face has this endpoint as the origin of the edge when stepping round the face. The required edge
  // is then the predecessor of this edge around the face.
  for (FaceIterator fi = facesBegin(); fi != facesEnd(); ++fi)
  {
    Face * face = *fi;
    MeshEdge * prev = face->getPredecessor(this);
    if (prev->hasEndpoint(endpoints[i]))  // found it!
      return prev;
  }

  return NULL;
}

double MeshEdge::getFaceAngle()
{

  if (numFaces() != 2)
  {
    return -1;
  }
  FaceIterator fi1 = facesBegin();
  FaceIterator fi2 = ++facesBegin();
  return Math::fastArcCos(((*fi1) -> getNormal()).dot((*fi2) -> getNormal()));
  // return Math::fastArcCos(2.0);
 
}

bool
MeshEdge::isSmooth(double threshold)
{
  double angle = getFaceAngle();
  if (angle < 0)
  {
    return false;
  }
  if (angle < threshold)
  {
    return true;
  }
  return false;
}
