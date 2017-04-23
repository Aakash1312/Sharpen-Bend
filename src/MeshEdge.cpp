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

double MeshEdge::getWeight()
{
  FaceIterator fi1 = facesBegin();
  FaceIterator fi2 = ++facesBegin();
  return (((*fi1) -> area()) + ((*fi2) -> area()));
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
    is_brown = true;
    return true;
  }
  return false;
}

int
MeshEdge::setNonSmoothFace(std::vector<MeshFace*> &v)
{
  FaceIterator f1 = facesBegin();
  FaceIterator f2 = ++facesBegin();
  if (((*f1) -> is_smooth) && ((*f2) -> is_smooth))
  {
    return 0;
  }
  if (!((*f1) -> is_smooth) && !((*f2) -> is_smooth))
  {
    return -1;
  }
  if (((*f1) -> is_smooth))
  {
    (*f2) -> is_smooth = true;
    v.push_back((*f2));
  }
  else
  {
    (*f1) -> is_smooth = true;
    v.push_back((*f1));
  }
  return 1;
}