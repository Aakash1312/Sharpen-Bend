#include "MeshVertex.hpp"
#include "MeshEdge.hpp"
#include "MeshFace.hpp"

MeshEdge *
MeshVertex::getEdgeTo(MeshVertex const * v)
{
  if (v == this) return NULL;

  for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
  {
    Edge * e = *ei;
    if (e->hasEndpoint(v)) return e;
  }

  return NULL;
}

bool
MeshVertex::isBoundary() const
{
  if (edges.empty()) return true;

  for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
    if ((*ei)->isBoundary()) return true;

  return false;
}

void
MeshVertex::addFace(Face * face, bool update_normal)
{
  faces.push_back(face);
  if (update_normal && !has_precomputed_normal)
    addFaceNormal(face->getNormal());
}

void
MeshVertex::removeFace(Face * face)
{
  for (FaceIterator fi = facesBegin(); fi != facesEnd(); )
  {
    if (*fi == face)
    {
      fi = faces.erase(fi);
      if (!has_precomputed_normal) removeFaceNormal(face->getNormal());

      // Keep going, just in case the face somehow got added twice
    }
    else
      ++fi;
  }
}

void
MeshVertex::updateNormal()
{
  if (!faces.empty())
  {
    Vector3 sum_normals = Vector3::zero();
    for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
      sum_normals += (*fi)->getNormal();  // weight by face area?

    normal_normalization_factor = sum_normals.length();
    setNormal(normal_normalization_factor < 1e-20f ? Vector3::zero() : sum_normals / normal_normalization_factor);
  }
  else
  {
    setNormal(Vector3::zero());
    normal_normalization_factor = 0;
  }

  has_precomputed_normal = false;
}

double
MeshVertex::getIncidentFaceAngle(MeshFace* f)
{
  std::vector<MeshEdge*> v;
  Vector3 me[2];
  f -> collectEdges(v);
  int size = v.size();
  for (int i = 0,k =0; i < size, k < 2; ++i)
  {
    if (v[i] -> getEndpoint(0) == this)
    {
      me[k++] = v[i] -> getEndpoint(1) -> getPosition();
    }
    else
    {
      if (v[i] -> getEndpoint(1) == this)
      {
        me[k++] = v[i] -> getEndpoint(0) -> getPosition();
      }      
    }
  }
  Vector3 a = (me[0] - getPosition());
  Vector3 b = (me[1] - getPosition());
  a.unitize();
  b.unitize();
  return Math::fastArcCos(a.dot(b));

}

bool
MeshVertex::isSmooth()
{ 
  for (EdgeIterator ei = edgesBegin(); ei != edgesEnd(); ei++)
  {
    if (!((*ei) -> is_brown))
    {
      return false;
    }
  }
  is_smooth = true;
  return true;
}

void
MeshVertex::markFaces(std::vector<MeshFace*> &v)
{
  for (FaceIterator vj = faces.begin(); vj != faces.end(); ++vj)
  {
    if (!((*vj) -> is_smooth))
    {
      v.push_back(*vj);
      (*vj) -> is_smooth = true;
    }
  }
}

Vector3
MeshVertex::computeSmoothNormal()
{
  Vector3 N(0, 0, 0);
  double sum = 0;
  for (FaceIterator fj = faces.begin(); fj != faces.end(); ++fj)
  {
    if (((*fj) -> is_smooth))
    {
      double weight = getIncidentFaceAngle((*fj));
      N +=  weight * ((*fj) -> getNormal());
      sum += weight;
    }
  }
  return N/sum;
}

bool
MeshVertex::isManifoldVertex() {
  int count = 0;
  for (MeshVertex::EdgeIterator it = edgesBegin(); it != edgesEnd(); ++it) {
    if ((*it)->is_sharp) {
      ++count;
    }
  }
  if (count == 2) {
    return true;
  } else {
    return false;
  }
}