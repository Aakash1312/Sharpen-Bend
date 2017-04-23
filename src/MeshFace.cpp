#include "MeshFace.hpp"
#include "MeshVertex.hpp"

void
MeshFace::updateNormal()
{
  // Assume the face is planar.
  VertexConstIterator vi2 = vertices.begin();
  VertexConstIterator vi0 = vi2++;
  VertexConstIterator vi1 = vi2++;

  if (vertices.size() > 3)
  {
    // vi1 might be a concave corner -- we need to add up the cross products at all vertices
    Vector3 sum_cross = Vector3::zero();
    for ( ; vi0 != vertices.end(); ++vi0, ++vi1, ++vi2)
    {
      if (vi1 == vertices.end()) vi1 = vertices.begin();
      if (vi2 == vertices.end()) vi2 = vertices.begin();

      Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
      Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
      sum_cross += e2.cross(e1);
    }

    setNormal(sum_cross.unit());
  }
  else
  {
    Vector3 e1 = (*vi0)->getPosition() - (*vi1)->getPosition();
    Vector3 e2 = (*vi2)->getPosition() - (*vi1)->getPosition();
    setNormal(e2.cross(e1).unit());  // counter-clockwise
  }
}

bool
MeshFace::contains(Vector3 const & p) const
{
  if (vertices.empty()) return false;

  // Generate a ray for the even-odd test, from p to the midpoint of the first halfedge. Ignore degenerate situations for
  // now.
  VertexConstIterator vi    =  verticesBegin();
  VertexConstIterator last  =  vi++;
  Vector3 u = 0.5 * ((*last)->getPosition() + (*vi)->getPosition()) - p;

  long count = 1;  // first halfedge is obviously intersected, since we generated the ray through its midpoint
  for ( ; last != verticesEnd(); ++vi)
  {
    if (vi == verticesEnd()) vi = verticesBegin();

    Vector3 v0 = (*last)->getPosition() - p;
    Vector3 v1 = (*vi)->getPosition()   - p;

    // If winding order is: vector to first vertex, ray, vector to second vertex, then intersects
    Vector3 c0 = v0.cross(u);
    Vector3 c1 = u.cross(v1);
    if (c0.dot(c1) > 0)  // intersects, now check forward or reverse
    {
      // Forward if the vector to the point nearest to p on the line containing the edge makes an acute angle with u.
      //
      // The point p' on line v + t * e closest to point p is v + t0 * e, where t0 = e.dot(p - v) / e.dot(e)
      // (see www.geometrictools.com/Documentation/DistancePointLine.pdf).
      //
      // We translate p to the origin for simpler computations.
      Vector3 edge = v1 - v0;
      Real t0 = -edge.dot(v0) / edge.dot(edge);
      Vector3 u0 = v0 + t0 * edge;

      if (u0.dot(u) > 0)
        count++;
    }

    last = vi;
  }

  return (count % 2 == 1);
}
bool
MeshFace::isFacingOrigin() const
{
  VertexConstIterator vi    =  verticesBegin();
  VertexConstIterator last  =  vi++;
  Vector3 c = Vector3(0, 0, 0);
  for ( ; last != verticesEnd(); ++vi)
  {
    c += ((*last) -> getPosition());
    last = vi;
  }
  if ((c.dot(getNormal())) > 0)
  {
    return false;
  }
return true;
}

double
MeshFace::signedVolume(Vector3 positions[]) const
{
  double sum = 0.0;
  sum -= positions[2].x() * positions[1].y() * positions[0].z();
  sum += positions[1].x() * positions[2].y() * positions[0].z();
  sum += positions[2].x() * positions[0].y() * positions[1].z();
  sum -= positions[0].x() * positions[2].y() * positions[1].z();
  sum -= positions[1].x() * positions[0].y() * positions[2].z();
  sum += positions[0].x() * positions[1].y() * positions[2].z();
  // if (sum < 0)
  // {
  //   sum *= -1.0;
  // }
  return sum/6.0;
}

double
MeshFace::triangleArea(Vector3 positions[]) const
{
  double a = (positions[0] - positions[1]).length();
  double b = (positions[1] - positions[2]).length();
  double c = (positions[0] - positions[2]).length();
  double s = (a + b + c)/2;
  return sqrt(s*(s-a)*(s-b)*(s-c));
}

double
MeshFace::area()
{
  VertexConstIterator vi    =  verticesBegin();
  VertexConstIterator last  =  vi++;
  Polygon3 p = Polygon3();
  for ( ; last != verticesEnd(); ++vi)
  {
    p.addVertex((*last) -> getPosition());
    last = vi;
  }
  std::vector<long> indices;
  double area = 0.0;
  long n = 3 * p.triangulate(indices);
  for (int i = 0; i < n;)
  {
    Vector3 positions[3];
    for (int j = 0; j < 3; ++i, ++j)
    {
      positions[i] = p.getVertex(indices[i]).position;
    }
    area += triangleArea(positions);
  }
  return area;  
}

double
MeshFace::volume() const
{
  VertexConstIterator vi    =  verticesBegin();
  VertexConstIterator last  =  vi++;
  Polygon3 p = Polygon3();
  for ( ; last != verticesEnd(); ++vi)
  {
    p.addVertex((*last) -> getPosition());
    last = vi;
  }
  std::vector<long> indices;
  double volume = 0.0;
  long n = 3 * p.triangulate(indices);
  for (int i = 0; i < n;)
  {
    Vector3 positions[3];
    for (int j = 0; j < 3; ++i, ++j)
    {
      positions[i] = p.getVertex(indices[i]).position;
    }
    volume += signedVolume(positions);
  }
  return volume;
}

void
MeshFace::smoothen(std::vector<MeshVertex*> &v, std::vector<MeshEdge*> &e)
{
  for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
  {
    e.push_back(*ej);
  }
  for (VertexIterator vj = vertices.begin(); vj != vertices.end(); ++vj)
  {
    v.push_back(*vj);
  }
}

void
MeshFace::collectEdges(std::vector<MeshEdge*> &e)
{
  for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
  {
      e.push_back(*ej);
  }
}