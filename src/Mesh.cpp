#include "Mesh.hpp"
#include "MeshVertex.hpp"
#include "MeshEdge.hpp"
#include "MeshFace.hpp"
#include "DGP/FilePath.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_map>

MeshEdge *
Mesh::mergeEdges(Edge * e0, Edge * e1)
{
  if (!e0) return e1;
  if (!e1) return e0;

  alwaysAssertM(e0->isCoincidentTo(*e1), std::string(getName()) + ": Edges to merge must have the same endpoints");

  // Transfer faces from e1 to e0
  for (Edge::FaceIterator fi = e1->facesBegin(); fi != e1->facesEnd(); ++fi)
  {
    if (!e0->hasIncidentFace(*fi))
      e0->addFace(*fi);

    (*fi)->replaceEdge(e1, e0);  // in case the face references both e1 and e0, we need to do the replacement even if the above
                                 // 'if' condition is false
  }

  Edge * edges_to_remove[2] = { e1, NULL };
  if (e0->numFaces() <= 0)
    edges_to_remove[1] = e0;

  // Need to cache these before e0 is potentially removed
  Vertex * u = e0->getEndpoint(0);
  Vertex * v = e0->getEndpoint(1);

  for (int i = 0; i < 2; ++i)
  {
    if (!edges_to_remove[i])
      continue;

    edges_to_remove[i]->getEndpoint(0)->removeEdge(edges_to_remove[i]);
    edges_to_remove[i]->getEndpoint(1)->removeEdge(edges_to_remove[i]);

    for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
      if (&(*ej) == edges_to_remove[i])
      {
        edges.erase(ej);
        break;
      }
  }

  Vertex * vertices_to_remove[2] = { NULL, NULL };
  if (u->numFaces() <= 0 && u->numEdges() <= 0) vertices_to_remove[0] = u;
  if (v->numFaces() <= 0 && v->numEdges() <= 0) vertices_to_remove[1] = v;

  for (int i = 0; i < 2; ++i)
  {
    if (!vertices_to_remove[i])
      continue;

    for (VertexIterator vj = vertices.begin(); vj != vertices.end(); ++vj)
      if (&(*vj) == vertices_to_remove[i])
      {
        vertices.erase(vj);
        break;
      }
  }

  return (edges_to_remove[1] == e0 ? NULL : e0);
}

MeshVertex *
Mesh::collapseEdge(Edge * edge)
{
  if (!edge)
    return NULL;

  Vertex * u = edge->getEndpoint(0);
  Vertex * v = edge->getEndpoint(1);

  if (u == v)
  {
    DGP_CONSOLE << getName() << ": Can't collapse edge that is self loop";
    return NULL;
  }

  // Check if u is a repeated vertex in any face. If so, preferentially remove it (one copy at a time)
  bool stop = false;
  for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
  {
    int num_occurrences = 0;
    for (MeshFace::VertexConstIterator vi = fi->verticesBegin(); vi != fi->verticesEnd(); ++vi)
    {
      if (*vi == u)
      {
        num_occurrences++;
        if (num_occurrences >= 2)  // u is the one to remove
        {
          std::swap(u, v);
          stop = true;
          break;
        }
      }
    }

    if (stop)
      break;
  }

  // Remove edge and v from adjacent faces
  std::vector<Face *> edge_faces(edge->facesBegin(), edge->facesEnd());
  for (size_t i = 0; i < edge_faces.size(); ++i)
  {
    Face * face = edge_faces[i];

    // The edge and vertex have to be removed in linked fashion from the loops around the face. This is tricky.
    while (true)
    {
      MeshFace::VertexIterator vi = face->verticesBegin();
      bool found = false;
      for (MeshFace::EdgeIterator ei = face->edgesBegin(); ei != face->edgesEnd(); ++ei, ++vi)
      {
        if (*ei != edge)
          continue;

        face->edges.erase(ei);
        if (*vi == v)
          face->vertices.erase(vi);
        else
        {
          ++vi;
          bool last = false;
          if (vi == face->verticesEnd()) { vi = face->verticesBegin(); last = true; }
          if (*vi == v)
          {
            face->vertices.erase(vi);

            if (last)  // we need to do a shift
            {
              MeshVertex * b = face->vertices.back();
              face->vertices.pop_back();
              face->vertices.push_front(b);
            }
          }
          else
          {
            DGP_ERROR << getName() << ": Vertex not found in expected location";
            return NULL;
          }
        }

        found = true;
        break;
      }

      if (!found)
        break;
    }

    // Wrap up with the easy stuff
    edge->removeFace(face);
    v->removeFace(face);
  }

  assert(edge->numFaces() == 0);

  // Faces adjacent to the edge do not reference it, or v, any more. 'edge' is isolated except for its references to u and v.
  // The mesh is in an inconsistent state since the face does not reference v whereas the next edge around the face does.

  for (Vertex::EdgeIterator vei = v->edgesBegin(); vei != v->edgesEnd(); ++vei)
  {
    (*vei)->replaceVertex(v, u);  // this also makes edge = (u, u)

    if (!u->hasIncidentEdge(*vei))
      u->addEdge(*vei);
  }

  v->edges.clear();

  assert(v->numEdges() == 0);

  // No edges reference v any more. The mesh is in an inconsistent state since faces incident to v but not containing 'edge' do
  // not contain u yet.

  for (Vertex::FaceIterator vfi = v->facesBegin(); vfi != v->facesEnd(); ++vfi)
  {
    (*vfi)->replaceVertex(v, u);  // no-op if we had removed v from this face above

    if (!u->hasIncidentFace(*vfi))
      u->addFace(*vfi);
  }

  v->faces.clear();

  assert(v->numFaces() == 0);

  // No faces reference v any more. The mesh is in a consistent state.

  for (EdgeIterator ei = edges.begin(); ei != edges.end(); ++ei)
    if (&(*ei) == edge)
    {
      edges.erase(ei);
      break;
    }

  u->removeEdge(edge);

  // No more edge. The mesh is in a consistent state

  for (VertexIterator vi = vertices.begin(); vi != vertices.end(); ++vi)
    if (&(*vi) == v)
    {
      vertices.erase(vi);
      break;
    }

  // No more v. The mesh is in a consistent state.

  for (size_t i = 0; i < edge_faces.size(); ++i)
  {
    Face * face = edge_faces[i];
    if (face->numVertices() >= 3)
      continue;

    removeFace(face);
  }

  // All faces shrunk to zero by the edge collapse have been removed.

  bool found = true;
  while (found)
  {
    found = false;

    for (Vertex::EdgeIterator vei = u->edgesBegin(); vei != u->edgesEnd(); ++vei)
    {
      for (Vertex::EdgeIterator vej = u->edgesBegin(); vej != vei; ++vej)
        if ((*vei)->isCoincidentTo(**vej))
        {
          mergeEdges(*vei, *vej);
          found = true;
          break;
        }

      if (found)
        break;
    }
  }

  // All double edges have been collapsed to single edges (this can happen either because faces were shrunk to zero, or because
  // the mesh has genus > 0).

  return u;
}

void
Mesh::draw(Graphics::RenderSystem & render_system, bool draw_edges, bool use_vertex_data, bool send_colors) const
{
  // Three separate passes over the faces is probably faster than using Primitive::POLYGON for each face

  if (draw_edges)
  {
    render_system.pushShapeFlags();
    render_system.setPolygonOffset(true, 1);
  }

  // First try to render as much stuff using triangles as possible
  render_system.beginPrimitive(Graphics::RenderSystem::Primitive::TRIANGLES);
    for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      if (fi->isTriangle()) drawFace(*fi, render_system, use_vertex_data, send_colors);
  render_system.endPrimitive();

  // Now render all quads
  render_system.beginPrimitive(Graphics::RenderSystem::Primitive::QUADS);
    for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
      if (fi->isQuad()) drawFace(*fi, render_system, use_vertex_data, send_colors);
  render_system.endPrimitive();

  // Finish off with all larger polygons
  for (FaceConstIterator fi = facesBegin(); fi != facesEnd(); ++fi)
    if (fi->numEdges() > 4)
    {
      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::POLYGON);
        drawFace(*fi, render_system, use_vertex_data, send_colors);
      render_system.endPrimitive();
    }

  if (draw_edges)
    render_system.popShapeFlags();

  if (draw_edges)
  {
    render_system.pushShader();
    render_system.pushColorFlags();

      render_system.setShader(NULL);
      render_system.setColor(ColorRGBA(0.2, 0.3, 0.7, 1));  // set default edge color

      render_system.beginPrimitive(Graphics::RenderSystem::Primitive::LINES);
        for (EdgeConstIterator ei = edgesBegin(); ei != edgesEnd(); ++ei)
        {
          render_system.sendVertex(ei->getEndpoint(0)->getPosition());
          render_system.sendVertex(ei->getEndpoint(1)->getPosition());
        }
      render_system.endPrimitive();

    render_system.popColorFlags();
    render_system.popShader();
  }
}

bool
Mesh::loadOFF(std::string const & path)
{
  std::ifstream in(path.c_str());
  if (!in)
  {
    DGP_ERROR << "Could not open '" << path << "' for reading";
    return false;
  }
  clear();

  std::string magic;
  if (!(in >> magic) || magic != "OFF")
  {
    DGP_ERROR << "Header string OFF not found at beginning of file '" << path << '\'';
    return false;
  }

  long nv, nf, ne;
  if (!(in >> nv >> nf >> ne))
  {
    DGP_ERROR << "Could not read element counts from OFF file '" << path << '\'';
    return false;
  }

  if (nv < 0 || nf < 0 || ne < 0)
  {
    DGP_ERROR << "Negative element count in OFF file '" << path << '\'';
    return false;
  }

  std::vector<Vertex *> indexed_vertices;
  Vector3 p;
  for (long i = 0; i < nv; ++i)
  {
    if (!(in >> p[0] >> p[1] >> p[2]))
    {
      DGP_ERROR << "Could not read vertex " << indexed_vertices.size() << " from '" << path << '\'';
      return false;
    }

    Vertex * v = addVertex(p);
    if (!v)
      return false;

    indexed_vertices.push_back(v);
  }

  std::vector<Vertex *> face_vertices;
  long num_face_vertices, vertex_index;
  for (long i = 0; i < nf; ++i)
  {
    if (!(in >> num_face_vertices) || num_face_vertices < 0)
    {
      DGP_ERROR << "Could not read valid vertex count of face " << faces.size() << " from '" << path << '\'';
      return false;
    }

    face_vertices.resize(num_face_vertices);
    for (size_t j = 0; j < face_vertices.size(); ++j)
    {
      if (!(in >> vertex_index))
      {
        DGP_ERROR << "Could not read vertex " << j << " of face " << faces.size() << " from '" << path << '\'';
        return false;
      }

      if (vertex_index < 0 || vertex_index >= (long)vertices.size())
      {
        DGP_ERROR << "Out-of-bounds index " << vertex_index << " of vertex " << j << " of face " << faces.size() << " from '"
                  << path << '\'';
        return false;
      }

      face_vertices[j] = indexed_vertices[(size_t)vertex_index];
    }
    addFace(face_vertices.begin(), face_vertices.end());  // ok if this fails, just skip the face with a warning

  }

  setName(FilePath::objectName(path));

  return true;
}

bool
Mesh::loadOBJ(std::string const & path)
{
  std::ifstream in(path.c_str());
  if (!in)
  {
    DGP_ERROR << "Could not open '" << path << "' for reading";
    return false;
  }
  clear();

  // std::string magic;
  // if (!(in >> magic) || magic != "OFF")
  // {
  //   DGP_ERROR << "Header string OFF not found at beginning of file '" << path << '\'';
  //   return false;
  // }

  // long nv, nf, ne;
  // if (!(in >> nv >> nf >> ne))
  // {
  //   DGP_ERROR << "Could not read element counts from OFF file '" << path << '\'';
  //   return false;
  // }

  // if (nv < 0 || nf < 0 || ne < 0)
  // {
  //   DGP_ERROR << "Negative element count in OFF file '" << path << '\'';
  //   return false;
  // }

  std::vector<Vertex *> indexed_vertices;
  Vector3 p;
  std::vector<Vertex *> face_vertices;
  // for (long i = 0; i < nv; ++i)
  // {
  //   if (!(in >> p[0] >> p[1] >> p[2]))
  //   {
  //     DGP_ERROR << "Could not read vertex " << indexed_vertices.size() << " from '" << path << '\'';
  //     return false;
  //   }

  //   Vertex * v = addVertex(p);
  //   if (!v)
  //     return false;

  //   indexed_vertices.push_back(v);
  // }
  std::string input;
  while(!(in.eof()))
  {
    getline(in, input);
    if (input.size() != 0)
    {
      if (input[0] == 'v')
      {
        std::cout <<"YES" << std::endl;
      }
    }
  // std::cout << input << std::endl;
  }
  // while(in >> input) 
  // {
  //     if (input == "v")
  //     {
  //       in >> p[0] >> p[1] >> p[2];
  //       Vertex* v = addVertex(p);
  //       indexed_vertices.push_back(v);
  //     }
  //     if (input == "f")
  //     {
  //       string fv;
  //       while(in >> fv) 
  //       {
          
  //       }
  //     }

  // }

  // std::vector<Vertex *> face_vertices;
  // long num_face_vertices, vertex_index;
  // for (long i = 0; i < nf; ++i)
  // {
  //   if (!(in >> num_face_vertices) || num_face_vertices < 0)
  //   {
  //     DGP_ERROR << "Could not read valid vertex count of face " << faces.size() << " from '" << path << '\'';
  //     return false;
  //   }

  //   face_vertices.resize(num_face_vertices);
  //   for (size_t j = 0; j < face_vertices.size(); ++j)
  //   {
  //     if (!(in >> vertex_index))
  //     {
  //       DGP_ERROR << "Could not read vertex " << j << " of face " << faces.size() << " from '" << path << '\'';
  //       return false;
  //     }

  //     if (vertex_index < 0 || vertex_index >= (long)vertices.size())
  //     {
  //       DGP_ERROR << "Out-of-bounds index " << vertex_index << " of vertex " << j << " of face " << faces.size() << " from '"
  //                 << path << '\'';
  //       return false;
  //     }

  //     face_vertices[j] = indexed_vertices[(size_t)vertex_index];
  //   }
  //   addFace(face_vertices.begin(), face_vertices.end());  // ok if this fails, just skip the face with a warning

  // }

  setName(FilePath::objectName(path));

  return true;
}

bool
Mesh::saveOFF(std::string const & path) const
{
  std::ofstream out(path.c_str(), std::ios::binary);
  if (!out)
  {
    DGP_ERROR << "Could not open '" << path << "' for writing";
    return false;
  }

  out << "OFF\n";
  out << numVertices() << ' ' << numFaces() << " 0\n";

  typedef std::unordered_map<Vertex const *, long> VertexIndexMap;
  VertexIndexMap vertex_indices;
  long index = 0;
  for (VertexConstIterator vi = vertices.begin(); vi != vertices.end(); ++vi, ++index)
  {
    Vector3 const & p = vi->getPosition();
    out << p[0] << ' ' << p[1] << ' ' << p[2] << '\n';

    vertex_indices[&(*vi)] = index;
  }

  for (FaceConstIterator fi = faces.begin(); fi != faces.end(); ++fi)
  {
    out << fi->numVertices();

    for (Face::VertexConstIterator vi = fi->verticesBegin(); vi != fi->verticesEnd(); ++vi)
    {
      VertexIndexMap::const_iterator existing = vertex_indices.find(*vi);
      if (existing == vertex_indices.end())
      {
        DGP_ERROR << "Face references vertex absent from mesh '" << path << '\'';
        return false;
      }

      out << ' ' << existing->second;
    }

    out << '\n';
  }

  return true;
}

bool
Mesh::load(std::string const & path)
{
  std::string path_lc = toLower(path);
  bool status = false;
  // loadOBJ(path);
  if (endsWith(path_lc, ".off"))
    status = loadOFF(path);
  else
  {
    DGP_ERROR << "Unsupported mesh format: " << path;
  }
  return status;
}

bool
Mesh::save(std::string const & path) const
{
  std::string path_lc = toLower(path);
  if (endsWith(path_lc, ".off"))
    return saveOFF(path);

  DGP_ERROR << "Unsupported mesh format: " << path;
  return false;
}

double
Mesh::getThreshold()
{
  double threshold = 0.0;
  double num_edges = 0.0;
  for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
  {
    double angle = (*ej).getFaceAngle();
    // double weight = (*ej).getWeight();
    if (angle >= 0)
    {
      threshold += angle;
      num_edges += 1;
    }
  }
  std::cout << 2*threshold/num_edges << std::endl;
  return 2*threshold/num_edges;
}

void
Mesh::getSmoothEdges(std::vector<MeshEdge*> &v)
{
  double threshold = getThreshold();
  for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
  {
    if ((*ej).isSmooth(threshold))
    {
      v.push_back(&(*ej));
    }
  }

}
void copy(std::vector<MeshEdge*> &A, std::vector<MeshEdge*> &B)
{
  A.clear();
  for (int i = 0; i < B.size(); ++i)
  {
    A.push_back(B[i]);
  }
}
////////////////////////////////////////////////////////////////////////////////////////
void
Mesh::getNonSmoothEdges(std::vector<MeshEdge*> &v)
{
  std::vector<MeshEdge*> smooth_edges;//brown edges
  std::vector<MeshVertex*> smooth_vertices;
  std::vector<MeshFace*> smooth_faces;
  double threshold = getThreshold();
  for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
  {
    if ((*ej).isSmooth(threshold))
    {
      smooth_edges.push_back(&(*ej));
    }
  }

  for (VertexIterator vj = vertices.begin(); vj != vertices.end(); ++vj)
  {
    if ((*vj).isSmooth())
    {
      (*vj).markFaces(smooth_faces);
    }
  }
  bool change = true;
  std::vector<MeshEdge*> temp_edges;
  while(change)
  {
    change = false;
    temp_edges.clear();
    for (int i = 0; i < smooth_edges.size(); ++i)
    {
      int res = smooth_edges[i] -> setNonSmoothFace(smooth_faces);
      if (res <= 0)
      {
        if (res == -1)
        {
          temp_edges.push_back(smooth_edges[i]);
        }
      }
      else
      {
        change = true;
      }
    }
    copy(smooth_edges,temp_edges);
  }
  smooth_edges.clear();
  for (int i = 0; i < smooth_faces.size(); ++i)
  {
    smooth_faces[i] -> smoothen(smooth_vertices, smooth_edges);
  }
  for (int i = 0; i < smooth_vertices.size(); ++i)
  {
    smooth_vertices[i] -> is_smooth = true;
  }
  for (int i = 0; i < smooth_edges.size(); ++i)
  {
    smooth_edges[i] -> is_smooth = true;
  }

  for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
  {
    if (!((*ej).is_smooth))
    {
      if (((*ej).getEndpoint(0)) -> is_smooth && ((*ej).getEndpoint(1)) -> is_smooth)
      {
        v.push_back(&(*ej));
      }
    }
  }

}
////////////////////////////////////////////////////////////////////////////////

void
Mesh::getBrownEdges(std::vector<MeshEdge*> &v)
{
  std::vector<MeshEdge*> smooth_edges;
  std::vector<MeshVertex*> smooth_vertices;
  std::vector<MeshFace*> smooth_faces;
  double threshold = getThreshold();
  for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
  {
    if ((*ej).isSmooth(threshold))
    {
      smooth_edges.push_back(&(*ej));
    }
  }
  copy(v,smooth_edges);
}

void
Mesh::getBrownVertices(std::vector<MeshVertex*> &brown_vertices)
{
    std::vector<MeshVertex*> smooth_vertices;
    double threshold = getThreshold();
    for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
    {
      if ((*ej).isSmooth(threshold))
      {
      }
    }

    for (VertexIterator vj = vertices.begin(); vj != vertices.end(); ++vj)
    {
      if ((*vj).isSmooth())
      {
        brown_vertices.push_back(&(*vj));
      }
    }
}

void
Mesh::getSmoothFaces(std::vector<MeshFace*> &smooth_faces, std::vector<MeshEdge*> smooth_edges)
{
  if (smooth_edges.empty())
  {
    std::vector<MeshVertex*> smooth_vertices;
    double threshold = getThreshold();
    for (EdgeIterator ej = edges.begin(); ej != edges.end(); ++ej)
    {
      if ((*ej).isSmooth(threshold))
      {
        smooth_edges.push_back(&(*ej));
      }
    }

    for (VertexIterator vj = vertices.begin(); vj != vertices.end(); ++vj)
    {
      if ((*vj).isSmooth())
      {
        (*vj).markFaces(smooth_faces);
      }
    }
    // std::vector<MeshEdge*> face_edges;
    // for (FaceIterator fj = faces.begin(); fj != faces.end(); ++fj)
    // {
    //   face_edges.clear();
    //   (*fj).collectEdges(face_edges);
    //   int i;
    //   for (i = 0; i < face_edges.size(); ++i)
    //   {
    //     if (!(face_edges[i] -> is_brown))
    //     {
    //       break;
    //     }
    //   }
    //   if (i == face_edges.size())
    //   {
    //     (*fj).is_smooth = true;
    //     smooth_faces.push_back(&(*fj));
    //   }
    // }
    // return;
    bool change = true;
    std::vector<MeshEdge*> temp_edges;
    change = false;
    temp_edges.clear();
    for (int i = 0; i < smooth_edges.size(); ++i)
    {
      int res = smooth_edges[i] -> setNonSmoothFace(smooth_faces);
      if (res <= 0)
      {
        if (res == -1)
        {
          temp_edges.push_back(smooth_edges[i]);
        }
      }
      else
      {
        change = true;
      }
    }
    copy(smooth_edges,temp_edges);
  }
  else
  {
    bool change = true;
    std::vector<MeshEdge*> temp_edges;
    change = false;
    temp_edges.clear();
    for (int i = 0; i < smooth_edges.size(); ++i)
    {
      int res = smooth_edges[i] -> setNonSmoothFace(smooth_faces);
      if (res <= 0)
      {
        if (res == -1)
        {
          temp_edges.push_back(smooth_edges[i]);
        }
      }
      else
      {
        change = true;
      }
    }
    copy(smooth_edges,temp_edges);
  }
}