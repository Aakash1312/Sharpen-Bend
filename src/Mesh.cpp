#include "Mesh.hpp"
#include "MeshVertex.hpp"
#include "MeshEdge.hpp"
#include "MeshFace.hpp"
#include "DGP/FilePath.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <unordered_map>
#include "DGP/Matrix.hpp"
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
Mesh::getChamferEdgeAndFace(std::set<MeshEdge*> &v, std::set<MeshFace*> &f, std::set<MeshFace*> &f1, std::set<MeshFace*> &f2)
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
  filter3bis();
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
        v.insert(&(*ej));
        (*ej).is_chamfer = true;
      }
    }
  }

  std::vector<MeshEdge*> temp_face_edges;
  for (FaceIterator fj = faces.begin(); fj != faces.end(); ++fj)
  {
    temp_face_edges.clear();
    (*fj).collectEdges(temp_face_edges);
    int counter = 0;
    for (int i = 0; i < temp_face_edges.size(); ++i)
    {
      if ((temp_face_edges[i] -> is_chamfer))
      {
        counter++;
      }
    }
    if (counter == temp_face_edges.size())
    {
      (*fj).is_chamfer = true;
      f.insert(&(*fj));
    }
    if (counter == 1)
    {
      f1.insert(&(*fj));
    }
    if (counter == 2)
    {
      f2.insert(&(*fj));
    }
  }
}
void
Mesh::triangulate(MeshFace* f, int num)
{
  if (num == 3)
  {
    MeshVertex* s[3];
    MeshVertex* m[4];
    MeshEdge* ea[3];
    int k = 0;
    for (MeshFace::EdgeIterator ej = f -> edges.begin(); ej != f -> edges.end(); ++ej)
    {
      ea[k] = &(*(*ej));
      k++;
      if (k == 3)
      {
        break;
      }
    }
    s[0] = (ea[0] -> getEndpoint(0));
    s[1] = (ea[0] -> getEndpoint(1));
    if (ea[1] -> hasEndpoint(s[1]))
    {
      s[2] = ea[1] -> getOtherEndpoint(s[1]);
    }
    else
    {
      s[2] = ea[1] -> getOtherEndpoint(s[0]);      
    }

    // m[0] = new MeshVertex((s[0] -> getPosition() + s[1] -> getPosition())/2);
    // m[1] = new MeshVertex((s[1] -> getPosition() + s[2] -> getPosition())/2);
    // m[2] = new MeshVertex((s[0] -> getPosition() + s[2] -> getPosition())/2);
    m[0] = getSharpVertex(ea[0], ea[0] -> orientation);
    m[1] = getSharpVertex(ea[1], ea[1] -> orientation);
    m[2] = getSharpVertex(ea[2], ea[2] -> orientation);
    Vector3 A;
    Vector3 B;
    Vector3 C;
    Vector3 N;
    Vector3 M;
    Vector3 S;
    A = s[0] -> getPosition();
    B = s[1] -> getPosition();
    C = s[2] -> getPosition();
    N = s[0] -> computeSmoothNormal();
    M = s[1] -> computeSmoothNormal();
    S = s[2] -> computeSmoothNormal();
    N.unitize();
    M.unitize();
    S.unitize();
    float a = A.dot(N);
    float b = B.dot(M);
    float c = C.dot(S);
    float result[] = {0, 0, 0};
    float constant[] = {a, b, c};
    float Nm[] = {N.x(), N.y(), N.z()};
    float Mm[] = {M.x(), M.y(), M.z()};
    float Sm[] = {S.x(), S.y(), S.z()};
    // Matrix3 E = Matrix3::fromRows(N, M, S);
    Matrix<float> E(3,3);
    E.setRow(0, Nm);
    E.setRow(1, Mm);
    E.setRow(2, Sm);
    E.invert();
    E.postmulVector(constant, result);
    MeshVertex* new_vertex = new MeshVertex(Vector3(result[0], result[1], result[2]));

    std::vector<MeshVertex*> mesh_vertices;

    mesh_vertices.push_back(s[0]);
    mesh_vertices.push_back(new_vertex);
    mesh_vertices.push_back(m[2]);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(m[2]);
    mesh_vertices.push_back(s[2]);
    mesh_vertices.push_back(new_vertex);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(m[1]);
    mesh_vertices.push_back(s[2]);
    mesh_vertices.push_back(new_vertex);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(m[1]);
    mesh_vertices.push_back(s[1]);
    mesh_vertices.push_back(new_vertex);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(s[1]);
    mesh_vertices.push_back(new_vertex);
    mesh_vertices.push_back(m[0]);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(s[0]);
    mesh_vertices.push_back(new_vertex);
    mesh_vertices.push_back(m[0]);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

  }
  if (num == 2)
  {
    MeshVertex* s[3];
    MeshVertex* m[2];
    MeshEdge* ea[2];
    int k = 0;
    for (MeshFace::EdgeIterator ej = f -> edges.begin(); ej != f -> edges.end(); ++ej)
    {
      if ((*ej) -> is_chamfer)
      {
        ea[k] = &(*(*ej));
        k++;
      }
      if (k == 2)
      {
        break;
      }
    }
    s[0] = (ea[0] -> getEndpoint(0));
    s[1] = (ea[0] -> getEndpoint(1));
    if (ea[1] -> hasEndpoint(s[1]))
    {
      s[2] = ea[1] -> getOtherEndpoint(s[1]);
    }
    else
    {
      s[1] = (ea[0] -> getEndpoint(0));
      s[0] = (ea[0] -> getEndpoint(1));
      s[2] = ea[1] -> getOtherEndpoint(s[1]);
    }

    int orientation = 1;
    m[0] = getSharpVertex(ea[0], ea[0] -> orientation);
    if ((ea[0] -> getEndpoint(0) == ea[1] -> getEndpoint(0)) || (ea[0] -> getEndpoint(1) == ea[1] -> getEndpoint(1)))
    {
      orientation = 1;
    }
    else
    {
      orientation = -1 * (ea[0] -> orientation);
      ea[1] -> orientation = orientation;
    }
    m[1] = getSharpVertex(ea[1], ea[1] -> orientation);

    std::vector<MeshVertex*> mesh_vertices;

    mesh_vertices.push_back(s[0]);
    mesh_vertices.push_back(m[1]);
    mesh_vertices.push_back(m[0]);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(m[0]);
    mesh_vertices.push_back(s[1]);
    mesh_vertices.push_back(m[1]);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(s[0]);
    mesh_vertices.push_back(m[1]);
    mesh_vertices.push_back(s[2]);
    addFace(mesh_vertices.begin(), mesh_vertices.end());
  }

  if (num == 1)
  {
    MeshVertex* s[3];
    MeshVertex* m;
    MeshEdge* se = NULL;
    MeshEdge* nse = NULL;
    for (MeshFace::EdgeIterator ej = f -> edges.begin(); ej != f -> edges.end(); ++ej)
    {
      if (se != NULL && nse != NULL)
      {
        break;
      }
      if ((*ej) -> is_chamfer)
      {
        se = (*ej);
      }
      else
      {
        nse =(*ej);
      }
    }
    s[0] = se -> getEndpoint(0);
    s[1] = se -> getEndpoint(1);
    if (nse->hasEndpoint(s[0])) {
      s[2] = nse->getOtherEndpoint(s[0]);
    } else {
      s[2] = nse->getOtherEndpoint(s[1]);
    }
    m = getSharpVertex(se, se -> orientation);
    std::vector<MeshVertex*> mesh_vertices;
    mesh_vertices.push_back(s[2]);
    mesh_vertices.push_back(s[0]);
    mesh_vertices.push_back(m);
    addFace(mesh_vertices.begin(), mesh_vertices.end());

    mesh_vertices.clear();
    mesh_vertices.push_back(s[2]);
    mesh_vertices.push_back(s[1]);
    mesh_vertices.push_back(m);
    addFace(mesh_vertices.begin(), mesh_vertices.end());
 
  }
}

MeshVertex*
Mesh::getSharpVertex(MeshEdge* ea, int &orientation)
{
  Vector3 A;
  Vector3 B;
  Vector3 N;
  Vector3 M;
  if (orientation == 1)
  {
    A = ea -> getEndpoint(0) -> getPosition();
    B = ea -> getEndpoint(1) -> getPosition();
    N = ea -> getEndpoint(0) -> computeSmoothNormal();
    M = ea -> getEndpoint(1) -> computeSmoothNormal();
  }
  if (orientation == -1)
  {
    A = ea -> getEndpoint(1) -> getPosition();
    B = ea -> getEndpoint(0) -> getPosition();
    N = ea -> getEndpoint(1) -> computeSmoothNormal();
    M = ea -> getEndpoint(0) -> computeSmoothNormal();
  }
    N.unitize();
    M.unitize();
    Vector3 AB = A - B;
    Vector3 H = AB.cross(M.cross(N));
    float h = AB.dot(N);
    float k = (2.0 * (M.dot(N)) * (AB.dot(N))) - (2.0 * AB.dot(M));
    Vector3 new_v = (A+B)/2.0 + (h/k)*H;
    MeshVertex* new_vertex = addVertex(new_v);
    new_vertex -> is_sharp = true;
    return new_vertex;
}
void
Mesh::sharpenMesh()
{
  std::set<MeshEdge*> v;
  std::set<MeshFace*> f;
  std::set<MeshFace*> f1;
  std::set<MeshFace*> f2;
  getChamferEdgeAndFace(v, f, f1, f2);
  for (std::set<MeshFace*>::iterator i = f1.begin(); i != f1.end(); ++i)
  {
    triangulate((*i), 1);
  }
  for (std::set<MeshFace*>::iterator i = f2.begin(); i != f2.end(); ++i)
  {
    triangulate((*i), 2);
  }
  for (std::set<MeshFace*>::iterator i = f.begin(); i != f.end(); ++i)
  {
    triangulate((*i), 3);
  }
  bender();
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


void
Mesh::filter3bis() {
  // Find non-brown edges with 2 adjacent red triangles and tag the endpoints as sharp.
  for (EdgeIterator it = edgesBegin(); it != edgesEnd(); ++it) {
    if (!(it->is_brown)) {
      bool red_ngbrs = true;
      for (MeshEdge::FaceIterator face_it = it->facesBegin(); face_it != it->facesEnd(); ++face_it) {
        if (!((*face_it)->is_smooth)) {
          red_ngbrs = false;
          break;
        }
      }
      if (red_ngbrs) {
        MeshVertex* v1 = it->getEndpoint(0);
        MeshVertex* v2 = it->getEndpoint(1);
        v1->is_sharp = true;
        v2->is_sharp = true;
      }
    }
  }

  // Tag all vertices having a non-manifold red neighbourhood
  // TODO
}

void
Mesh::tagAllSharpEdges() {
  for (EdgeIterator it = edgesBegin(); it != edgesEnd(); ++it) {
    // Sharp edge if both vertices are sharp
    if (it->getEndpoint(0)->is_sharp && it->getEndpoint(1)->is_sharp) {
      it->is_sharp = true;
    }
    // All boundary edges are sharp edges
    if (it->isBoundary()) {
      it->is_sharp = true;
    }
  }
}

void
Mesh::bender() {
  tagAllSharpEdges();
  butterflySubdivide(faces);
}

void
Mesh::butterflySubdivide(std::list<MeshFace> faces) {
  std::map<MeshEdge*, Vector3> new_vertices;
  int count = 0;
  for (std::list<MeshFace> :: iterator it = faces.begin(); it != faces.end(); ++it) {
    for (MeshFace::EdgeIterator it2 = it->edgesBegin(); it2 != it->edgesEnd(); ++it2) {
      if (!((*it2)->bfly_divided)) {
        Vector3 v = divideEdge(*it2);
        new_vertices[*it2] = v;
        (*it2)->bfly_divided = true;
      }
    }
  }
  count = 0;
  for (std::list<MeshFace> :: iterator it = faces.begin(); it != faces.end(); ++it) {
    std::vector<MeshVertex*> vec(3); 
    for (MeshFace::EdgeIterator it2 = it->edgesBegin(); it2 != it->edgesEnd(); ++it2) {

    }    
  }
}

void
Mesh::addFaces(MeshVertex* v[3]) {

}

Vector3
Mesh::divideEdge(MeshEdge* edge) {
  // Degree check
  int k0 = edge->getEndpoint(0)->degree();
  int k1 = edge->getEndpoint(1)->degree();
  Vector3 vertex;
  if (k0 ==  6 && k1 == 6) {
    vertex = divideRegularEdge(edge);
  } else if (k0 != 6 && k1 != 6) {
    Vector3 v1 = divideExtraordinaryEdge(edge, edge->getEndpoint(0));
    Vector3 v2 = divideExtraordinaryEdge(edge, edge->getEndpoint(0));
    vertex = (v1+v2)/2.0;
  } else {
    if (k0 == 6) {
      vertex = divideExtraordinaryEdge(edge, edge->getEndpoint(0));
    } else {
      vertex = divideExtraordinaryEdge(edge, edge->getEndpoint(1));
    }
  }
  return vertex;
}

Vector3
Mesh::divideExtraordinaryEdge(MeshEdge* edge, MeshVertex* vertex) {

}

Vector3
Mesh::divideRegularEdge(MeshEdge* edge) {
  MeshEdge *e1[4];
  MeshEdge *e2[4][2];
  std::vector<MeshFace*> faces;
  Vector3 final_position;
  int sz, i;
  // Assuming manifold mesh
  for (MeshEdge::FaceIterator it = edge->facesBegin(); it != edge->facesEnd(); ++it) {
    faces.push_back(*it);
  }
  sz = faces.size();
  if (sz == 1 || edge->is_sharp) {
    // Boundary Edge or sharp edge
    MeshVertex* v1 = edge->getEndpoint(0);
    MeshVertex* v2 = edge->getEndpoint(2);
    if (v1->isManifoldVertex() && v2->isManifoldVertex()) {
      MeshFace *f1, *f2, *f3, *f4;
      e1[0] = faces[0]->getSuccessor(edge);
      e1[1] = faces[0]->getSuccessor(e1[0]);
      e2[0][0] = findCorrespondingEdge(faces[0], e1[0], edge, f1);
      e2[0][1] = findCorrespondingEdge(faces[0], e1[1], edge, f2);
      e2[1][0] = findCorrespondingEdge(f1, e2[0][0], edge, f3);
      e2[1][1] = findCorrespondingEdge(f2, e2[0][1], edge, f4);
      final_position = computeStencilR2(edge, e2[1][0], e2[1][1]);
    } else if (!(v1->isManifoldVertex() || v2->isManifoldVertex())) {
      final_position = (v1->getPosition() + v2->getPosition())/2.0;
    } else {

    }
  } else if (sz == 2) {
    // 2-faced non-sharp edge
    e1[0] = faces[0]->getSuccessor(edge);
    e1[1] = faces[0]->getSuccessor(e1[0]);
    e1[2] = faces[1]->getSuccessor(edge);
    e1[3] = faces[1]->getSuccessor(e1[0]);
    if (!(e1[0]->isConnectedTo(*e1[3]))) {
      e1[3] = faces[1]->getSuccessor(edge);
      e1[2] = faces[1]->getSuccessor(e1[0]);
    }
    MeshFace* e2_face;
    for (i = 0; i < 4; ++i) {
      // Can possibly be NULL
      e2[i][0] = findCorrespondingEdge(faces[i/2], e1[i], edge, e2_face);
      if (e2_face == NULL) {
        e2[i][1] = NULL;
      } else if (e2_face->getSuccessor(e2[i][0]) != e1[i]) {
        e2[i][1] = e2_face->getSuccessor(e2[i][0]);
      } else {
        e2[i][1] = e2_face->getPredecessor(e2[i][0]);
      }
    }
    bool b1 = isPairOfEdgesSharp(e1[0], e2[3][0]);
    bool b2 = isPairOfEdgesSharp(e2[1][0], e1[2]);
    bool b3 = isPairOfEdgesSharp(e2[0][0], e1[3]);
    bool b4 = isPairOfEdgesSharp(e1[1], e2[2][0]);
    if (b1 && b2) {
      final_position = computeStencilR3(edge); 
    } else if (b3 && b4) {
      final_position = computeStencilR3(edge); 
    } else if (b1 && b4) {
      final_position = computeStencilR5(edge, e2[3][1], e2[2][1]);
    } else if (b2 && b3) {
      final_position = computeStencilR5(edge, e2[0][1], e2[1][1]);
    } else if (b1 || b2 || b3 || b4) {
      if (b1) {
        final_position = computeStencilR4(edge, e2[1][1], e2[3][1], e2[2][1], e1[0], e2[3][0]);
      } else if (b2) {
        final_position = computeStencilR4(edge, e2[3][1], e2[1][1], e2[0][1], e2[1][0], e1[2]);
      } else if (b3) {
        final_position = computeStencilR4(edge, e2[2][1], e2[1][1], e2[0][1], e2[0][0], e1[3]);
      } else if (b4) {
        final_position = computeStencilR4(edge, e2[0][1], e2[3][1], e2[2][1], e1[1], e2[2][0]);
      }
    } else {
      final_position = computeStencilR1(edge, e2[0][1], e2[1][1], e2[3][1], e2[2][1]);
    }
  } else {
    // Non-manifold case
  }
  return final_position;
}

MeshEdge*
Mesh::findCorrespondingEdge(MeshFace* adjacent_face, MeshEdge* adjacent_edge, MeshEdge* connected_edge, MeshFace* face) {
  // Will choose the edge on some face adjacent to adjacent_face through adjacent_edge and connected to connected_edge
  for (MeshEdge::FaceIterator it = adjacent_edge->facesBegin(); it != adjacent_edge->facesEnd(); ++it) {
    if (*it != adjacent_face) {
      face = *it;
      break;
    }
  }
  MeshEdge* edge;
  if (face != NULL && face->getSuccessor(adjacent_edge) != NULL) {
    edge = face->getSuccessor(adjacent_edge);
    if (!(edge->isConnectedTo(*connected_edge))) {
      edge = face->getPredecessor(adjacent_edge);
    }
  }
  return edge;
}

bool
Mesh::isPairOfEdgesSharp(MeshEdge* e1, MeshEdge* e2) {
  if (e1 == NULL || e2 == NULL) {
    return false;
  } else {
    return e1->is_sharp && e2->is_sharp;
  }
}

Vector3
Mesh::computeStencilR3(MeshEdge* edge) {
  Vector3 v1 = edge->getEndpoint(0)->getPosition();
  Vector3 v2 = edge->getEndpoint(1)->getPosition();
  return (v1+v2)/2.0;
}

Vector3
Mesh::computeStencilR5(MeshEdge* edge, MeshEdge* adj1, MeshEdge* adj2) {
  Vector3 v1 = edge->getEndpoint(0)->getPosition();
  Vector3 v2 = edge->getEndpoint(1)->getPosition();
  Vector3 v3, v4, v5;
  if (adj2->hasEndpoint(adj1->getEndpoint(0))) {
    v3 = adj1->getEndpoint(1)->getPosition();
    v4 = adj1->getEndpoint(0)->getPosition();
    v5 = adj2->getOtherEndpoint(adj1->getEndpoint(0))->getPosition();
  } else {
    v3 = adj1->getEndpoint(0)->getPosition();
    v4 = adj1->getEndpoint(1)->getPosition();
    v5 = adj2->getOtherEndpoint(adj1->getEndpoint(1))->getPosition();
  }
  return (v1/2.0) + (v2/2.0) + (v3/-8.0) + (v4/4.0) + (v5/-8.0);
}

Vector3
Mesh::computeStencilR4(MeshEdge* edge, MeshEdge* up, MeshEdge* down1, MeshEdge* down2, MeshEdge* sharp_up, MeshEdge* sharp_down) {
  MeshVertex* mv1 = edge->getEndpoint(0);
  MeshVertex* mv3 = down1->getEndpoint(0);
  MeshVertex* mv4 = down1->getEndpoint(1);
  MeshVertex* mv6 = up->getEndpoint(0);
  Vector3 v[8];
  if (sharp_up->hasEndpoint(mv1)) {
    v[1] = edge->getEndpoint(0)->getPosition();
    v[2] = edge->getEndpoint(1)->getPosition();
  } else {
    v[1] = edge->getEndpoint(1)->getPosition();
    v[2] = edge->getEndpoint(0)->getPosition();
  }
  if (sharp_up->hasEndpoint(mv6)) {
    v[6] = up->getEndpoint(0)->getPosition();
    v[7] = up->getEndpoint(1)->getPosition();
  } else {
    v[6] = up->getEndpoint(1)->getPosition();
    v[7] = up->getEndpoint(0)->getPosition();
  }
  if (sharp_down->hasEndpoint(mv3)) {
    v[3] = down1->getEndpoint(0)->getPosition();
    v[4] = down1->getEndpoint(1)->getPosition();
    v[5] = down2->getOtherEndpoint(mv4)->getPosition();
  } else {
    v[3] = down1->getEndpoint(1)->getPosition();
    v[4] = down1->getEndpoint(0)->getPosition();
    v[5] = down2->getOtherEndpoint(mv3)->getPosition();
  }
  return (3.0*v[1]/8.0) + (5.0*v[2]/8.0) + (-1.0*v[3]/16.0) + (3.0*v[4]/16.0) + (-1.0*v[5]/8.0) + (v[6]/16.0) + (-1.0*v[7]/16.0);
}

Vector3
Mesh::computeStencilR1(MeshEdge* edge, MeshEdge* up1, MeshEdge* up2, MeshEdge* down1, MeshEdge* down2) {
  Vector3 v[9];
  v[1] = edge->getEndpoint(0)->getPosition();
  v[2] = edge->getEndpoint(1)->getPosition();
  if (up2->hasEndpoint(up1->getEndpoint(0))) {
    v[3] = up1->getEndpoint(1)->getPosition();
    v[4] = up1->getEndpoint(0)->getPosition();
    v[5] = up2->getOtherEndpoint(up1->getEndpoint(0))->getPosition();
  } else {
    v[3] = up1->getEndpoint(0)->getPosition();
    v[4] = up1->getEndpoint(1)->getPosition();
    v[5] = up2->getOtherEndpoint(up1->getEndpoint(1))->getPosition();
  }
  if (down2->hasEndpoint(down1->getEndpoint(0))) {
    v[6] = down1->getEndpoint(1)->getPosition();
    v[7] = down1->getEndpoint(0)->getPosition();
    v[8] = down2->getOtherEndpoint(down1->getEndpoint(0))->getPosition();
  } else {
    v[6] = down1->getEndpoint(0)->getPosition();
    v[7] = down1->getEndpoint(1)->getPosition();
    v[8] = down2->getOtherEndpoint(down1->getEndpoint(1))->getPosition();
  }
  return (v[1]/2.0) + (v[2]/2.0) + (-1.0*v[3]/16.0) + (v[4]/8.0) + (-1.0*v[5]/16.0) + (-1.0*v[6]/16.0) + (v[7]/8.0) + (-1.0*v[8]/16.0);
}

Vector3
Mesh::computeStencilR2(MeshEdge* edge, MeshEdge* left, MeshEdge* right) {
  Vector3 v1, v2, v3, v4;
  MeshVertex* mv1, *mv2;
  mv1 = edge->getEndpoint(0);
  mv2 = edge->getEndpoint(1);
  v1 = mv1->getPosition();
  v2 = mv2->getPosition();
  if (left->hasEndpoint(mv1)) {
    v3 = left->getOtherEndpoint(mv1)->getPosition();
  } else {
    v3 = left->getOtherEndpoint(mv2)->getPosition();
  }
  if (right->hasEndpoint(mv1)) {
    v3 = right->getOtherEndpoint(mv1)->getPosition();
  } else {
    v3 = right->getOtherEndpoint(mv2)->getPosition();
  }
  return (9.0*v1/16.0) + (9.0*v2/16.0) + (-1.0*v3/16.0) + (-1.0*v4/16.0);
}
