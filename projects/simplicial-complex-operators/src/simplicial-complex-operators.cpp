// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <cstddef>

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)

}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    // Each triplet is a (value, row, column) tuple, and SparseMatrix need a NNZ to show how much edges are connected to current
    size_t num_vertices = mesh->nVertices();
    size_t num_edges = mesh->nEdges();
    
    Eigen::SparseMatrix<size_t> A0(num_vertices, num_edges);

    typedef Eigen::Triplet<size_t> T;
    std::vector<T> triplets;
    triplets.reserve(2*num_edges); // reserve space for all edges, each edge has two vertices

    for (const Edge& e : mesh->edges()) {
        size_t edgeIndex = geometry->edgeIndices[e];
        const Vertex &v1 = e.firstVertex();
        const Vertex &v2 = e.secondVertex();

        size_t v1Index = geometry->vertexIndices[v1];
        size_t v2Index = geometry->vertexIndices[v2];

        triplets.emplace_back(v1Index, edgeIndex, 1); // v1 is connected to edge
        triplets.emplace_back(v2Index, edgeIndex, 1); // v2 is connected to edge
    }
    // Create the sparse matrix from the triplets
    A0.setFromTriplets(triplets.begin(), triplets.end());

    return A0;


}


/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    geometry->requireFaceIndices();
    geometry->requireEdgeIndices();

    size_t num_Faces = mesh->nFaces();
    size_t num_Edges = mesh->nEdges();

    Eigen::SparseMatrix<size_t> A1(num_Faces, num_Edges);
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> triplets;
    triplets.reserve(3 * num_Faces); // reserve space for all faces, each face has three edges

    for (const Face& f : mesh->faces()) {
        size_t faceIndex = geometry->faceIndices[f];
        const Halfedge& he = f.halfedge();
        size_t edge1Index = geometry->edgeIndices[he.edge()];
        size_t edge2Index = geometry->edgeIndices[he.next().edge()];
        size_t edge3Index = geometry->edgeIndices[he.next().next().edge()];
        triplets.emplace_back(faceIndex, edge1Index, 1); // face is connected to edge1
        triplets.emplace_back(faceIndex, edge2Index, 1); // face is connected to edge2
        triplets.emplace_back(faceIndex, edge3Index, 1); // face is connected to edge3
    }

    A1.setFromTriplets(triplets.begin(), triplets.end());
    return A1;

}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    set<size_t> vertices = subset.vertices;
    size_t num_vertices = mesh->nVertices();
    
    Vector<size_t> vertexVector = Vector<size_t>::Zero(num_vertices);
    
    for (size_t v : vertices) {
        if (v < num_vertices) {
            vertexVector[v] = 1; // Mark the vertex as present in the subset
        }
    }
    
    return vertexVector;

}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    set<size_t> edges = subset.edges;
    size_t num_edges = mesh->nEdges();

    Vector<size_t> edgeVector = Vector<size_t>::Zero(num_edges);

    for (size_t e : edges) {
        if (e < num_edges) {
            edgeVector[e] = 1; // Mark the edge as present in the subset
        }
    }

    return edgeVector;

}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    set<size_t> faces = subset.faces;
    size_t num_faces = mesh->nFaces();

    Vector<size_t> faceVector = Vector<size_t>::Zero(num_faces);

    for (size_t f : faces) {
        if (f < num_faces) {
            faceVector[f] = 1; // Mark the face as present in the subset
        }
    }

    return faceVector;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    set<size_t> starVertices = subset.vertices;
    set<size_t> starEdges = subset.edges;
    set<size_t> starFaces = subset.faces;

    for (size_t vidx: subset.vertices){
        Vertex v = mesh->vertex(vidx);
        Halfedge start = v.halfedge();
        Halfedge he = v.halfedge();
        if (he == Halfedge()) continue; // skip isolated vertices
        do{
            starEdges.insert(he.edge().getIndex());
            if (he.face() != Face()) {
                starFaces.insert(he.face().getIndex());
            }
            he = he.twin().next(); // reverse direction: find another incident edge
        }
        while (he != start);
    } // end for vertices

    for (size_t eidx: subset.edges){
        Edge e = mesh->edge(eIdx);
        Face f1 = e.halfedge().face(); // get the face of the 
        Face f2 = e.halfedge().twin().face();
        if (f1 != Face()) starFaces.insert(f1.getIndex());
        if (f2 != Face()) starFaces.insert(f2.getIndex());
    } // end for edges


    for (size_t eIdx: starEdges) {
        Edge e = mesh->edge(eIdx);
        Vertex v1 = e.firstVertex();
        Vertex v2 = e.secondVertex();
        starVertices.insert(v1.getIndex());
        starVertices.insert(v2.getIndex());
    } // end for closure vrtices

    for (size_t fIdx: starFaces) {
        Face f = mesh->face(fIdx);
        Halfedge he = f.halfedge();
        Halfedge startHe = he; // store the starting halfedge
        do {
            starEdges.insert(he.edge().getIndex());
            starVertices.insert(he.vertex().getIndex());
            he = he.next();
        } while (he != startHe); // loop through all halfedges of the face
    } // end for closure faces


    return MeshSubset(starVertices, starEdges, starFaces);
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    std::set<size_t> closureVertices = subset.vertices;
    std::set<size_t> closureEdges = subset.edges;
    std::set<size_t> closureFaces = subset.faces;

    // Add all vertices of the edges in the subset
    for (size_t eIdx : subset.edges) {
        Edge e = mesh->edge(eIdx);
        closureVertices.insert(e.firstVertex().getIndex());
        closureVertices.insert(e.secondVertex().getIndex());
    }
    // Add all vertices of the faces in the subset
    for (size_t fIdx : subset.faces) {
        Face f = mesh->face(fIdx);
        Halfedge he = f.halfedge();
        Halfedge startHe = he; // store the starting halfedge
        do {
            closureVertices.insert(he.vertex().getIndex());
            closureEdges.insert(he.edge().getIndex());
            he = he.next();
        } while (he != startHe); // loop through all halfedges of the face
    }

    return MeshSubset(closureVertices, closureEdges, closureFaces);
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    const MeshSubset& closureSubset = closure(subset);
    const MeshSubset& starSubset = star(subset);

    MeshSubset linkSubset = 
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}