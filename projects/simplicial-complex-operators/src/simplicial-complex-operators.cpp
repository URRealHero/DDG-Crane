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

	std::set<size_t> vertices = subset.vertices;
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

	std::set<size_t> edges = subset.edges;
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

	std::set<size_t> faces = subset.faces;
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
    // The star of S is the set of simplices T such that some simplex in S is a face of T.
    // We start with the initial set, as every simplex is a face of itself.
    std::set<size_t> starVertices = subset.vertices;
    std::set<size_t> starEdges = subset.edges;
    std::set<size_t> starFaces = subset.faces;

    // For each vertex in the subset, add all incident edges and faces.
    for (size_t vIdx : subset.vertices) {
        Vertex v = mesh->vertex(vIdx);
        for (Edge e : v.adjacentEdges()) {
            starEdges.insert(e.getIndex());
        }
        for (Face f : v.adjacentFaces()) {
            starFaces.insert(f.getIndex());
        }
    }

    // For each edge in the subset, add all incident faces.
    for (size_t eIdx : subset.edges) {
        Edge e = mesh->edge(eIdx);
        // An edge is a face of at most two faces.
        if (e.halfedge().face() != Face()) {
            starFaces.insert(e.halfedge().face().getIndex());
        }
        if (e.halfedge().twin().face() != Face()) {
            starFaces.insert(e.halfedge().twin().face().getIndex());
        }
    }

    // Faces in the subset don't have any higher-dimensional simplices containing them (in a 2D mesh).

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
	// The link is the closure of star of the subset minus the star of the closure of the subset.
	
	MeshSubset starSubset = star(subset);
	MeshSubset closureSubset = closure(subset);
	MeshSubset closureStarSubset = closure(starSubset);
	MeshSubset starClosureSubset = star(closureSubset);

	// The link is the closure of the subset minus the star of the subset
	MeshSubset linkSubset = closureStarSubset;
	linkSubset.deleteSubset(starClosureSubset);

	return linkSubset;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

	// A simplicial complex is a set of simplices such that every face of a simplex is also in the set.
	//
	
    return subset.equals(closure(subset)); // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

	// Pure simplicial complex: every < k simplicial complex is in at least one k-simplicial complex.
	if (!isComplex(subset))
		return -1;
		
	
	// Check the complexes' dimension
	int d = -1;
	if (!subset.faces.empty()){
		d = 2;
	}
	else if (!subset.edges.empty()) {
		d = 1;
	}
	else if (!subset.vertices.empty()) {
		d = 0;
	}
    else {
        return -1;
    }

    if (d == 2) {
        // Every edge must be part of a selected face.
        for (const size_t eIdx : subset.edges) {
            Edge e = mesh->edge(eIdx);
            bool isAttached = false;
            if (e.halfedge().face() != Face() && subset.faces.count(e.halfedge().face().getIndex())) {
                isAttached = true; // If the edge is part of a selected face
            }
            if (!isAttached && e.halfedge().twin().face() != Face() && subset.faces.count(e.halfedge().twin().face().getIndex())) {
                isAttached = true; // If the twin(reverse) of the edge is part of a selected face
            }
            if (!isAttached) return -1;
        }

        // Every vertex must be part of a selected face.
        for (const size_t vIdx : subset.vertices) {
            Vertex v = mesh->vertex(vIdx);
            bool isAttached = false;
            for (Face f : v.adjacentFaces()) {
                if (subset.faces.count(f.getIndex())) {
                    isAttached = true;
                    break;
                }
            }
            if (!isAttached) return -1;
        }
    } else if (d == 1) {
        // Every vertex must be part of a selected edge.
        for (const size_t vIdx : subset.vertices) {
            Vertex v = mesh->vertex(vIdx);
            bool isAttached = false;
            for (Edge e : v.adjacentEdges()) {
                if (subset.edges.count(e.getIndex())) {
                    isAttached = true;
                    break;
                }
            }
            if (!isAttached) return -1;
        }
    }

    return d;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    std::set<size_t> bdVertices;
    std::set<size_t> bdEdges;

    // Loop through the edges in the subset
    for (size_t eIdx : subset.edges){
        Edge e = mesh->edge(eIdx);
        HalfEdge he = e.halfedge();
        int d = 0;

        Face f1 = he.face();
        if (f1 != Face() && subset.faces.count(f1.getIndex())){
            d++;
        }
        Face f2 = he.twin().face();
        if (f2 != Face() && subset.faces.count(f2.getIndex())){
            d++;
        }
        if (d == 1){
            bdEdges.insert(eIdx); // Edge is in the boundary if it is part of exactly one face
        }
    }

    MeshSubset bdedgeSubset = MeshSubset({}, bdEdges, {}); // Create a MeshSubset for the boundary edges
    return closure(bdedgeSubset); 

}
