#ifndef ORBIFOLD_MESH_SLICER_H
#define ORBIFOLD_MESH_SLICER_H

#include <MeshDefinition.h>
#include <Topology/MeshSlicer.h>
#include <map>

class OrbifoldMeshSlicer {

public:
	OrbifoldMeshSlicer(SurfaceMesh &mesh);
	SurfaceMesh CutAndSelectSingularities(int n_cones=3);
	std::vector<OpenMesh::VertexHandle> GetConeVertices() { return cone_vertices_; }
	std::vector<std::vector<OpenMesh::VertexHandle>> GetSegments() { return segments_vts_; }

protected:
	SurfaceMesh &mesh_;
	SurfaceMesh sliced_mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vertices_;
	OpenMesh::VertexHandle slice_root_vtx_;
	std::vector<std::vector<OpenMesh::VertexHandle>> segments_vts_;
	
protected:
	void CutBoundaryToSegments();
};
#endif // !ORBIFOLD_MESH_SLICER_H

