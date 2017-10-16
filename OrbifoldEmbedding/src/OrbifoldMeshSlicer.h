#ifndef ORBIFOLD_MESH_SLICER_H
#define ORBIFOLD_MESH_SLICER_H

#include <MeshDefinition.h>
#include <MeshSlicer.h>
#include <map>

class OrbifoldMeshSlicer {

public:
	OrbifoldMeshSlicer(SurfaceMesh &mesh);
	SurfaceMesh CutAndSelectSingularities();
	std::vector<OpenMesh::VertexHandle> GetConeVertices() { return cone_vertices_; }

protected:
	SurfaceMesh &mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vertices_;
};
#endif // !ORBIFOLD_MESH_SLICER_H

