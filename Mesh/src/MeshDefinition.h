#ifndef MESHDEFINITION_H
#define MESHDEFINITION_H

#include <OpenMesh/Core/Geometry/VectorT.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>


struct SurfaceMeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	typedef OpenMesh::Vec2d TexCoord2D;

	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::TexCoord2D);
	FaceAttributes(OpenMesh::Attributes::Status);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

	VertexTraits
	{
	public:
		VertexT() :is_singularity(false) {}
		typename Refs::VertexHandle equivalent_vertex() { return equivalent_vertex_; }
		void set_equivalent_vertex(typename Refs::VertexHandle v) { equivalent_vertex = v; }
		bool is_singularity() { return singularity; }
		void set_singularity(bool s) { is_singularity = s; }
	private:
		typename Refs::VertexHandle equivalent_vertex_;
		bool is_singularity;
	};

	HalfedgeTraits
	{
	private:
		double weight_;
	public:
		HalfedgeT() : weight_(0.0) {}
		double weight() { return weight_; }
		void set_weight(double w) { weight_ = w; }
	};

	EdgeTraits
	{
	private:
		double length_;
	public:
		EdgeT():length_(0.0), wedge_(0) {}
		double length() { return length_; }
		void set_length(double l) { length_ = l; }

	};

	FaceTraits
	{
	private:
		
	public:
		FaceT() {};
	};
};

typedef OpenMesh::PolyMesh_ArrayKernelT<SurfaceMeshTraits> BaseSurfaceMesh;

class SurfaceMesh : public BaseSurfaceMesh
{
public:
	SurfaceMesh();
};


#endif