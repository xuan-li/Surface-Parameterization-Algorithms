#ifndef HYPERBOLIC_COVERING_SPACE_H_
#define HYPERBOLIC_COVERING_SPACE_H_

#include <MeshDefinition.h>
#include <list>
#include <queue>
#include "..\Geometry\HyperbolicGeometry.h"
#include <iostream>

#include "..\MeshFormConverter.h"

struct HyperbolicSegment {
	OpenMesh::Vec2d start_coord;
	OpenMesh::Vec2d end_coord;
	OpenMesh::VertexHandle start;
	OpenMesh::VertexHandle end;
	bool valid;

	double dist() {
		auto middle = (start_coord + end_coord) / 2;
		return HyperbolicDistance(Complex(0,0), Complex(middle[0], middle[1]));
	}

	OpenMesh::Vec2d middle() {
		return (start_coord + end_coord) / 2;
	}

	double operator - (HyperbolicSegment &s)
	{
		auto middle1 = (start_coord + end_coord) / 2;
		auto middle2 = (s.start_coord + s.end_coord) / 2;
		return HyperbolicDistance(Complex(middle1[0], middle1[1]), Complex(middle2[0], middle2[1]));
	}

};



struct hyperboliccomparator {
	bool operator()(std::list<HyperbolicSegment>::iterator it1, std::list<HyperbolicSegment>::iterator it2) {
		return it1->dist() > it2->dist();
	}
};



class HyperbolicCoveringSpaceComputer {
public:
	HyperbolicCoveringSpaceComputer(SurfaceMesh &mesh, std::vector<OpenMesh::VertexHandle> cones);
	void Compute();
	void GenerateMeshMatrix(Eigen::MatrixXd &V, Eigen::MatrixXd &NV, Eigen::MatrixXi &F, Eigen::MatrixXd &NF);
protected:
	SurfaceMesh &mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vts_;
	std::list<HyperbolicSegment> boundary_segs_;
	std::vector<std::vector<OpenMesh::Vec2d>> orbits_;
	std::priority_queue<std::list<HyperbolicSegment>::iterator, std::vector<std::list<HyperbolicSegment>::iterator>, hyperboliccomparator> min_heap_;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> next_cone_vtx;

	double max_dist_ = 3.5;

protected:
	void Init();
	void ExtendOrbits(std::function<OpenMesh::Vec2d(OpenMesh::VertexHandle)>);
	void StitchCommonSegment(std::list<HyperbolicSegment>::iterator it1, std::list<HyperbolicSegment>::iterator it2);
	bool Update();
	std::list<HyperbolicSegment>::iterator FindLastValid(std::list<HyperbolicSegment>::iterator it);
	std::list<HyperbolicSegment>::iterator FindNextValid(std::list<HyperbolicSegment>::iterator it);
};

#endif // !EUCLIDEAN_COVERING_SPACE_H_
