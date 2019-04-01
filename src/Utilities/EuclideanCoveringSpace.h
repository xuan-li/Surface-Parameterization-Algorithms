#ifndef EUCLIDEAN_COVERING_SPACE_H_
#define EUCLIDEAN_COVERING_SPACE_H_

#include <MeshDefinition.h>
#include <list>
#include <queue>
#include "EuclideanGeometry2D.h"
#include <iostream>

#include "MeshFormConverter.h"

struct Segment {
	OpenMesh::Vec2d start_coord;
	OpenMesh::Vec2d end_coord;
	OpenMesh::VertexHandle start;
	OpenMesh::VertexHandle end;
	bool valid;

	double dist() {
		return ((start_coord + end_coord) / 2).norm();
	}

	OpenMesh::Vec2d middle() {
		return (start_coord + end_coord) / 2;
	}

} ;



struct comparator {
	bool operator()(std::list<Segment>::iterator it1, std::list<Segment>::iterator it2) {
		return it1->dist() > it2->dist();
	}
};



class EuclideanCoveringSpaceComputer {
public:
	EuclideanCoveringSpaceComputer(SurfaceMesh &mesh, std::vector<OpenMesh::VertexHandle> cones);
	void Compute();
	void GenerateMeshMatrix(Eigen::MatrixXd &V, Eigen::MatrixXd &NV, Eigen::MatrixXi &F, Eigen::MatrixXd &NF);
protected:
	SurfaceMesh &mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vts_;
	std::list<Segment> boundary_segs_;
	std::vector<std::vector<OpenMesh::Vec2d>> orbits_;
	std::priority_queue<std::list<Segment>::iterator, std::vector<std::list<Segment>::iterator>, comparator> min_heap_;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> next_cone_vtx;

	double max_dist_ = 5;
	
protected:
	void Init();
	void ExtendOrbits(std::function<OpenMesh::Vec2d(OpenMesh::VertexHandle)>);
	void StitchCommonSegment(std::list<Segment>::iterator it1, std::list<Segment>::iterator it2);
	bool Update();
	std::list<Segment>::iterator FindLastValid(std::list<Segment>::iterator it);
	std::list<Segment>::iterator FindNextValid(std::list<Segment>::iterator it);
};

#endif // !EUCLIDEAN_COVERING_SPACE_H_
