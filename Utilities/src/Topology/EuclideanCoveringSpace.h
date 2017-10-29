#ifndef EUCLIDEAN_COVERING_SPACE_H_
#define EUCLIDEAN_COVERING_SPACE_H_

#include <MeshDefinition.h>
#include <list>
#include <queue>
#include "..\Geometry\EuclideanGeometry2D.h"

struct Segment {
	OpenMesh::Vec2d start_coord;
	OpenMesh::Vec2d end_coord;
	OpenMesh::VertexHandle start;
	OpenMesh::VertexHandle end;
	bool valid;
	
	/*bool operator()(Segment seg1,  Segment seg2) {
		return ((seg1.start_coord + seg1.end_coord)/2).norm() > ((seg2.start_coord + seg2.end_coord) / 2).norm();
	}*/

	bool operator<(Segment seg) {
		return dist() < seg.dist();
	}

	double dist() {
		return ((start_coord + end_coord) / 2).norm();
	}

	OpenMesh::Vec2d middle() {
		return (start_coord + end_coord) / 2;
	}

} ;

std::ostream & operator<<(std::ostream & os, const Segment & seg)
{
	os << ((seg.start_coord + seg.end_coord) / 2).norm();
	return os;
}

struct comparator {
	bool operator()(std::list<Segment>::iterator it1, std::list<Segment>::iterator it2) {
		it1->dist() > it2->dist();
	}
};



class EuclideanCoveringSpaceComputer {
public:
	EuclideanCoveringSpaceComputer(SurfaceMesh &mesh, std::vector<OpenMesh::VertexHandle> cones);
	void Compute();
protected:
	SurfaceMesh &mesh_;
	std::vector<OpenMesh::VertexHandle> cone_vts_;
	std::list<Segment> boundary_segs_;
	std::vector<std::vector<OpenMesh::Vec2d>> orbits_;
	std::priority_queue<std::list<Segment>::iterator, std::vector<std::list<Segment>::iterator>, comparator> min_heap_;
	OpenMesh::VPropHandleT<OpenMesh::VertexHandle> next_cone_vtx;

	double max_dist_ = 10;
	
protected:
	void Init();
	void ExtendOrbits(std::function<OpenMesh::Vec2d(OpenMesh::VertexHandle)>);
	void StitchCommonSegment(std::list<Segment>::iterator it1, std::list<Segment>::iterator it2);
	bool Update();
};

#endif // !EUCLIDEAN_COVERING_SPACE_H_
