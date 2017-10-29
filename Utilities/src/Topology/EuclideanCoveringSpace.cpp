#include "EuclideanCoveringSpace.h"

EuclideanCoveringSpaceComputer::EuclideanCoveringSpaceComputer(SurfaceMesh & mesh, std::vector<OpenMesh::VertexHandle> cones)
	:mesh_(mesh), cone_vts_(cones)
{
	
}

void EuclideanCoveringSpaceComputer::Compute()
{
	Init();
	while (Update()) {
		std::cout << "Min Dist:" << *(min_heap_.top()) << std::endl;
	}
}

void EuclideanCoveringSpaceComputer::Init()
{
	if (!next_cone_vtx.is_valid())
		mesh_.add_property(next_cone_vtx);
	boundary_segs_.clear();
	for (int i = 0; i < cone_vts_.size(); ++i) {
		Segment seg;
		seg.start = cone_vts_[i];
		seg.end = cone_vts_[(i + 1) % cone_vts_.size()];
		mesh_.property(next_cone_vtx, seg.start) = seg.end;
		seg.start_coord = mesh_.texcoord2D(seg.start);
		seg.end_coord = mesh_.texcoord2D(seg.end);
		seg.valid = true;
		boundary_segs_.push_back(seg);
	}
	for (auto it = boundary_segs_.begin(); it != boundary_segs_.end(); ++it) {
		min_heap_.push(it);
	}
	orbits_.clear();
	orbits_.resize(mesh_.n_vertices());
}

void EuclideanCoveringSpaceComputer::ExtendOrbits(std::function<OpenMesh::Vec2d(OpenMesh::VertexHandle)> transformation)
{
	using namespace OpenMesh;
	for (auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		orbits_[v.idx()].push_back(transformation(v));
	}
}

void EuclideanCoveringSpaceComputer::StitchCommonSegment(std::list<Segment>::iterator it1, std::list<Segment>::iterator it2)
{
	if ((it1->middle() - it2->middle()).norm() < 1e-5) {
		it1->valid = false;
		it2->valid = false;
		boundary_segs_.erase(it1);
		boundary_segs_.erase(it2);
	}
}

bool EuclideanCoveringSpaceComputer::Update()
{
	using namespace OpenMesh;
	while (!min_heap_.top()->valid) {
		min_heap_.pop();
	}

	auto it = min_heap_.top();
	min_heap_.pop();
	it->valid = false;
	if (it->dist() > max_dist_) return false;

	VertexHandle start = it->start;
	VertexHandle end = it->end;

	VertexHandle start_equiv;
	VertexHandle end_equiv;

	if (mesh_.data(start).equivalent_vertex().is_valid())
		start_equiv = mesh_.data(start).equivalent_vertex();
	else
		start_equiv = start;

	if (mesh_.data(end).equivalent_vertex().is_valid())
		end_equiv = mesh_.data(end).equivalent_vertex();
	else
		end_equiv = end;

	Vec2d start_equiv_coord = mesh_.texcoord2D(start_equiv);
	Vec2d end_equiv_coord = mesh_.texcoord2D(end_equiv);
	auto transfer_to_complex = [](Vec2d p)->Complex {return Complex(p[0], p[1]); };
	
	auto transformation_complex = ComputeRigidTransformation(
		transfer_to_complex(start_equiv_coord),
		transfer_to_complex(end_equiv_coord),
		transfer_to_complex(it->start_coord),
		transfer_to_complex(it->end_coord)
	);
	
	auto transformation_vertex = [=](VertexHandle v)->Vec2d {
		Vec2d p = mesh_.texcoord2D(v);
		auto result_complex = transformation_complex(transfer_to_complex(p));
		return Vec2d(result_complex.real(), result_complex.imag());
	};

	auto prev_it = std::prev(it);

	VertexHandle viter;
	VertexHandle viter_bound;
	if (mesh_.property(next_cone_vtx, start_equiv) == end_equiv) {
		viter = end_equiv;
		viter_bound = start_equiv;
	}
	else {
		viter = start_equiv;
		viter_bound = end_equiv;
	}
	
	while (viter != end_equiv)
	{
		Segment seg;
		seg.valid = true;
		seg.start = viter;
		viter = mesh_.property(next_cone_vtx, viter);
		seg.end = viter;

		seg.start_coord = transformation_vertex(seg.start);
		seg.end_coord = transformation_vertex(seg.end);
		
		boundary_segs_.insert(it, seg);
		min_heap_.push(std::prev(it));

	}
	it = boundary_segs_.erase(it);
	StitchCommonSegment(prev_it, std::next(prev_it));
	StitchCommonSegment(std::prev(it), it);
	

	ExtendOrbits(transformation_vertex);

	return true;
}
