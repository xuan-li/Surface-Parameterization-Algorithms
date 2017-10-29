#include "EuclideanCoveringSpace.h"

EuclideanCoveringSpaceComputer::EuclideanCoveringSpaceComputer(SurfaceMesh & mesh, std::vector<OpenMesh::VertexHandle> cones)
	:mesh_(mesh), cone_vts_(cones)
{
	
}

void EuclideanCoveringSpaceComputer::Compute()
{
	Init();

	while (Update()) {
		while (!min_heap_.top()->valid) {
			auto it = min_heap_.top();
			min_heap_.pop();
			boundary_segs_.erase(it);
		}
		//std::cout << "Min Dist:" << (min_heap_.top())->dist() << std::endl;
	}
}

void EuclideanCoveringSpaceComputer::GenerateMeshMatrix(Eigen::MatrixXd & V, Eigen::MatrixXd & NV, Eigen::MatrixXi & F, Eigen::MatrixXd & NF)
{
	OpenMeshToMatrix(mesh_, V, NV, F, NF);
	
	int n_copies = orbits_[0].size();

	Eigen::MatrixXd new_V(V.rows()* n_copies, 3);
	Eigen::MatrixXi new_F(F.rows() * n_copies, 3);
	Eigen::MatrixXd new_NV(NV.rows() *n_copies, 3);
	Eigen::MatrixXd new_NF(NF.rows() * n_copies, 3);
	new_V.setZero();
	for (int i = 0; i < n_copies; ++i) {
		new_NV.block(i * V.rows(), 0, V.rows(), 3) = NV;
		new_NF.block(i * F.rows(), 0, F.rows(), 3) = NF;
		for(auto viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter)
		{
			OpenMesh::VertexHandle v = *viter;
			for (int k = 0; k < 2; ++k)
				new_V(i * V.rows() + v.idx(), k) = orbits_[v.idx()][i][k];
		}
		new_F.block(i * F.rows(), 0, F.rows(), 3) = F + Eigen::MatrixXi::Constant(F.rows(), 3, i * V.rows());
	}
	V = new_V;
	F = new_F;
	NV = new_NV;
	NF = new_NF;

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

	auto identity = [&mesh_=mesh_](OpenMesh::VertexHandle v)->OpenMesh::Vec2d {return mesh_.texcoord2D(v); };
	ExtendOrbits(identity);
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
	}
}

bool EuclideanCoveringSpaceComputer::Update()
{
	using namespace OpenMesh;

	// pop it and find its neighbors in the segment list;
	auto it = min_heap_.top();
	min_heap_.pop();
	it->valid = false;
	auto prev_it = FindLastValid(it);
	auto next_it = FindNextValid(it);
	

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
	
	auto transformation_vertex = [=, &mesh_=mesh_](VertexHandle v)->Vec2d {
		Vec2d p = mesh_.texcoord2D(v);
		auto result_complex = transformation_complex(transfer_to_complex(p));
		return Vec2d(result_complex.real(), result_complex.imag());
	};

	

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
		seg.start = viter;
		viter = mesh_.property(next_cone_vtx, viter);
		seg.end = viter;

		seg.start_coord = transformation_vertex(seg.start);
		seg.end_coord = transformation_vertex(seg.end);
		seg.valid = true;
		boundary_segs_.insert(it, seg);
		min_heap_.push(std::prev(it));

	}

	boundary_segs_.erase(it);
	StitchCommonSegment(FindLastValid(next_it), next_it);
	StitchCommonSegment(FindNextValid(prev_it), prev_it);
	

	ExtendOrbits(transformation_vertex);

	return true;
}

std::list<Segment>::iterator EuclideanCoveringSpaceComputer::FindLastValid(std::list<Segment>::iterator it)
{
	while (true) {
		if (it == boundary_segs_.begin())
			it = std::prev(boundary_segs_.end());
		else
			it = std::prev(it);
		if (it->valid)
			break;
	}
	return it;
}

std::list<Segment>::iterator EuclideanCoveringSpaceComputer::FindNextValid(std::list<Segment>::iterator it)
{
	while (true) {
		if (std::next(it) == boundary_segs_.end())
			it = boundary_segs_.begin();
		else
			it = std::next(it);
		if (it->valid)
			break;
	}
	return it;
}
