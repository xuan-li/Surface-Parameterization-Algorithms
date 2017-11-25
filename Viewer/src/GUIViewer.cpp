#include "GUIViewer.h"

OTEViewer::OTEViewer()
{
}

void OTEViewer::Init()
{
	InitMenu();
	InitKeyboard();
	InitMouse();
}


void OTEViewer::InitMenu()
{
	callback_init = [this](igl::viewer::Viewer& viewer)
	{
		viewer.ngui->addWindow(Eigen::Vector2i(220, 0), "Control Panel");

		viewer.ngui->addGroup("Mesh IO");
		viewer.ngui->addButton("Load Mesh", [this]() {this->LoadMesh(); });
		viewer.ngui->addButton("Load Texture", [this]() {this->LoadTexture(); });
		viewer.ngui->addButton("Save Mesh", [this]() {this->SaveMesh(); });

		viewer.ngui->addGroup("Cutting System");
		viewer.ngui->addButton("Load Marker", [this] { this->LoadMarker(); });
		viewer.ngui->addButton("Save Marker", [this] {this->SaveMarker(); });
		viewer.ngui->addVariable("Cone Angle", cone_angle_);
		viewer.ngui->addButton("Add Cone", [this]() {this->SetSingularity(cone_angle_); UpdateMeshViewer(); });
		viewer.ngui->addButton("Add Slice", [this]() {this->SetSlice(); UpdateMeshViewer(); });
		viewer.ngui->addButton("Reset Marker", [this]() {this->marker_.ResetMarker(); selected_verts_.clear(); UpdateMeshViewer();  });


		viewer.ngui->addGroup("Core Functions");
		viewer.ngui->addButton("Euclidean Orbifold", [this]() {
			EuclideanOrbifoldSolver solver(this->mesh_, marker_.GetSingularityFlag(), marker_.GetConeAngleFlag(), marker_.GetSliceFlag());
			this->sliced_mesh_ = solver.Compute();
			this->cone_vts_ = solver.ConeVertices();
			hyperbolic_ = false;
			euclidean_ = true;
		});
		viewer.ngui->addButton("Hyperbolic Orbifold", [this]() {
			HyperbolicOrbifoldSolver solver(this->mesh_, marker_.GetSingularityFlag(), marker_.GetSliceFlag());
			this->sliced_mesh_ = solver.Compute();
			this->cone_vts_ = solver.ConeVertices();
			euclidean_ = false;
			hyperbolic_ = true;
		});
		
		viewer.ngui->addButton("Hilbert BFF", [this]() {
			BFFSolver solver(this->mesh_, marker_.GetSingularityFlag(), marker_.GetConeAngleFlag(), marker_.GetSliceFlag());
			this->sliced_mesh_ = solver.Compute(0);
			this->cone_vts_ = solver.ConeVertices();
			euclidean_ = false;
			hyperbolic_ = false;
		});

		viewer.ngui->addButton("Harmonic BFF", [this]() {
			BFFSolver solver(this->mesh_, marker_.GetSingularityFlag(), marker_.GetConeAngleFlag(), marker_.GetSliceFlag());
			this->sliced_mesh_ = solver.Compute(1);
			this->cone_vts_ = solver.ConeVertices();
			euclidean_ = false;
			hyperbolic_ = false;
		});

		viewer.ngui->addButton("Free Boundary BFF", [this]() {
			BFFSolver solver(this->mesh_, marker_.GetSingularityFlag(), marker_.GetConeAngleFlag(), marker_.GetSliceFlag());
			this->sliced_mesh_ = solver.Compute(2);
			this->cone_vts_ = solver.ConeVertices();
			euclidean_ = true;
			hyperbolic_ = false;
		});

		viewer.ngui->addGroup("Viewer Options");


		viewer.ngui->addVariable<ShowOption>(
			"Show Option",
			[this](const ShowOption & v) {  this->show_option_ = v; this->UpdateMeshViewer(); },
			[this]()->ShowOption { return this->show_option_; }
		)->setItems({ "Original", "Sliced", "Embedding", "Covering Space" });
		
		viewer.ngui->addVariable<bool>(
			"Show Boundaries",
			[this](const bool &v) {this->show_boundaries_ = v; this->UpdateMeshViewer(); },
			[this]() -> bool {return this->show_boundaries_; }
		);

		viewer.ngui->addVariable<bool>(
			"Show Slices and Cones",
			[this](const bool &v) {this->show_slice_ = v; this->UpdateMeshViewer(); },
			[this]() -> bool {return this->show_slice_; }
		);

		viewer.ngui->addVariable<bool>(
			"Show Vertex Labels",
			[this](const bool &v) {this->show_vertex_labels_ = v; this->UpdateMeshViewer(); },
			[this]() -> bool {return this->show_vertex_labels_; }
		);


		viewer.screen->performLayout();
		return false; 
	};

	
}

void OTEViewer::InitKeyboard()
{
	callback_key_down = [this](igl::viewer::Viewer& viewer, unsigned char key, int modifier) 
	{
		if (key == 'S') {
			selection_mode_ = true;
		}
		return false;
	};
	callback_key_up = [this](igl::viewer::Viewer& viewer, unsigned char key, int modifier)
	{
		if (key == 'S') {
			selection_mode_ = false;
		}
		return false;
	};
}

void OTEViewer::InitMouse()
{
	callback_mouse_down = [this](igl::viewer::Viewer& viewer, int button, int modifier)->bool
	{
		if ( selection_mode_ && show_option_ == ORIGINAL && button == GLFW_MOUSE_BUTTON_1) {
			Eigen::Vector4f viewport = viewer.core.viewport;
			Eigen::Vector2f center(viewport(2) / 2., viewport(3)/2.);
			double x = viewer.down_mouse_x;
			double y = viewer.down_mouse_y;
			x = (viewer.down_mouse_x - center(0)) / (viewport(2) / 2.);
			y = -(viewer.down_mouse_y - center(1)) / (viewport(3) / 2.);
			this->FindIntersection(x, y);
			this->UpdateMeshViewer();
		} 

		return false;
	};
}

void OTEViewer::LoadMesh()
{
	std::string fname = igl::file_dialog_open();
	if (fname.length() == 0)
		return;
	OpenMesh::IO::Options opt;
	opt += OpenMesh::IO::Options::VertexTexCoord;
	OpenMesh::IO::read_mesh(mesh_, fname, opt);
	NormalizeMesh(mesh_);
	UpdateMeshData(mesh_);
	UpdateTextureCoordData(mesh_);
	marker_.SetObject(mesh_);
	marker_.ResetMarker();
}

void OTEViewer::LoadTexture()
{
	std::string fname = igl::file_dialog_open();
	if (fname.length() == 0)
		return;
	igl::png::texture_from_png(fname, R_, G_, B_, A_);
	data.set_texture(R_, G_, B_, A_);
}

void OTEViewer::SaveMesh()
{
	std::string fname = igl::file_dialog_save();
	if (fname.length() == 0)
		return;
	OpenMesh::IO::Options opt;
	opt += OpenMesh::IO::Options::VertexTexCoord;

	if (show_option_ == ORIGINAL) {
		OpenMesh::IO::write_mesh(mesh_, fname, opt);
	}

	else if (show_option_ == SLICED) {
		OpenMesh::IO::write_mesh(sliced_mesh_, fname, opt);
	}


	
}

void OTEViewer::UpdateMeshData(SurfaceMesh &mesh)
{
	OpenMeshToMatrix(mesh, V_,V_normal_, F_, F_normal_);
	TC_.resize(F_.rows(), 3);
	TC_.col(0) = Eigen::VectorXd::Constant(TC_.rows(), 1.0);
	TC_.col(1) = Eigen::VectorXd::Constant(TC_.rows(), 1.0);
	TC_.col(2) = Eigen::VectorXd::Constant(TC_.rows(), 1.0);
	if (V_.rows() == 0) return;
	data.clear();
	data.set_mesh(V_, F_);
	data.set_colors(TC_);
	data.set_normals(V_normal_);
	data.set_normals(F_normal_);
}


void OTEViewer::UpdateTextureCoordData(SurfaceMesh &mesh)
{
	OpenMeshCoordToMatrix(mesh, UV_);
	data.set_uv(UV_);
	if (R_.rows() > 0)
		data.set_texture(R_, G_, B_);
}


void OTEViewer::ShowUV()
{
	if (UV_.rows() == 0) return;
	UV_Z0_.resize(UV_.rows(), 3);
	UV_Z0_.setZero();
	UV_Z0_.block(0, 0, UV_.rows(), UV_.cols()) = UV_;
	data.clear();
	data.set_mesh(UV_Z0_, F_);
	data.set_uv(UV_);
	data.set_colors(TC_);
	data.set_texture(R_, G_, B_);

	Eigen::MatrixXd V_normal = V_normal_;
	V_normal.col(2) = (V_normal.col(2)).cwiseAbs();
	Eigen::MatrixXd F_normal = F_normal_;
	F_normal.col(2) = F_normal.col(2).cwiseAbs();
	data.set_normals(V_normal);
	data.set_normals(F_normal);
}

void OTEViewer::ShowCoveringSpace()
{
	if (euclidean_) {
		EuclideanCoveringSpaceComputer covering_computer(sliced_mesh_, cone_vts_);
		covering_computer.Compute();
		covering_computer.GenerateMeshMatrix(V_, V_normal_, F_, F_normal_);
		data.clear();
		data.set_mesh(V_, F_);
		V_normal_.col(2) = (V_normal_.col(2)).cwiseAbs();
		F_normal_.col(2) = (F_normal_.col(2)).cwiseAbs();
		data.set_normals(V_normal_);
		data.set_normals(F_normal_);
		data.set_uv(V_.block(0, 0, V_.rows(),2));
	}

	if (hyperbolic_) {
		HyperbolicCoveringSpaceComputer covering_computer(sliced_mesh_, cone_vts_);
		covering_computer.Compute();
		covering_computer.GenerateMeshMatrix(V_, V_normal_, F_, F_normal_);
		data.clear();
		data.set_mesh(V_, F_);
		V_normal_.col(2) = (V_normal_.col(2)).cwiseAbs();
		F_normal_.col(2) = (F_normal_.col(2)).cwiseAbs();
		data.set_normals(V_normal_);
		data.set_normals(F_normal_);
		data.set_uv(V_.block(0, 0, V_.rows(), 2));
	}

}

void OTEViewer::ShowHalfedges(SurfaceMesh &mesh, std::vector<OpenMesh::HalfedgeHandle> h_vector)
{
	// This algorithm to generate cylinder is from the library libhedra
	using namespace OpenMesh;

	if (h_vector.size() == 0) return;

	Eigen::MatrixXd lineV, lineTC;
	Eigen::MatrixXi lineT;
	Eigen::MatrixXd P1, P2;
	Eigen::MatrixXd lineColors;
	double radius = 0.005;

	HalfedgesToMatrix(mesh, h_vector, P1, P2);
	lineColors.resize(P1.rows(), 3);
	lineColors.setZero();
	lineColors.col(0) = Eigen::VectorXd::Constant(P1.rows(), 1);
	LineCylinders(
		P1,
		P2,
		radius, lineColors,
		10,
		false,
		lineV,
		lineT,
		lineTC);

	MergeMeshMatrix(V_, F_, TC_, lineV, lineT, lineTC, V_, F_, TC_);

	
	data.clear();
	data.set_mesh(V_, F_);
	data.set_colors(TC_);
}

void OTEViewer::ShowBoundaries(SurfaceMesh &mesh)
{
	mesh.RequestBoundary();
	auto boundaries = mesh.GetBoundaries();
	ShowHalfedges(mesh,boundaries.front());
	
}

void OTEViewer::ShowSliceAndCones()
{

	// For edges;
	Eigen::MatrixXd lineV, lineTC;
	Eigen::MatrixXi lineT;
	Eigen::MatrixXd P1, P2;
	Eigen::MatrixXd lineColors;

	// For vertices;
	Eigen::MatrixXd P;
	Eigen::MatrixXi T;
	Eigen::MatrixXd V, TC;
	Eigen::MatrixXd vertexColors;


	marker_.GenerateMatrix(P, P1, P2);
	
	if (P1.rows() > 0) {

		lineColors.resize(P1.rows(), 3);
		lineColors.setZero();
		lineColors.col(0) = Eigen::VectorXd::Constant(P1.rows(), 1);


		LineCylinders(
			P1,
			P2,
			0.005, lineColors,
			10,
			false,
			lineV,
			lineT,
			lineTC);
		MergeMeshMatrix(V_, F_, TC_, lineV, lineT, lineTC, V_, F_, TC_);
	}

	if (P.rows() > 0) {

		vertexColors.resize(P.rows(), 3);
		vertexColors.setZero();
		vertexColors.col(1) = Eigen::VectorXd::Constant(P.rows(), 1);


		PointSpheres(P,
			0.03,
			vertexColors,
			10,
			false,
			V,
			T,
			TC);
		MergeMeshMatrix(V_, F_, TC_, V, T, TC, V_, F_, TC_);
	}
	//show
	data.clear();
	data.set_mesh(V_, F_);
	data.set_colors(TC_);
}

void OTEViewer::ShowVertexLabels()
{
	
	using namespace OpenMesh;
	if (mesh_.n_edges() == 0) return;
	for (SurfaceMesh::VertexIter viter = mesh_.vertices_begin(); viter != mesh_.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec3d p = mesh_.point(v);
		data.add_label(Eigen::Vector3d(p[0], p[1], p[2]), std::to_string(v.idx()));
	}
	
}



void OTEViewer::UpdateMeshViewer()
{
	if (show_option_ == ORIGINAL) {
		UpdateMeshData(mesh_);
		UpdateTextureCoordData(mesh_);
		if(show_slice_)
			ShowSliceAndCones();
		if (show_vertex_labels_)
			ShowVertexLabels();
		ShowSelction();
	}
	else if (show_option_ == SLICED) {
		UpdateMeshData(sliced_mesh_);
		UpdateTextureCoordData(sliced_mesh_);
		if (show_boundaries_) {
			ShowBoundaries(sliced_mesh_);
		}
	}
	else if (show_option_ == EMBEDDING) {
		UpdateMeshData(sliced_mesh_);
		UpdateTextureCoordData(sliced_mesh_);
		ShowUV();
	}

	else if (show_option_ == COVERING_SPACE) {
		ShowCoveringSpace();
	}

}

void OTEViewer::SetSlice()
{
	using namespace OpenMesh;
	VertexHandle v1 = selected_verts_.front();
	VertexHandle v2 = selected_verts_.back();
	marker_.ComputeAndSetSlice(v1, v2);
	selected_verts_.clear();
}

void OTEViewer::SetSingularity(double cone_angle)
{
	using namespace OpenMesh;
	SurfaceMesh &mesh = mesh_;
	VertexHandle v = selected_verts_.front();
	marker_.SetSingularity(v, cone_angle);
	selected_verts_.clear();
}

void OTEViewer::LoadMarker()
{
	std::string fname = igl::file_dialog_open();
	marker_.LoadFromFile(fname);
}

void OTEViewer::SaveMarker()
{
	std::string fname = igl::file_dialog_save();
	marker_.SaveToFile(fname);
}

void OTEViewer::FindIntersection(double x, double y)
{
	using namespace OpenMesh;
	using namespace Eigen;
	SurfaceMesh &mesh = mesh_;

	Eigen::Matrix4f model = core.model;
	Eigen::Matrix4f view = core.view;
	Eigen::Matrix4f proj = core.proj;

	Eigen::Matrix4f world = view * model;
	
	double min_dist = 10000;
	VertexHandle nearest_vertex;
	for (auto viter = mesh.vertices_begin(); viter != mesh.vertices_end(); ++viter) {
		VertexHandle v = *viter;
		Vec3d p = mesh.point(v);
		Vec3d n = mesh.normal(v);
		Eigen::Vector4f p_h(p[0], p[1], p[2], 1.);
		Eigen::Vector4f n_h(n[0], n[1], n[2], 0.);
		
		Eigen::Vector4f current_n = world * n_h;
		Eigen::Vector4f current_p = world * p_h;
		Eigen::Vector4f eye_direction(-current_p(0), -current_p(1), 5 - current_p(2), 0);
		
		if ((current_n.transpose() * eye_direction)(0) < 0) continue;

		Eigen::Vector4f p_proj = proj * current_p;
		p_proj /= p_proj(3);
		
		Vec2f mouse(x, y);
		Vec2f pp(p_proj(0), p_proj(1));

		double dist = (mouse - pp).norm();
		
		if (dist < min_dist) {
			min_dist = dist;
			nearest_vertex = v;
		}
		
	}
	bool exist = false;
	for (auto it = selected_verts_.begin(); it != selected_verts_.end(); ++it) {
		if (nearest_vertex == *it) {
			selected_verts_.erase(it);
			exist = true;
			break;
		}
	}
	if(!exist)
		selected_verts_.push_back(nearest_vertex);
}

void OTEViewer::ShowSelction()
{
	using namespace OpenMesh;
	// For edges;
	Eigen::MatrixXd bigV, bigTC;
	Eigen::MatrixXi bigT;

	// For vertices;
	Eigen::MatrixXd P;
	Eigen::MatrixXi T;
	Eigen::MatrixXd V, TC;
	Eigen::MatrixXd vertexColors;


	bigV = V_;
	bigT = F_;
	bigTC = TC_;

	P.resize(selected_verts_.size(), 3);
	
	if (P.rows() > 0) {
		int i = 0;
		for (auto it = selected_verts_.begin(); it != selected_verts_.end(); ++it, ++i) {
			VertexHandle v = *it;
			Vec3d p = mesh_.point(v);
			P.row(i) = Eigen::Vector3d(p[0], p[1], p[2]);
		}

		vertexColors.resize(P.rows(), 3);
		vertexColors.setZero();
		vertexColors.col(1) = Eigen::VectorXd::Constant(P.rows(), 1);


		PointSpheres(P,
			0.03,
			vertexColors,
			10,
			false,
			V,
			T,
			TC);
		MergeMeshMatrix(bigV, bigT, bigTC, V, T, TC, bigV, bigT, bigTC);
		data.clear();
		data.set_mesh(bigV, bigT);
		data.set_colors(bigTC);
	}
	

}
