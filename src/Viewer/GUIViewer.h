#ifndef OTE_VIEWER
#define OTE_VIEWER

#include <iostream>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/png/texture_from_png.h>

#include <Eigen/Dense>

#include <MeshDefinition.h>

#include <MeshFormConverter.h>
#include <LineCylinder.h>
#include <MeshMerger.h>
#include <PointSphere.h>

#include <EuclideanOrbifoldSolver.h>
#include <HyperbolicOrbifoldSolver.h>
#include <BFF.h>

#include <EuclideanCoveringSpace.h>
#include <HyperbolicCoveringSpace.h>

#include <MeshMarker.h>


enum ShowOption { ORIGINAL, SLICED, EMBEDDING, COVERING_SPACE };

class OTEViewer: public igl::opengl::glfw::Viewer
{
public: 
	OTEViewer();
	void Init();
private:
	SurfaceMesh mesh_;
	SurfaceMesh sliced_mesh_;
	MeshMarker marker_;

	std::vector<OpenMesh::VertexHandle> cone_vts_;

	Eigen::MatrixXd V_;
	Eigen::MatrixXd V_normal_;
	Eigen::MatrixXd UV_Z0_; // z equal to zero
	Eigen::MatrixXi F_;
	Eigen::MatrixXd F_normal_;
	Eigen::MatrixXd TC_;
	Eigen::MatrixXd UV_;

	

	// texture RGB channels;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R_;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G_;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B_;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A_;
	
	// Flags
	ShowOption show_option_ = ORIGINAL;
	bool show_boundaries_ = false;
	bool show_slice_ = false;
	bool show_vertex_labels_ = false;
	bool euclidean_ = false;
	bool hyperbolic_ = false;

	double cone_angle_ = 0.;

	bool selection_mode_ = false;


	std::list<OpenMesh::VertexHandle> selected_verts_;


protected:

	// Init functions
	void InitMenu();
	void InitKeyboard();
	void InitMouse();

	// IO functions
	void LoadMesh();
	void LoadTexture();
	void SaveMesh();

	void UpdateMeshData(SurfaceMesh &mesh);  
	void UpdateTextureCoordData(SurfaceMesh &mesh);
	void ShowUV();
	void ShowCoveringSpace();

	void ShowHalfedges(SurfaceMesh &mesh, std::vector<OpenMesh::HalfedgeHandle> h_vector);
	void ShowBoundaries(SurfaceMesh &mesh);
	void ShowSliceAndCones();
	void ShowVertexLabels();
	void UpdateMeshViewer();


	// Setting slices and singularities
	void SetSlice();
	void SetSingularity(double cone_angle);
	
	void LoadMarker();
	void SaveMarker();

	// These two functions are used to achieve vertex selection.
	void FindIntersection(double x, double y);
	void ShowSelction();

};

#endif