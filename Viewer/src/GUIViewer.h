#ifndef OTE_VIEWER
#define OTE_VIEWER

#include <iostream>

#include <igl/viewer/Viewer.h>
#include <igl/png/texture_from_file.h>

#include <nanogui/formhelper.h>
#include <nanogui/screen.h>

#include <Eigen/Dense>

#include <MeshDefinition.h>

#include <MeshFormConverter.h>
#include <LineCylinder.h>
#include <MeshMerger.h>

#include <EuclideanOrbifoldSolver.h>



enum ShowOption { ORIGINAL, SLICED, EMBEDDING };

class OTEViewer: public igl::viewer::Viewer
{
public: 
	OTEViewer();
	void Init();
private:
	SurfaceMesh mesh_;
	SurfaceMesh sliced_mesh_;

	Eigen::MatrixXd V_;
	Eigen::MatrixXd UV_Z0_; // z equal to zero
	Eigen::MatrixXi F_;
	Eigen::MatrixXd TC_;
	Eigen::MatrixXd UV_;

	// texture RGB channels;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R_;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> G_;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> B_;
	Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> A_;
	
	// Init functions
	void InitMenu();
	void InitKeyboard();

	// IO functions
	void LoadMesh();
	void LoadTexture();
	void SaveMesh();

	void UpdateMeshData(SurfaceMesh &mesh);  
	void UpdateTextureCoordData(SurfaceMesh &mesh);
	void ShowUV();

	void ShowHalfedges(SurfaceMesh &mesh, std::vector<OpenMesh::HalfedgeHandle> h_vector);
	void ShowBoundaries(SurfaceMesh &mesh);

	void UpdateMeshViewer();

protected: // Flags
	ShowOption show_option_ = ORIGINAL;
	bool show_boundaries_ = false;
};

#endif