#include "GUIViewer.h"

OTEViewer::OTEViewer()
{
}

void OTEViewer::Init()
{
	InitMenu();
	InitKeyboard();
}


void OTEViewer::InitMenu()
{
	callback_init = [this](igl::viewer::Viewer& viewer) 
	{
		viewer.ngui->addWindow(Eigen::Vector2i(220, 0), "Control Panel");
		
		viewer.ngui->addGroup("Mesh IO");
		viewer.ngui->addButton("Load Mesh", [this]() {this->LoadMesh(); });
		viewer.ngui->addButton("Load Texture", []() {printf("Load Texture"); });
		viewer.ngui->addButton("Save Mesh", []() {printf("Save Mesh"); });
		
		viewer.ngui->addGroup("Core Functions");
		viewer.ngui->addButton("Euclidean Orbifold Embedding", []() {printf("Euclidean Orbifold Embedding"); });

		viewer.screen->performLayout();
		return false; 
	};

	
}

void OTEViewer::InitKeyboard()
{
	callback_key_down = [](igl::viewer::Viewer& viewer, unsigned char key, int modifier) 
	{
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
	UpdateMeshData(mesh_);
	UpdateTextureCoordData(mesh_);
	ShowMesh();

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
	OpenMesh::IO::write_mesh(mesh_, fname, opt);
}

void OTEViewer::UpdateMeshData(SurfaceMesh &mesh)
{
	OpenMeshToMatrix(mesh, V_, F_);
}


void OTEViewer::UpdateTextureCoordData(SurfaceMesh &mesh)
{
	OpenMeshCoordToMatrix(mesh, UV_);
	UV_Z0_.resize(UV_.rows(), 3);
	UV_Z0_.setZero();
	UV_Z0_.block(0, 0, UV_.rows(), UV_.cols()) = UV_;
	data.set_uv(UV_);
}

void OTEViewer::ShowMesh()
{
	data.set_mesh(V_, F_);
}

void OTEViewer::ShowUV()
{
	data.set_mesh(UV_Z0_, F_);
}
