#ifndef OTE_VIEWER
#define OTE_VIEWER
#include <iostream>
#include <igl/viewer/Viewer.h>
#include <Eigen/Dense>
#include <nanogui/formhelper.h>
#include <nanogui/screen.h>

class OTEViewer
{
public: 
	OTEViewer();
private:
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	igl::viewer::Viewer viewer;
	bool init()
	{
		// Add new group
		viewer.ngui->addGroup("New Group");

		// Add a button
		viewer.ngui->addButton("Print Hello", []() { std::cout << "Hello\n"; });

		// Generate menu
		viewer.screen->performLayout();

		return false;
	};
};

#endif