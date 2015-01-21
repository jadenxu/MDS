#ifndef MY_QLABEL_H
#define MY_QLABEL_H

#include <QLabel>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QMessageBox>

#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<Eigen/dense>
#include<iostream>

#include "mds_core/pcoa.h"
#include "mds_core/smacof.h"
#include "mds_core/NLDR.h"

using namespace std;
using namespace Eigen;
using namespace cv;

class my_qlabel : public QLabel
{
	Q_OBJECT

public:
	my_qlabel(QWidget *parent);
	~my_qlabel();
	void show_image();
	void show_graph();
	void draw_circle(Point x, double num);
	void my_clear();
	void show_result(int dim);
	void draw_origin();
	void mousePressEvent(QMouseEvent *ev);
	void mouseMoveEvent(QMouseEvent *ev);
	void mouseReleaseEvent(QMouseEvent *ev);
	void pcoa(int dim, bool geodesic);
	void smacof(int dim, bool geodesic);
	void lle();
	void pca_local();
	void enable_pca_local();
	void smacof_local();
	void enable_smacof_local();
	int pos;
	bool pca_l;
	bool smacof_l;

private:
	Mat image;
	bool draw;
	vector<Vector2i> pts;
	vector<vector<Vector2f> > my_graph;
	VectorXd result;
	MatrixXd result_2;
	Vector2f my_dir;
};

#endif // MY_QLABEL_H
