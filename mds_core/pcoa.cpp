#include "pcoa.h"

void pcoa_process(int dim, bool geodesic, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph, VectorXd& result, MatrixXd& result_2)
{
	MatrixXd del_m(pts.size(), pts.size());
	if(geodesic)
	{
		build_graph(pts, my_graph);
		show_graph(pts, my_graph);
		cal_del(del_m, pts, my_graph);
	}
	else
	{
		for(int i = 0 ; i < pts.size(); i++)
		{
			for(int j = 0; j < pts.size(); j++)
			{
				del_m(i,j) = sqrt((double)(pts[i] - pts[j]).squaredNorm());
			}
		}
	}

	MatrixXd Q_m(pts.size(), pts.size());
	for(int i = 0 ; i < pts.size(); i++)
	{
		for(int j = 0; j < pts.size(); j++)
		{
			Q_m(i,j) = -0.5 * del_m(i,j) * del_m(i,j);
			//Q_m(i,j) = -0.5 * (abs(pts[i](0) - pts[j](0)) + abs(pts[i](1) - pts[j](1))) * (abs(pts[i](0) - pts[j](0)) + abs(pts[i](1) - pts[j](1)));
		}
	}

	MatrixXd H_m = MatrixXd::Identity(pts.size(), pts.size());
	VectorXd one_v = VectorXd::Ones(pts.size());
	double ratio = 1.0 / pts.size();
	H_m = H_m - ratio * one_v * one_v.transpose();

	MatrixXd B_m = H_m * Q_m * H_m;

	SelfAdjointEigenSolver<MatrixXd> eigensolver(B_m);
	VectorXd value = eigensolver.eigenvalues();
	MatrixXd matrix = eigensolver.eigenvectors();
	if(dim == 1)
	{
		result = matrix.col(matrix.cols()-1) * sqrt(value(value.rows()-1));
	}
	else
	{
		vector<VectorXd> tem;
		tem.resize(2);
		result_2 = MatrixXd(pts.size(),2);
		tem[0] = matrix.col(matrix.cols()-1) * sqrt(value(value.rows()-1));
		tem[1] = matrix.col(matrix.cols()-2) * sqrt(value(value.rows()-2));
		for(int i = 0; i < pts.size(); i++)
		{
			for(int j = 0; j < 2; j++)
				result_2(i,j) = tem[j](i);
		}
	}
}
