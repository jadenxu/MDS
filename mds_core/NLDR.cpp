#include "NLDR.h"

Vector2f find_near(Vector2i cur, vector<Vector2i>& pts)
{
	Vector2f tem(-1,FLT_MAX);
	for(int i = 0; i < pts.size(); i++)
	{
		double dis = sqrt((double)(cur - pts[i]).squaredNorm());
		if(dis < tem(1))
		{
			tem(0) = i;
			tem(1) = dis;
		}
	}

	return tem;
}

void my_lle(vector<Vector2i>& pts, VectorXd& result)
{
	vector<Vector3i> pts_check;
	pts_check.resize(pts.size());
	for(int i = 0; i < pts.size(); i++)
	{
		pts_check[i] = Vector3i(pts[i](0), pts[i](1), 0);
	}

	vector<vector<int> > my_neighbor;
	my_neighbor.resize(pts.size());

	int num;
	if(pts.size() >= 2)
		num = 2;	
	else
		num = pts.size();

	for(int i = 0; i < pts.size(); i++)
	{
		for(int j = 0; j < num; j++)
		{
			Vector2f tem_v = find_min(pts_check, i, pts);
				my_neighbor[i].push_back(tem_v(0));
		}

		for(int j = 0; j < num; j++)
		{
			pts_check[my_neighbor[i][j]](2) = 0;
		}
	}

	MatrixXd W_m = MatrixXd::Zero(pts.size(), pts.size());
	
	for(int i = 0; i < pts.size(); i++)
	{
		MatrixXd dis_m(2, num);
		for(int j = 0; j < num; j++)
		{
			Vector2i tem = pts[my_neighbor[i][j]] - pts[i];
			dis_m(0,j) = tem(0);
			dis_m(1,j) = tem(1);
		}

		MatrixXd cov_m = dis_m.transpose() * dis_m;
		while(cov_m.determinant() == 0)
		{
			cov_m += (1e-3 * MatrixXd::Identity(num,num));
		}
		VectorXd w_v = cov_m.inverse() * VectorXd::Ones(num);

		double w_sum = 0;
		for(int j = 0; j < num; j++)
			w_sum += w_v(j);
		for(int j = 0; j < 2; j++)
		{
			W_m(i, my_neighbor[i][j]) = w_v(j) / w_sum;
		}
	}

	MatrixXd tem_m = MatrixXd::Identity(pts.size(),pts.size()) - W_m;
	MatrixXd M = tem_m.transpose() * tem_m;

	SelfAdjointEigenSolver<MatrixXd> eigensolver(M);
	VectorXd value = eigensolver.eigenvalues();
	MatrixXd matrix = eigensolver.eigenvectors();

	result = matrix.col(1);
}

Vector2d local_pca(int ind, vector<Vector2i>& pts, VectorXd& result)
{
	vector<Vector3i> pts_check;
	pts_check.resize(pts.size());
	for(int i = 0; i < pts.size(); i++)
	{
		pts_check[i] = Vector3i(pts[i](0), pts[i](1), 0);
	}

	vector<int> my_neighbor;
	Vector2d center(pts[ind](0), pts[ind](1));
	my_neighbor.push_back(ind);

	int num;
	if(pts.size() >= 4)
		num = 4;	
	else
		num = pts.size();

	MatrixXd X_m(num+1,2);
	X_m(0,0) = pts[ind](0);
	X_m(0,1) = pts[ind](1);

	for(int j = 0; j < num; j++)
	{
		Vector2f tem_v = find_min(pts_check, ind, pts);
		my_neighbor.push_back(tem_v(0));
		X_m(j+1,0) = pts[tem_v(0)](0);
		X_m(j+1,1) = pts[tem_v(0)](1);
		center(0) += pts[tem_v(0)](0);
		center(1) += pts[tem_v(0)](1); 
	}
	
	center = center / my_neighbor.size();

	for(int j = 0; j < X_m.rows(); j++)
	{
		X_m(j,0) -= center(0);
		X_m(j,1) -= center(1);
	}

	MatrixXd C_m = X_m.transpose() * X_m;
	C_m *= (1.0) / X_m.rows();

	SelfAdjointEigenSolver<MatrixXd> eigensolver(C_m);
	VectorXd value = eigensolver.eigenvalues();
	MatrixXd matrix = eigensolver.eigenvectors();
	//cout<<matrix<<endl<<endl;

	VectorXd result_v = X_m * matrix.col(matrix.cols()-1);
	//cout<<matrix.col(matrix.cols()-1)<<endl;
	result = VectorXd::Ones(pts.size());
	result = 1.5 * result;
	double ma = result_v(0);
	double mi = result_v(0);
	for(int i = 1; i < result_v.rows(); i++)
	{
		if(result_v(i) < mi)
			mi = result_v(i);

		if(result_v(i) > ma)
			ma = result_v(i);
	}

	double range = ma - mi;

	for(int i = 0; i < my_neighbor.size(); i++)
	{
		result[my_neighbor[i]] = (result_v(i)-mi) / range;
	}

	return matrix.col(matrix.cols()-1);
}

void dynamic_pca(int ind, Vector2f& pre, vector<Vector2i>& pts, VectorXd& result)
{
	vector<Vector3i> pts_check;
	pts_check.resize(pts.size());
	for(int i = 0; i < pts.size(); i++)
	{
		pts_check[i] = Vector3i(pts[i](0), pts[i](1), 0);
	}

	vector<int> my_neighbor;
	Vector2d center(pts[ind](0), pts[ind](1));
	my_neighbor.push_back(ind);

	int num;
	if(pts.size() >= 4)
		num = 4;	
	else
		num = pts.size();

	MatrixXd X_m(num+1,2);
	X_m(0,0) = pts[ind](0);
	X_m(0,1) = pts[ind](1);

	for(int j = 0; j < num; j++)
	{
		Vector2f tem_v = find_min(pts_check, ind, pts);
		my_neighbor.push_back(tem_v(0));
		X_m(j+1,0) = pts[tem_v(0)](0);
		X_m(j+1,1) = pts[tem_v(0)](1);
		center(0) += pts[tem_v(0)](0);
		center(1) += pts[tem_v(0)](1); 
	}
	
	center = center / (double)my_neighbor.size();

	for(int j = 0; j < X_m.rows(); j++)
	{
		X_m(j,0) -= center(0);
		X_m(j,1) -= center(1);
	}

	MatrixXd C_m = X_m.transpose() * X_m;
	C_m *= (1.0) / X_m.rows();

	SelfAdjointEigenSolver<MatrixXd> eigensolver(C_m);
	VectorXd value = eigensolver.eigenvalues();
	VectorXd b_v = value(value.rows()-1) / 2.0 * Vector2d(-pre(0),-pre(1));
	MatrixXd C_inv = (C_m - MatrixXd::Identity(C_m.rows(), C_m.cols())).inverse();
	VectorXd dir = C_inv * b_v;
	dir.normalize();
	pre = Vector2f(dir(0),dir(1));

	VectorXd result_v = X_m * dir;
	//cout<<matrix.col(matrix.cols()-1)<<endl;
	result = VectorXd::Ones(pts.size());
	result = 1.5 * result;
	double ma = result_v(0);
	double mi = result_v(0);
	for(int i = 1; i < result_v.rows(); i++)
	{
		if(result_v(i) < mi)
			mi = result_v(i);

		if(result_v(i) > ma)
			ma = result_v(i);
	}

	double range = ma - mi;

	for(int i = 0; i < my_neighbor.size(); i++)
	{
		result[my_neighbor[i]] = (result_v(i)-mi) / range;
	}
}
