#include "smacof.h"

double smacof_stress(VectorXd X_v, MatrixXd del_m, MatrixXd W_m)
{
	double stress = 0;
	for(int i = 0; i < X_v.rows(); i++)
	{
		for(int j = i+1; j < X_v.rows(); j++)
		{
			stress += W_m(i,j) * (abs(X_v(i) - X_v(j)) - del_m(i,j)) * (abs(X_v(i) - X_v(j)) - del_m(i,j));
		}
	}

	return stress;
}

MatrixXd my_pinv(MatrixXd& X_m)
{
	/*
	MatrixXd one_m = MatrixXd::Ones(X_m.rows(), X_m.cols());
	double tem_d = (1.0 / (double)(X_m.rows() * X_m.rows()));
	MatrixXd X_inv = (X_m + one_m).inverse() - tem_d * one_m;
	*/
	
	JacobiSVD<MatrixXd> svd(X_m, ComputeThinU | ComputeThinV);
	VectorXd m_sigma_inv = svd.singularValues(); 
	double tol = 1e-5;
	for(int i = 0; i < m_sigma_inv.rows(); i++)
	{
		if(m_sigma_inv(i) > tol)
			m_sigma_inv(i) = 1.0 / m_sigma_inv(i);
		else
			m_sigma_inv(i) = 0;
	}

	MatrixXd X_inv = svd.matrixV() * (m_sigma_inv.asDiagonal()) * (svd.matrixU().transpose());
	/*cout<<svd.matrixV().rows()<<" "<<svd.matrixV().cols()<<endl;
	cout<<m_sigma_inv.asDiagonal().rows()<<" "<<m_sigma_inv.asDiagonal().cols()<<endl;
	cout<<svd.matrixU().rows()<<" "<<svd.matrixU().cols()<<endl;*/

	return X_inv;
}

void smacof_process(bool geodesic, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph, VectorXd& result)
{
	srand(time(0));
	result = VectorXd(pts.size());
	for(int i = 0; i < result.rows(); i++)
	{
		//result(i) = rand();
		result(i) = pts[i](1);
	}
	
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

	//cout<<del_m<<endl;

	for(int k = 0; k < 100; k++)
	{
		double first_stress = smacof_stress(result, del_m, MatrixXd::Ones(pts.size(), pts.size()));
		MatrixXd B_m(pts.size(), pts.size());
		for(int i = 0; i < pts.size(); i++)
		{
			double sum = 0;
			for(int j = 0; j < pts.size(); j++)
			{
				if(i == j)
					B_m(i,j) = 0;
				else if(result(i) - result(j) == 0)
				{
					B_m(i,j) = 0;
				}
				else
				{
					B_m(i,j) = -del_m(i,j) / abs(result(i) - result(j));
				}
				sum += -B_m(i,j);
			}
			B_m(i,i) = sum;
		}
		//cout<<B_m<<endl;

		VectorXd tem_v= result;
		result = (1.0 / pts.size()) * B_m * tem_v;
		double second_stress = smacof_stress(result, del_m, MatrixXd::Ones(pts.size(), pts.size()));

		cout<<first_stress<<" "<<second_stress<<endl;

		if(abs(first_stress - second_stress) < 1e-6)
			break;
	}
	cout<<endl;
}

Vector2f smacof_dynamic(int ind, Vector2f pre_dir, vector<Vector2i>& pts, VectorXd& result)
{
	vector<Vector3i> pts_check;
	pts_check.resize(pts.size());
	for(int i = 0; i < pts.size(); i++)
	{
		pts_check[i] = Vector3i(pts[i](0), pts[i](1), 0);
	}

	vector<int> my_neighbor;
	my_neighbor.push_back(ind);

	int num;
	if(pts.size() >= 4)
		num = 4;	
	else
		num = pts.size();

	for(int j = 0; j < num; j++)
	{
		Vector2f tem_v = find_min(pts_check, ind, pts);
		my_neighbor.push_back(tem_v(0));
	}

	srand(time(0));
	VectorXd result_v = VectorXd(my_neighbor.size());
	for(int i = 0; i < result_v.rows(); i++)
	{
		result_v(i) = rand();
		//result(i) = pts[my_neighbor[i]](1);
	}
	
	MatrixXd del_m(my_neighbor.size(), my_neighbor.size());
	for(int i = 0 ; i < my_neighbor.size(); i++)
	{
		for(int j = 0; j < my_neighbor.size(); j++)
		{
			del_m(i,j) = sqrt((double)(pts[my_neighbor[i]] - pts[my_neighbor[j]]).squaredNorm());
		}
	}

	MatrixXd W_m(my_neighbor.size(), my_neighbor.size());
	for(int i = 0; i < my_neighbor.size(); i++)
	{
		for(int j = 0; j < my_neighbor.size(); j++)
		{
			Vector2f p;
			if(i == j)
				p = Vector2f(0,0);
			else
			{
				p = Vector2f(pts[my_neighbor[i]](0) - pts[my_neighbor[j]](0), pts[my_neighbor[i]](1) - pts[my_neighbor[j]](1)).normalized();
			}

			double dot_value = Vector2f(p(0),p(1)).dot(pre_dir);
			W_m(i,j) = abs(dot_value);// * exp(-del_m(i,j)/100.0);
		}
	}
	//cout<<W_m<<endl<<endl;

	MatrixXd V_m(my_neighbor.size(), my_neighbor.size());
	for(int i = 0; i < my_neighbor.size(); i++)
	{
		double sum = 0;
		for(int j = 0; j < my_neighbor.size(); j++)
		{
			if(i == j)
				V_m(i,j) = 0;
			else
				V_m(i,j) = -W_m(i,j);
			sum +=  -V_m(i,j);
		}
		V_m(i,i) = sum;
	}

	MatrixXd V_inv = my_pinv(V_m);

	for(int k = 0; k < 100; k++)
	{
		double first_stress = smacof_stress(result_v, del_m, W_m);
		MatrixXd B_m(my_neighbor.size(), my_neighbor.size());
		for(int i = 0; i < my_neighbor.size(); i++)
		{
			double sum = 0;
			for(int j = 0; j < my_neighbor.size(); j++)
			{
				if(i == j)
					B_m(i,j) = 0;
				else if(result_v(i) - result_v(j) == 0)
				{
					B_m(i,j) = 0;
				}
				else
				{
					B_m(i,j) = -del_m(i,j) * W_m(i,j) / abs(result_v(i) - result_v(j));
				}
				sum += -B_m(i,j);
			}
			B_m(i,i) = sum;
		}
		//cout<<B_m<<endl;

		VectorXd tem_v= result_v;
		result_v = V_inv * B_m * tem_v;
		double second_stress = smacof_stress(result_v, del_m, W_m);

		cout<<first_stress<<" "<<second_stress<<endl;

		if(abs(first_stress - second_stress) < 1e-6)
			break;
	}
	cout<<endl;

	int pre = 0;
	double dis = FLT_MAX;
	for(int i = 1; i < result_v.rows(); i++)
	{
		if(result_v(i) < result_v(0) & abs(result_v(i) - result_v(0)) < dis)
		{
			pre = i;
			dis = abs(result_v(i) - result_v(0));
		}
	}

	int  next = 0;
	dis = FLT_MAX;
	for(int i = 1; i < result_v.rows(); i++)
	{
		if(result_v(i) > result_v(0) & abs(result_v(i) - result_v(0)) < dis)
		{
			next = i;
			dis = abs(result_v(i) - result_v(0));
		}
	}

	cout<<my_neighbor[pre]<<" "<<my_neighbor[next]<<endl;
	Vector2f dir;
	if(next == 0)
		dir = Vector2f(pts[pre](0) - pts[ind](0), pts[pre](1) - pts[ind](1)).normalized();
	else if(pre == 0)
		dir = Vector2f(pts[ind](0) - pts[next](0), pts[ind](1) - pts[next](1)).normalized();
	else
	{
		dir = Vector2f(pts[pre](0) - pts[ind](0), pts[pre](1) - pts[ind](1)).normalized() + Vector2f(pts[ind](0) - pts[next](0), pts[ind](1) - pts[next](1)).normalized();
		dir = dir / 2.0;
		dir.normalize();
	}

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

	return dir;
}

double smacof_stress_2(MatrixXd X_m, MatrixXd del_m)
{
	double stress = 0;

	for(int i = 0; i < del_m.rows(); i++)
	{
		for(int j = i+1; j < del_m.rows(); j++)
		{
			double dis = sqrt((X_m(i,0) - X_m(j,0)) * (X_m(i,0) - X_m(j,0)) + (X_m(i,1) - X_m(j,1)) * (X_m(i,1) - X_m(j,1)));
			stress += (dis - del_m(i,j)) * (dis - del_m(i,j));
		}
	}

	return stress;
}

void smacof_process_2(bool geodesic, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph, MatrixXd& result_2)
{
	srand(time(0));
	result_2 = MatrixXd(pts.size(),2);
	if(geodesic)
	{
		for(int i = 0; i < pts.size(); i++)
		{
			result_2(i,0) = pts[i](0);
			result_2(i,1) = pts[i](1);
		}
	}
	else
	{
		for(int i = 0; i < pts.size(); i++)
		{
			for(int j = 0; j < 2; j++)
				result_2(i,j) = rand();
		}
	}
	
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

	//cout<<del_m<<endl;
	int num;
	if(geodesic)
		num = 5;
	else
		num = 1000;
	
	for(int k = 0; k < num; k++)
	{
		double first_stress = smacof_stress_2(result_2, del_m);
		MatrixXd B_m(pts.size(), pts.size());
		for(int i = 0; i < pts.size(); i++)
		{
			double sum = 0;
			for(int j = 0; j < pts.size(); j++)
			{
				double dis = sqrt((result_2(i,0) - result_2(j,0)) * (result_2(i,0) - result_2(j,0)) + (result_2(i,1) - result_2(j,1)) * (result_2(i,1) - result_2(j,1)));
				if(i == j)
					B_m(i,j) = 0;
				else if(dis == 0)
				{
					B_m(i,j) = 0;
				}
				else
				{
					B_m(i,j) = -del_m(i,j) / dis;
				}
				sum += -B_m(i,j);
			}
			B_m(i,i) = sum;
		}

		MatrixXd tem_m = result_2;
		result_2 = (1.0 / pts.size()) * B_m * tem_m;
		double second_stress = smacof_stress_2(result_2, del_m);

		cout<<first_stress<<" "<<second_stress<<endl;

		if(abs(first_stress - second_stress) < 1e-3)
			break;
	}
	cout<<endl;
}
