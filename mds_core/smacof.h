#include<iostream>
#include<time.h>
#include<Eigen/dense>
#include "my_graph.h"

using namespace std;
using namespace Eigen;

double smacof_stress(VectorXd X_v, MatrixXd del_m, MatrixXd W_m);
double smacof_stress_2(MatrixXd X_m, MatrixXd del_m);
void smacof_process_2(bool geodesic, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph, MatrixXd& result_2);
void smacof_process(bool geodesic, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph, VectorXd& result);
Vector2f smacof_dynamic(int ind, Vector2f pre_dir, vector<Vector2i>& pts, VectorXd& result);
MatrixXd my_pinv(MatrixXd& X_m);
