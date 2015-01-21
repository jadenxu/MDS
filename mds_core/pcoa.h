#include<iostream>
#include<Eigen/dense>
#include "my_graph.h"

using namespace std;
using namespace Eigen;

void pcoa_process(int dim, bool geodesic, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph, VectorXd& result, MatrixXd& result_2);
