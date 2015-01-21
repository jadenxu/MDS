#include<iostream>
#include<vector>
#include<Eigen/dense>
#include "my_graph.h"

using namespace std;
using namespace Eigen;

void my_lle(vector<Vector2i>& pts, VectorXd& result);
Vector2d local_pca(int ind, vector<Vector2i>& pts, VectorXd& result);
Vector2f find_near(Vector2i cur, vector<Vector2i>& pts);
void dynamic_pca(int ind, Vector2f& pre, vector<Vector2i>& pts, VectorXd& result);
