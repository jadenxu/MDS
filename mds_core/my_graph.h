#ifndef MY_GRAPH
#define MY_GRAPH

#include "my_node.h"
#include<queue>
#include<iostream>
#include<Eigen/dense>
#include<opencv2/core/core.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>

using namespace cv;
using namespace std;
using namespace Eigen;

Vector2f find_min(vector<Vector3i>& pts_check, int ind, vector<Vector2i>& pts);
bool check_graph(int position, int insert, vector<vector<Vector2f> >& my_graph);
void my_dijkstra(MatrixXd& del_m, int ind, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph);
void build_graph(vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph);
void show_graph(vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph);
void cal_del(MatrixXd& del_m, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph);

#endif
