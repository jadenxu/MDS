#include"my_graph.h"

Vector2f find_min(vector<Vector3i>& pts_check, int ind,vector<Vector2i>& pts)
{
	Vector2f tem_v(-1, FLT_MAX);
	for(int i = 0; i < pts_check.size(); i++)
	{
		if(i == ind)
			continue;
		
		double dis = sqrt((double)(pts[i] - pts[ind]).squaredNorm());
		if(dis < tem_v(1) && pts_check[i](2) == 0)
		{
			tem_v(0) = i;
			tem_v(1) = dis;
		}
	}

	pts_check[tem_v(0)](2) = 1;
	return tem_v;
}

bool check_graph(int position, int insert, vector<vector<Vector2f> >& my_graph)
{
	bool exist = false;

	for(int i = 0; i < my_graph[position].size(); i++)
	{
		if(my_graph[position][i](0) == insert)
		{
			exist = true;
			break;
		}
	}

	return exist;
}

void my_dijkstra(MatrixXd& del_m, int ind, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph)
{
	vector<my_node> nodes;
	nodes.resize(pts.size());
	for(int i = 0; i < pts.size(); i++)
	{
		nodes[i] = my_node(i,1e30,false);
	}

	nodes[ind].value = 0;
	priority_queue<my_node, vector<my_node>, cmp> q;
	q.push(nodes[ind]);

	while(!q.empty())
	{
		my_node tem = q.top();
		q.pop();

		if(tem.visit)
			continue;

		tem.visit = true;
		del_m(ind,tem.index) = tem.value;

		for(int i = 0; i < my_graph[tem.index].size(); i++)
		{
			int c_ind = my_graph[tem.index][i](0);
			if(nodes[c_ind].value > tem.value + my_graph[tem.index][i](1))
			{
				nodes[c_ind].value = tem.value + my_graph[tem.index][i](1);
				q.push(nodes[c_ind]);
			}
		}
	}
}

void build_graph(vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph)
{
	vector<Vector3i> pts_check;
	pts_check.resize(pts.size());
	for(int i = 0; i < pts.size(); i++)
	{
		pts_check[i] = Vector3i(pts[i](0), pts[i](1), 0);
	}

	my_graph.resize(pts.size());
	int num;
	if(pts.size() >= 4)
	{
		num = 4;	
	}
	else
	{
		num = pts.size();
	}

	for(int i = 0; i < pts.size(); i++)
	{
		
		for(int j = 0; j < 4; j++)
		{
			Vector2f tem_v = find_min(pts_check, i, pts);
			if(!check_graph(i, tem_v(0), my_graph))
				my_graph[i].push_back(tem_v);

			Vector2f tem_v1(i, tem_v(1));
			if(!check_graph(tem_v(0),i, my_graph))
				my_graph[tem_v(0)].push_back(tem_v1);
		}

		for(int j = 0; j < num; j++)
		{
			pts_check[my_graph[i][j](0)](2) = 0;
		}
	}
}

void show_graph(vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph)
{
	Mat new_img(500, 500, CV_8UC3);
	new_img = Scalar(255,255,255);

	for(int i = 0; i < pts.size(); i++)
	{
		circle(new_img, Point(pts[i](0),pts[i](1)), 5, Scalar(0,255,0), -1, 8);
	}

	for(int i = 0; i < my_graph.size(); i++)
	{
		for(int j = 0; j < my_graph[i].size(); j++)
		{
			line(new_img, Point(pts[i](0), pts[i](1)), Point(pts[my_graph[i][j](0)](0), pts[my_graph[i][j](0)](1)), Scalar(0,255,255), 1, CV_AA);
		}
	}

	imshow("ok",new_img);
}

void cal_del(MatrixXd& del_m, vector<Vector2i>& pts, vector<vector<Vector2f> >& my_graph)
{
	for(int i = 0; i < pts.size(); i++)
	{
		my_dijkstra(del_m, i, pts, my_graph);
	}
}
