#include "my_qlabel.h"

my_qlabel::my_qlabel(QWidget *parent)
	: QLabel(parent)
{
	Mat new_img(500, 500, CV_8UC3);
	new_img = Scalar(255,255,255);
	image = new_img;
	show_image();

	draw = false;
	pca_l = false;
	smacof_l = false;
	pos = -1;
	my_dir = Vector2f(-1,-1);
}

void my_qlabel::show_image()
{
	Mat tem;
	cvtColor(image,tem,CV_BGR2RGB);
	QImage img;
	img = QImage(tem.data, tem.cols, tem.rows, tem.cols*3, QImage::Format_RGB888);
	this->setPixmap(QPixmap::fromImage(img));
	this->resize(this->pixmap()->size());
}

void my_qlabel::draw_circle(Point x, double num)
{
	int thickness = -1;
	int lineType = 8;
	int radius = 5;
	Mat tem_hsv(1,1, CV_8UC3);
	Mat tem_bgr;
	tem_hsv.at<Vec3b>(0,0) = Vec3b(num * 90, 240, 240);
	cvtColor(tem_hsv, tem_bgr, CV_HSV2BGR);
	circle(image, x , radius , Scalar(tem_bgr.at<Vec3b>(0,0)[0],tem_bgr.at<Vec3b>(0,0)[1],tem_bgr.at<Vec3b>(0,0)[2]), thickness, lineType);

	show_image();
}

void my_qlabel::draw_origin()
{
	Mat new_img(500, 500, CV_8UC3);
	new_img = Scalar(255,255,255);
	image = new_img;

	for(int i = 0; i < pts.size(); i++)
	{
		circle(image, Point(pts[i](0), pts[i](1)), 5, Scalar(0,255,0), -1, 8);
	}
	show_image();
}

void my_qlabel::mousePressEvent(QMouseEvent* ev)
{
	draw = true;
	draw_circle(Point(ev->x(), ev->y()), 2.0/3.0);
	pts.push_back(Vector2i(ev->x(),ev->y()));
}

void my_qlabel::mouseMoveEvent(QMouseEvent *ev)
{
	if(draw)
	{
		draw_circle(Point(ev->x(), ev->y()), 2.0/3.0);
		pts.push_back(Vector2i(ev->x(),ev->y()));
	}
	/*else if(local)
	{
		Vector2f tem = find_near(Vector2i(ev->x(), ev->y()), pts);
		if(tem(1) < 5)
		{
			pca_local(tem(0));
		}
		else
		{
			draw_origin();
		}
	}*/
}

void my_qlabel::mouseReleaseEvent(QMouseEvent* ev)
{
	draw = false;
}

void my_qlabel::pcoa(int dim, bool geodesic)
{
	pcoa_process(dim, geodesic, pts, my_graph, result, result_2);
	if(dim == 1)
		show_result(1);
	else
		show_result(2);
}

void my_qlabel::smacof(int dim, bool geodesic)
{
	if(dim == 1)
	{
		smacof_process(geodesic, pts, my_graph, result);
		show_result(1);
	}
	else
	{
		smacof_process_2(geodesic, pts, my_graph, result_2);
		show_result(2);
	}
}

void my_qlabel::lle()
{
	my_lle(pts, result);
	show_result(1);
}

void my_qlabel::pca_local()
{
	if(pos < 0)
		pos += pts.size();
	else if(pos >= pts.size())
		pos -= pts.size();

	/*if(my_dir == Vector2f(-1,-1))
	{
		Vector2d tem = local_pca(0, pts, result);
		my_dir = Vector2f(tem(0),tem(1));
	}

	dynamic_pca(pos, my_dir, pts, result);*/

	Vector2d dir = local_pca(pos, pts, result);

	Mat new_img(500, 500, CV_8UC3);
	new_img = Scalar(255,255,255);
	image = new_img;
	for(int i = 0; i < pts.size(); i++)
	{
		draw_circle(Point(pts[i](0), pts[i](1)), result(i));
	}
	line(image, Point(pts[pos](0), pts[pos](1)), Point(pts[pos](0) + 10, pts[pos](1) + 10.0 * (dir(1)/dir(0))), Scalar(0,0,0), 2, CV_AA);
	show_image();
}

void my_qlabel::enable_pca_local()
{
	pca_l = true;
	pos = 0;
	pca_local();
}

void my_qlabel::enable_smacof_local()
{
	smacof_l = true;
	pos = 0;
	smacof_local();
}

void my_qlabel::smacof_local()
{
	if(pos < 0)
		pos += pts.size();
	else if(pos >= pts.size())
		pos -= pts.size();

	if(my_dir == Vector2f(-1,-1))
	{
		Vector2d tem = local_pca(0, pts, result);
		my_dir = Vector2f(tem(0),tem(1));
	}

	my_dir = smacof_dynamic(pos, my_dir, pts, result);

	Mat new_img(500, 500, CV_8UC3);
	new_img = Scalar(255,255,255);
	image = new_img;
	for(int i = 0; i < pts.size(); i++)
	{
		draw_circle(Point(pts[i](0), pts[i](1)), result(i));
	}
	//line(image, Point(pts[pos](0), pts[pos](1)), Point(pts[pos](0) + 20, pts[pos](1) + 20.0 * (my_dir(1)/my_dir(0))), Scalar(0,0,0), 2, CV_AA);
	show_image();
}

void my_qlabel::show_result(int dim)
{
	if(dim == 1)
	{
		double max = result(0);
		double min = result(0);
		for(int i = 1; i < result.rows(); i++)
		{
			if(result(i) < min)
			{
				min = result(i);
			}

			if(result(i)> max)
			{
				max = result(i);
			}
		}
		

		double range = max - min;
		for(int i = 0; i < result.rows(); i++)
		{
			result(i) = (result(i) - min) / range;
		}

		Mat new_img(500, 500, CV_8UC3);
		new_img = Scalar(255,255,255);
		image = new_img;
		for(int i = 0; i < pts.size(); i++)
		{
			draw_circle(Point(pts[i](0), pts[i](1)), result(i));
		}
		show_image();
	}
	else
	{
		for(int k = 0; k < 2; k++)
		{	
			double max = result_2(0,k);
			double min = result_2(0,k);
			for(int i = 1; i < pts.size(); i++)
			{
				if(result_2(i,k) < min)
				{
					min = result_2(i,k);
				}

				if(result_2(i,k) > max)
				{
					max = result_2(i,k);
				}
			}

			double range = max - min;
			for(int i = 0; i < pts.size(); i++)
			{
				result_2(i,k) = (result_2(i,k) - min) / range;
			}
		}
		Mat new_img(500, 500, CV_8UC3);
		new_img = Scalar(255,255,255);
		for(int i = 0; i < pts.size(); i++)
		{
			int thickness = -1;
			int lineType = 8;
			int radius = 5;
			circle(new_img, Point(result_2(i,0) * 400 + 50, result_2(i,1) * 400 + 50) , radius , Scalar(0,255,0), thickness, lineType);
		}
		imshow("2D-CASE", new_img);

		for(int i =0; i < pts.size(); i++)
		{
			pts[i] = Vector2i(result_2(i,0) * 400, result_2(i,1) * 400);
		}
	}
}

void my_qlabel::my_clear()
{
	Mat new_img(500, 500, CV_8UC3);
	new_img = Scalar(255,255,255);
	image = new_img;
	show_image();

	pts.clear();
	my_graph.clear();
	draw = false;
	pca_l = false;
	smacof_l = false;
	pos = -1;
	my_dir = Vector2f(-1,-1);
}

my_qlabel::~my_qlabel()
{

}
