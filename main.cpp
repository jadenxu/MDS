#include "mds.h"
#include <QtGui/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	MDS w;
	w.show();
	return a.exec();
}
