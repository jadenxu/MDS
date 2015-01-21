#ifndef MDS_H
#define MDS_H

#include <QtGui/QMainWindow>
#include "ui_mds.h"
#include "my_qlabel.h"

class MDS : public QMainWindow
{
	Q_OBJECT

public:
	MDS(QWidget *parent = 0, Qt::WFlags flags = 0);
	~MDS();
	void create_actions();

protected:
	void keyPressEvent(QKeyEvent* ev);

public slots:
	void my_forward1();
	void my_restart();
	void my_forward2();

private:
	Ui::MDSClass ui;
};

#endif // MDS_H
