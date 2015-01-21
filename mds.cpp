#include "mds.h"

MDS::MDS(QWidget *parent, Qt::WFlags flags)
	: QMainWindow(parent, flags)
{
	ui.setupUi(this);

	create_actions();
	setCentralWidget(ui.label);
	this->setWindowTitle("MDS");
	this->resize(500,500);
}

void MDS::create_actions()
{
	connect(ui.actionForward, SIGNAL(triggered()), this, SLOT(my_forward1()));
	connect(ui.actionForward1, SIGNAL(triggered()), this, SLOT(my_forward2()));
	connect(ui.actionRestart, SIGNAL(triggered()), this, SLOT(my_restart()));
}

void MDS::my_forward1()
{
	ui.label->pcoa(1, true);
}

void MDS::my_forward2()
{
	//ui.label->smacof(2, false);
	//ui.label->lle();
	ui.label->enable_pca_local();
	//ui.label->enable_smacof_local();
}

void MDS::keyPressEvent(QKeyEvent* ev)
{
	
	if(ui.label->pca_l)
	{
		switch(ev->key())
		{
			case Qt::Key_Left:
				ui.label->pos--;
				ui.label->pca_local();
				break;
			case Qt::Key_Right:
				ui.label->pos++;
				ui.label->pca_local();
				break;
		}
	}
	else if(ui.label->smacof_l)
	{
		switch(ev->key())
		{
			case Qt::Key_Left:
				ui.label->pos--;
				ui.label->smacof_local();
				break;
			case Qt::Key_Right:
				ui.label->pos++;
				ui.label->smacof_local();
				break;
		}
	}
}

void MDS::my_restart()
{
	ui.label->my_clear();
}

MDS::~MDS()
{

}
