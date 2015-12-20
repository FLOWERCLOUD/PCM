#include "main_window.h"
#include <QMessageBox>
#include <QFileDialog>
#include <QLabel>
#include <QStatusBar>
#include <QSettings>
#include <QCloseEvent>
#include <QPlainTextEdit>
#include <QAbstractItemModel>
#include <QStandardItemModel>
#include <QGroupBox>
#include <QColorDialog>
#include <QComboBox>
#include "file_system.h"
#include "sample.h"
#include "sample_set.h"
#include "file_io.h"
#include "color_table.h"
#include "time.h"
#include "trajectory_classifier.h"
#include "tracer.h"
#include "DeformaleRegistration.h"
#include "spectral_clustering.h"

#include "sample_properity.h"
using namespace qglviewer;
using namespace std;

#pragma optimize( "", off )

main_window::main_window(QWidget *parent)
	: QMainWindow(parent),
	cur_select_sample_idx_(-1),
	last_select_sample_idx_(-1)
{


	ui.setupUi(this);

	QGLFormat format = QGLFormat::defaultFormat();
	format.setSampleBuffers(true);
	format.setSamples(8);
	main_canvas_ = new PaintCanvas(format, this);

	setCentralWidget(main_canvas_);

	setWindowTitle("PCM");

	setContextMenuPolicy(Qt::CustomContextMenu);
//	setWindowState(Qt::WindowMaximized);

	setFocusPolicy(Qt::ClickFocus);

	main_canvas_->camera()->setType(Camera::PERSPECTIVE);

	//Codes for initializing treeWidget
	QStringList headerLabels;
	headerLabels <<"Index"<< "Name" << "Label" << "#Point";
	ui.treeWidget->setHeaderLabels(headerLabels);
	ui.treeWidget->setColumnWidth( 0,50 );
	connect(ui.treeWidget, SIGNAL(itemClicked( QTreeWidgetItem * , int  ) ),
		this, SLOT(selectedSampleChanged( QTreeWidgetItem * , int  )) );


	createAction();
	createStatusBar();

	if(m_linkageUi || m_graphCutUi || m_planFitUi)
	{
		m_linkageUi  = NULL;
		m_graphCutUi = NULL;
		m_planFitUi = NULL;
	}

	//srand(time(NULL));

}


void main_window::resetSampleSet()
{
	cur_import_files_attr_.clear();
	cur_select_sample_idx_ = last_select_sample_idx_ = -1;
	SampleSet::get_instance().clear();
}

void main_window::createAction()
{
	createFileMenuAction();
	createPaintSettingAction();
	createAlgorithmAction();
	createToolAction();

}


void main_window::createAlgorithmAction()
{
	connect(ui.actionClustering, SIGNAL(triggered()), this, SLOT(doClustering()));
	connect(ui.actionRegister,SIGNAL(triggered()),this,SLOT(doRegister()));
	//connect(ui.actionSpectral_Cluster,SIGNAL(triggered()),this,SLOT(doSpectralCluster()));
	connect(ui.actionGraphCut,SIGNAL(triggered()),this,SLOT(doGraphCut()));
	connect(ui.actionCalculateNorm,SIGNAL(triggered()),this,SLOT(computeSampleNormal()));
	connect(ui.actionClusterAll,SIGNAL(triggered()),this,SLOT(batchTrajClustering()));
	connect(ui.actionVisDistortion,SIGNAL(triggered()),this,SLOT(visDistor()));

	connect(ui.actionGCopti,SIGNAL(triggered()),this,SLOT(doGCOptimization()));

	connect(ui.actionPlanFit,SIGNAL(triggered()), this ,SLOT(doPlanFit() ) );

	connect(ui.actionOrderLabels,SIGNAL(triggered() ),this,SLOT(doOrder()) );
	connect(ui.actionRefineSigFrame, SIGNAL(triggered() ), this, SLOT( doRefineSigFrame() ));
}

void main_window::createPaintSettingAction()
{
	connect(ui.actionSet_Visible, SIGNAL(triggered()), this, SLOT(setSampleVisible()) );
	connect(ui.actionSet_Invisible, SIGNAL(triggered()), this, SLOT(setSampleInvisible()));
	connect(ui.actionObject_Color, SIGNAL(triggered()), this, SLOT(setObjectColorMode()));
	connect(ui.actionVertex_Color, SIGNAL(triggered()), this, SLOT(setVertexColorMode()));
	connect(ui.actionLabel_Color, SIGNAL(triggered()), this, SLOT(setLabelColorMode()));
	connect(ui.actionShow_Tracjectory, SIGNAL(triggered()), this, SLOT(showTracer()));
	connect(ui.actionDont_Trace,SIGNAL(triggered()), this, SLOT(clearTracer()));
}

void main_window::createToolAction()
{
	connect( ui.actionSelect_Mode, SIGNAL(triggered()), this, SLOT(setSelectToolMode()) );
	connect( ui.actionScene_Mode, SIGNAL(triggered()),this, SLOT(setSceneToolMode()));
}

void main_window::setObjectColorMode()
{
	main_canvas_->which_color_mode_ = PaintCanvas::OBJECT_COLOR;
	main_canvas_->updateGL();
}

void main_window::setVertexColorMode()
{
	main_canvas_->which_color_mode_ = PaintCanvas::VERTEX_COLOR;
	main_canvas_->updateGL();
}

void main_window::setLabelColorMode()
{
	main_canvas_->which_color_mode_ = PaintCanvas::LABEL_COLOR;
	main_canvas_->updateGL();
}

void main_window::setSelectToolMode()
{
	if (cur_select_sample_idx_==-1)
	{
		return;
	}

	if (main_canvas_->single_operate_tool_)
	{
		delete main_canvas_->single_operate_tool_;
	}
	main_canvas_->single_operate_tool_ = new SelectTool(main_canvas_);
	main_canvas_->single_operate_tool_->set_tool_type(Tool::SELECT_TOOL);
	main_canvas_->single_operate_tool_->set_cur_smaple_to_operate(cur_select_sample_idx_);

	main_canvas_->updateGL();

}

void main_window::setSceneToolMode()
{
	main_canvas_->single_operate_tool_->set_tool_type(Tool::EMPTY_TOOL);
	main_canvas_->updateGL();
}

void main_window::createFileMenuAction()
{
	connect(ui.actionImportFiles, SIGNAL(triggered() ),this, SLOT( openFiles() )  );
	connect(ui.actionSaveSnapshot, SIGNAL(triggered() ), this, SLOT(saveSnapshot() ) );
}

// void main_window::createFileMenuAction()
// {
// 	connect(ui.actionImportFiles, SIGNAL(triggered()),this, SLOT(openFiles()));
// }

bool main_window::openFile()
{
	QString fileName = QFileDialog::getOpenFileName(this,
		tr("Import point cloud from file"), ".",
		tr(
		"Ascii Point cloud (*.lab *.xyz *.pwn *.pcd)\n"
		"All files (*.*)")
		);

	if (fileName.isEmpty())
		return false;

	return true;
}


bool main_window::setSampleVisible()
{
	SampleSet::get_instance()[cur_select_sample_idx_].set_visble(true);
	main_canvas_->updateGL();
	return true;
}

bool main_window::setSampleInvisible()
{
	SampleSet::get_instance()[cur_select_sample_idx_].set_visble(false);
	main_canvas_->updateGL();
	return true;
}

void main_window::selectedSampleChanged(QTreeWidgetItem * item, int column)
{
	last_select_sample_idx_ = cur_select_sample_idx_;
	cur_select_sample_idx_ = item->text( 0 ).toInt();

	//change the active frame
	if (m_linkageUi)
	{
	  m_linkageUi->setCurrentSmp(cur_select_sample_idx_);//
	}

	if (m_graphCutUi)
	{
	  m_graphCutUi->setCurrentSmp(cur_select_sample_idx_);
	}


	if (last_select_sample_idx_!= -1)
	{
		SampleSet::get_instance()[last_select_sample_idx_].set_selected(false);
	}
	if ( cur_select_sample_idx_ != -1)
	{
		SampleSet::get_instance()[cur_select_sample_idx_].set_selected(true);
	}
	main_canvas_->updateGL();

}



void main_window::showTracer()
{
	main_canvas_->setTracerShowOrNot(true);
	main_canvas_->updateGL();
}

void main_window::clearTracer()
{
	main_canvas_->setTracerShowOrNot(false);
	Tracer::get_instance().clear_records();
	main_canvas_->updateGL();
}

bool main_window::openFiles()
{
	QString dir = QFileDialog::getExistingDirectory(this,tr("Import point cloud files"),".");
	if (dir.isEmpty())
		return false;

	resetSampleSet();

	QDir file_dir(dir);
	if ( !file_dir.exists() )
	{
		return false;
	}
	file_dir.setFilter(QDir::Files);

	QFileInfoList file_list = file_dir.entryInfoList();
	IndexType sample_idx = 0;
	for (IndexType file_idx = 0; file_idx < file_list.size(); file_idx++)
	{
		QFileInfo file_info = file_list.at(file_idx);
		FileIO::FILE_TYPE file_type;

		if (file_info.suffix() == "xyz")
		{
			file_type = FileIO::XYZ;
		}
		else if(file_info.suffix() == "ply")
		{
			file_type = FileIO::PLY;
		}
		else
		{
			continue;
		}

		string file_path = file_info.filePath().toStdString();
		cur_import_files_attr_.push_back( make_pair(FileSystem::base_name(file_path), 
													FileSystem::extension(file_path)) );

		Sample* new_sample = FileIO::load_point_cloud_file(file_path, file_type,sample_idx);
		if (new_sample != nullptr)
		{
			SampleSet::get_instance().push_back(new_sample);
			sample_idx++;
		}

	}


	createTreeWidgetItems();

	return true;
}

void main_window::createStatusBar()
{
	coord_underMouse_label_ = new QLabel(this);
	coord_underMouse_label_->setAlignment(Qt::AlignLeft);

	vtx_idx_underMouse_label_ = new QLabel(this);
	coord_underMouse_label_->setAlignment(Qt::AlignRight);
	
	ui.statusBar->addWidget( coord_underMouse_label_ , 1 );
	ui.statusBar->addWidget( vtx_idx_underMouse_label_, 0 );
}

void main_window::createTreeWidgetItems()
{
	ui.treeWidget->clear();
	
	SampleSet& set = SampleSet::get_instance();
	for ( int sample_idx=0; sample_idx < set.size(); sample_idx++ )
	{
		QTreeWidgetItem* item = new QTreeWidgetItem(ui.treeWidget); 


		ColorType color = set[ sample_idx ].color();

		item->setData(0, Qt::DisplayRole, sample_idx);
		item->setData(1,Qt::DecorationRole, QColor(color(0)*255, color(1)*255, color(2)*255) );
		item->setData(2, Qt::DisplayRole, set[sample_idx].num_vertices() );

		ui.treeWidget->insertTopLevelItem(sample_idx, item);
	}
}

void main_window::showCoordinateAndIndexUnderMouse( const QPoint& point )
{
	//Mouse point info come from canvas
	bool found = false;
    qglviewer::Vec v = main_canvas_->camera()->pointUnderPixel(point, found);
	if ( !found )
	{
		v = qglviewer::Vec();
	}
	QString coord_str = QString("XYZ = [%1, %2, %3]").arg(v.x).arg(v.y).arg(v.z);
	coord_underMouse_label_->setText(coord_str);

	IndexType idx;
	IndexType label;
	if ( !found || cur_select_sample_idx_==-1 )
	{
		idx = -1;
		label = -1;
	}
	else
	{
		Sample& cur_selected_sample = SampleSet::get_instance()[cur_select_sample_idx_];
		Vec4 v_pre(v.x - Paint_Param::g_step_size(0)*cur_select_sample_idx_,
			v.y - Paint_Param::g_step_size(1)*cur_select_sample_idx_,
			v.z - Paint_Param::g_step_size(2)*cur_select_sample_idx_ ,1.);
		//Necessary to do this step, convert view-sample space to world-sample space
		v_pre = cur_selected_sample.matrix_to_scene_coord().inverse() * v_pre;
		idx = cur_selected_sample.closest_vtx( PointType(v_pre(0), v_pre(1), v_pre(2)) );
		label = cur_selected_sample[idx].label();
	}
	QString idx_str = QString("VERTEX INDEX = [%1],LABEL = [%2]").arg(idx).arg(label);
	vtx_idx_underMouse_label_->setText( idx_str );

	return;
}

void main_window::doClustering()
{

// 	TrajectoryClassifier* classifier = new TrajectoryClassifier(cur_select_sample_idx_);
// 	connect(classifier, SIGNAL(finish_compute()), this, SLOT(finishClustering()));
// 	connect(classifier, SIGNAL(finished()), classifier, SLOT(deleteLater() ));
// 	classifier->start();

// 	JLinkageUI *dlg = new JLinkageUI(cur_select_sample_idx_);
// 	dlg->init();
// 	dlg->show();

	m_linkageUi = new JLinkageUI(cur_select_sample_idx_);
	m_linkageUi->init();
	m_linkageUi->show();



}
void main_window::doRegister()
{
	DeformableRegistration* register_ = new DeformableRegistration();
	connect(register_,SIGNAL(finish_compute()),this,SLOT(finishRegister()));
	connect(register_,SIGNAL(finished()),register_,SLOT(deleteLater()));
	register_->start();
}
void main_window::doGraphCut()
{
	GraphNodeCtr* graphCut = new GraphNodeCtr();
	connect(graphCut, SIGNAL(finish_compute()), this, SLOT(finishGraphCut()));
	connect(graphCut,SIGNAL(finished()),graphCut,SLOT(deleteLater()));
	graphCut->start();
}

void main_window::doGCOptimization()
{
// 	GCop* gc= new GCop();
// 	connect(gc,SIGNAL(finish_compute()),this,SLOT(finishGCotimization()));
// 	connect(gc,SIGNAL(finished()),gc,SLOT(deleteLater()));
// 	gc->start();

 	m_graphCutUi = new GraphCutUI(cur_select_sample_idx_);
 	m_graphCutUi->init();
 	m_graphCutUi->show();

}

void main_window::doPlanFit()
{
// 	PlanClassifier* planF = new PlanClassifier(cur_select_sample_idx_);
// 	connect(planF, SIGNAL(finish_compute()), this, SLOT(finishDoPlanFit()) );
// 	connect(planF, SIGNAL(finished()),planF,SLOT(deleteLater()));
// 	planF->start();

	m_planFitUi = new PlanFitUI(cur_select_sample_idx_);

	m_planFitUi->init();

	m_planFitUi->show();

}

void main_window::doOrder()
{
		GCop* gc= new GCop();
		gc->setCurSmpId(cur_select_sample_idx_);
		gc->orderLabelsOnly();
}

void main_window::doRefineSigFrame()
{
			GCop* gc= new GCop();
			gc->setCurSmpId(cur_select_sample_idx_);
			gc->refineSegm();
}

// void main_window::doSpectralCluster()
// {
// 	SpectralClusteringThread* specCla = new SpectralClusteringThread();
// 	connect(specCla,SIGNAL(finish_compute()),this,SLOT(finishSpectralCluster()));
// 	connect(specCla,SIGNAL(finished()),specCla,SLOT(deleteLater()));
// 	specCla->start();
// }
void main_window::finishClustering()
{

}
void main_window::finishRegister()
{

}
void main_window::finishGraphCut()
{

}

void main_window::finishGCotimization()
{

}
// void main_window::finishDoPlanFit()
// {
// 
// }

// void main_window::finishSpectralCluster()
// {
// 
// }
main_window::~main_window()
{
	SampleSet::get_instance().clear();


}
void main_window::computeSampleNormal()
{
	auto camera_look_at = main_canvas_->camera()->viewDirection();
	//SampleManipulation::computerMinMax(cur_select_sample_idx_);
	SampleManipulation::compute_normal( cur_select_sample_idx_ ,NormalType(-camera_look_at.x, -camera_look_at.y, -camera_look_at.z));
	//SampleManipulation::compute_normal( cur_select_sample_idx_ ,NormalType(-1.0,0.0,0.0));
}

void main_window::batchTrajClustering()
{
	iterate_sample_idx_ = 1;
	//iterate_sample_idx_ = 100;
	iterateTrajClustering();
}

void main_window::visDistor()
{
		TrajectoryClassifier* classifier = new TrajectoryClassifier(cur_select_sample_idx_);
		classifier->visDistor();
}

void main_window::iterateTrajClustering()
{
	if ( iterate_sample_idx_>=SampleSet::get_instance().size() )
	{
		return ;
	}
	TrajectoryClassifier* cluster = new TrajectoryClassifier(iterate_sample_idx_++);
	connect( cluster, SIGNAL(finished()), cluster, SLOT(deleteLater()) );
	connect(cluster,SIGNAL(finish_compute()), this, SLOT(iterateTrajClustering()));
	cluster->start();
}

bool main_window::saveSnapshot()
{
	SaveSnapshotDialog dialog(this);

	dialog.setValues(getCanvas()->ss);

	if (dialog.exec()==QDialog::Accepted)
	{
	getCanvas()->ss=dialog.getValues();
	getCanvas()->saveSnapshot();

	// if user ask to add the snapshot to raster layers
	/*
	if(dialog.addToRasters())
	{
	  QString savedfile = QString("%1/%2%3.png")
	.arg(GLA()->ss.outdir).arg(GLA()->ss.basename)
	.arg(GLA()->ss.counter,2,10,QChar('0'));

	  importRaster(savedfile);
	}
	*/
	return true;
	}

	return false;
}