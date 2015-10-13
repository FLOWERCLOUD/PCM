 
 
 #ifndef _DLG_FITPLAN_H
 #define _DLG_FITPLAN_H

 #include <QtWidgets/QDialog>
 #include "ui_dlg_plan.h"
 #include "basic_types.h"

 #include "plans_fit.h"
 
 class PlanFitUI : public QDialog
 {
   //Q_OBJECT
 public:
     PlanFitUI(IndexType _selectSmpId ):selectSmpId(_selectSmpId)
 	{
 		ui.setupUi(this);
 	}
 
 	void init()
 	{
 		connect(ui.smpLev, SIGNAL( textEdited( const QString&)), this, SLOT(getResolution(const QString&) ));  // getInitalValue(const QString&)) 
 		
		connect(ui.modelTime, SIGNAL(textEdited( const QString&) ), this, SLOT( getModels(const QString&) ));

		connect(ui.hierThershold, SIGNAL(valueChanged(int)), this, SLOT( updateThreshold(int)));

 		connect(ui.pushButton, SIGNAL(clicked (bool) ),this, SLOT(runT(bool)) );
 	
 	}
 
 public slots:
 public:
 	void runT(bool b)
 	{
 		PlanClassifier* pFit = new PlanClassifier(selectSmpId);
 
 		pFit->setParem(selectSmpId,m_reso,m_models,m_tradeOff);
 
 		connect(pFit,SIGNAL(finished()) ,pFit, SLOT(deleteLater() ));
 
 		pFit->start();
 	}
 
 	void getInitalValue(const QString& txt)
 	{
 		m_tradeOff = txt.toFloat();
 	}

	void updateThreshold(int _thereshold)
	{
		m_tradeOff = (0.01) *_thereshold;
		ui.ThresholdValue->setText(QString::number(m_tradeOff));
	}
 
	void getModels(const QString& txt)
	{
		m_models = txt.toInt();
	}

	void getResolution(const QString& txt)
	{
		m_reso = txt.toInt();
	}

 private:
 	ScalarType m_tradeOff;
 
 	IndexType selectSmpId;

	IndexType m_models;
	IndexType m_reso;
 
 	Ui::Planfit ui;
 
 };
 #endif
 
 
 
