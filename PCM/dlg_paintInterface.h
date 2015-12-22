#ifndef _DLG_PAINT_INTERFACE_H
#define _DLG_PAINT_INTERFACE_H

#include <QtWidgets/QDialog>
#include "basic_types.h"

#include "ui_dlg_paintInterface.h"
#include "paint_canvas.h"
#include "paint_interface.h"

class PaintUi : public QDialog
{
	Q_OBJECT
public:
	PaintUi()
	{
		ui.setupUi(this);
		canvas == NULL;
		cur_sample = -1;
	}

	void init(PaintCanvas* canvas_,IndexType cur_sample_id_)
	{
		canvas = canvas_;
		cur_sample = cur_sample_id_;

		connect(ui.init, SIGNAL( clicked() ), this, SLOT(initStroke() ) );
        connect(ui.doPaint, SIGNAL( clicked() ), this, SLOT(drawStrokes() ));
		
	};


public slots: 

	void initStroke()
	{	
		if (canvas->single_operate_tool_)
		{
			delete canvas->single_operate_tool_;
		}
		canvas->single_operate_tool_ = new PaintStroke(canvas);

		canvas->single_operate_tool_->set_tool_type(Tool::SELECT_STROKE);

		canvas->single_operate_tool_->set_cur_smaple_to_operate(cur_sample);

		//canvas->updateGL();
	}

	void drawStrokes()
	{
		PaintStroke* paintInter = new PaintStroke(canvas);
		paintInter->draw();
		canvas->updateGL();
	}
private:
	Ui::PaintInterface ui;
	IndexType cur_sample;
	PaintCanvas* canvas;
};

#endif