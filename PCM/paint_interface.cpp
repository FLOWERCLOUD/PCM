#include "paint_interface.h"
#include "sample_set.h"
#include "windows.h"
#include <gl/gl.h>
#include <gl/glu.h>
#include "paint_canvas.h"
#include "globals.h"
#include "color_table.h"

void PaintStroke::move(QMouseEvent *e)
{
	if (left_mouse_button_ == true)
	{
		drag(e);
	}

}

void PaintStroke::drag(QMouseEvent *e)
{

	rectangle_.setBottomRight( e->pos() );

	//
	crtSketch.push_back(e->pos() );

	canvas_->updateGL();

}

void PaintStroke::release(QMouseEvent *e)
{
	if (left_mouse_button_ == true)
	{
		// Possibly swap left/right and top/bottom to make rectangle_ valid.
		rectangle_ = rectangle_.normalized();
		//rectangle_ = QRect(e->pos(), e->pos());//delete current rect and propre to draw next rect;
		canvas_->updateGL();
		left_mouse_button_ = false;
	}

}

void PaintStroke::press(QMouseEvent* e)
{
	if (e->button() == Qt::LeftButton)
	{
		left_mouse_button_ = true;
		rectangle_ = QRect(e->pos(), e->pos());

		//
		crtSketch.push_back(e->pos() );

		canvas_->updateGL();
	}

}

void PaintStroke::draw_stroke()
{
	QPainter painter;
	if (crtSketch.size() == 0) return;
	QBrush thisBrush(QColor(255, 0, 0, 100));
	double strokeWidth;
	strokeWidth = 2;
	QPen thisPen(thisBrush, strokeWidth, Qt::SolidLine, Qt::RoundCap);
	thisPen.setJoinStyle(Qt::RoundJoin);
	painter.save();
	painter.setPen(thisPen);
	QPainterPath thisPath;
	thisPath.moveTo(crtSketch[0].x(), crtSketch[0].y());
	for (int i = 1; i < crtSketch.size(); ++i)
	{
		thisPath.lineTo(crtSketch[i].x(), crtSketch[i].y() );
	}
	painter.drawPath(thisPath);
	painter.restore();
}

void PaintStroke::draw_rectangle()
{
	canvas_->startScreenCoordinatesSystem();

	glDisable(GL_LIGHTING);

	glLineWidth(2.0);
	glColor4f(0.0f, 1.0f, 1.0f, 0.5f);
	glBegin(GL_LINE_LOOP);
	glVertex2i(rectangle_.left(),  rectangle_.top());
	glVertex2i(rectangle_.right(), rectangle_.top());
	glVertex2i(rectangle_.right(), rectangle_.bottom());
	glVertex2i(rectangle_.left(),  rectangle_.bottom());
	glEnd();	

	//draw a mask

	// 	glEnable(GL_BLEND);
	// 	glDepthMask(GL_FALSE);
	// 	glColor4f(0.0, 0.0, 0.4f, 0.3f);
	// 	glBegin(GL_QUADS);
	// 	glVertex2i(rectangle_.left(),  rectangle_.top());
	// 	glVertex2i(rectangle_.right(), rectangle_.top());
	// 	glVertex2i(rectangle_.right(), rectangle_.bottom());
	// 	glVertex2i(rectangle_.left(),  rectangle_.bottom());
	// 	glEnd();
	// 	glDisable(GL_BLEND);
	// 	glDepthMask(GL_TRUE);


	// draw a line

	// 	glLineWidth(2.0);
	// 	glColor4f(0.0f, 1.0f, 1.0f, 0.5f);
	//  	glBegin(GL_LINES);
	//  	glVertex2i(rectangle_.left(),  rectangle_.top());
	//     glVertex2i(rectangle_.right(), rectangle_.bottom());
	//     glEnd();

	glEnable(GL_LIGHTING);


	canvas_->stopScreenCoordinatesSystem();
}
void PaintStroke::draw()
{
	draw_rectangle();
	draw_stroke();

	//draw hightlight vertex
	Sample& sample = SampleSet::get_instance()[cur_sample_to_operate_];
	LOCK(sample);
	Matrix44 adjust_mat = sample.matrix_to_scene_coord();

	glPointSize(Paint_Param::g_point_size);
	glBegin(GL_POINTS);
	for ( IndexType v_idx = 0; v_idx < sample.num_vertices(); v_idx++ )
	{
		ColorType c;
		if (sample[v_idx].is_selected() == true )
		{
			c = SELECTED_COLOR;
		}
		else
			c = HIGHTLIGHTED_COLOR;
		glColor4f(  c(0), c(1), c(2),c(3) );
		sample[v_idx].draw_without_color(adjust_mat);
	}
	glEnd();
	UNLOCK(sample);
}
// void PaintStroke::paintEvent(QPaintEvent* e)
//{
// 	Q_UNUSED(e);
// 	QPainter painter;
// 	painter.begin(this);
// 	painter.setRenderHint(QPainter::Antialiasing);

//}