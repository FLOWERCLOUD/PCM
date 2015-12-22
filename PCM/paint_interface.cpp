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
	canvas_->updateGL();

}

void PaintStroke::release(QMouseEvent *e)
{
	if (left_mouse_button_ == true)
	{
		// Possibly swap left/right and top/bottom to make rectangle_ valid.
		rectangle_ = rectangle_.normalized();
		rectangle_ = QRect(e->pos(), e->pos());
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
		canvas_->updateGL();
	}

}

void PaintStroke::draw_stroke()
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

	glLineWidth(2.0);
	glColor4f(0.0f, 1.0f, 1.0f, 0.5f);
 	glBegin(GL_LINES);
 	glVertex2i(rectangle_.left(),  rectangle_.top());
    glVertex2i(rectangle_.right(), rectangle_.bottom());
    glEnd();

	glEnable(GL_LIGHTING);


	canvas_->stopScreenCoordinatesSystem();
}

void PaintStroke::draw()
{
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