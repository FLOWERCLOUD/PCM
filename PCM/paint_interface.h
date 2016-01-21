#ifndef  _PAINT_INTERFACE_
#define  _PAINT_INTERFACE_

#include "tool.h"
#include "basic_types.h"
#include <vector>

class PaintCanvas;

/* Draw sketch tool*/
class PaintStroke : public Tool
{
public:
	PaintStroke(PaintCanvas* canvas): Tool(canvas){}
	~PaintStroke(){}

	virtual void move(QMouseEvent *e);
	virtual void drag(QMouseEvent *e);
	virtual void release(QMouseEvent *e);
	virtual void press(QMouseEvent* e);
	virtual void draw();
	//virtual void paintEvent(QPaintEvent* e);

private:
	void draw_stroke();
	void draw_rectangle();

private:
	QPoint mouse_pressed_pos;
	QPoint mouse_move_pos;
	QRect  rectangle_;

private:
	std::vector<QPoint> crtSketch;

};


#endif