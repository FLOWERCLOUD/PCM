#ifndef _SCANNER_H
#define _SCANNER_H
#include <OpenNI.h>
#include <QThread>
#include "basic_types.h"

class Sample;

class Scanner: public QThread{
		Q_OBJECT
public:


	void run() Q_DECL_OVERRIDE;
	inline Sample* generate_point_cloud( const IndexType sample_idx );
	void build_sample_normal( Sample& smp );

private:
	openni::Device		device_;
	openni::VideoStream  stream_depth_;
	openni::VideoStream  stream_color_;
	openni::VideoFrameRef cur_frame_depth_;
	openni::VideoFrameRef cur_frame_color_;

	static IndexType	min_scan_dist_;
	static IndexType max_scan_dist_;
	static IndexType scan_frame_step_size_;
	static IndexType min_frame_;
	static IndexType max_frame_;

public:
signals:
	void finished_scan();
};

#endif	_SCANNER_H