#include "scanner.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "basic_types.h"
#include "sample.h"
#include "sample_set.h"

using namespace openni;


IndexType Scanner::min_scan_dist_ = 10;
IndexType Scanner::max_scan_dist_ = 800;
IndexType Scanner::min_frame_ = 100;
IndexType Scanner::max_frame_ = 1000;
IndexType Scanner:: scan_frame_step_size_ = 10;


void Scanner::run()
{
	

	//Initialize OpenNI environment 
	OpenNI::initialize();
	SampleSet::get_instance().clear();
	device_.open( ANY_DEVICE );

	stream_color_.create( device_, SENSOR_COLOR );
	stream_depth_.create( device_, SENSOR_DEPTH );

	VideoMode	mode_depth;
	mode_depth.setResolution( 640, 480 );
	mode_depth.setFps( 30 );
	mode_depth.setPixelFormat( PIXEL_FORMAT_DEPTH_1_MM );
	stream_depth_.setVideoMode( mode_depth );

	VideoMode mode_color;
	mode_color.setResolution( 640, 480 );
	mode_color.setFps( 30 );
	mode_color.setPixelFormat( PIXEL_FORMAT_RGB888 );
	stream_color_.setVideoMode( mode_color );

	//if( device_.isImageRegistrationModeSupported(IMAGE_REGISTRATION_DEPTH_TO_COLOR) )
	//{
	//	device_.setImageRegistrationMode( IMAGE_REGISTRATION_DEPTH_TO_COLOR );
	//}

	//Open data stream
	stream_color_.start();
	stream_depth_.start();


	//cv::namedWindow( "depth image", CV_WINDOW_AUTOSIZE );
	cv::namedWindow( "color image", CV_WINDOW_AUTOSIZE );

	IndexType max_depth = stream_depth_.getMaxPixelValue();
	IndexType cur_frame_idx = 0;
	IndexType ith_step = 0;

	SampleSet::get_instance().clear();

	cv::VideoWriter video_writer;
	video_writer.open("..\\sample_video\\video.avi", 0, 30 , cv::Size(640,480), true);

	Logger<<"Begin to scan point cloud, Press q at scan window to terminate ......"<<std::endl;

	while (true)
	{
		stream_depth_.readFrame( &cur_frame_depth_ );
		stream_color_.readFrame( &cur_frame_color_ );


		//depth data convert to OpenCV format
		const cv::Mat image_depth(  cur_frame_depth_.getHeight(), cur_frame_depth_.getWidth(),CV_16UC1,
													( void* )cur_frame_depth_.getData()  ); 
		//make sure depth image display apparently, convert CV_16U1 to CV_8U
		cv::Mat scaled_depth;
		image_depth.convertTo( scaled_depth, CV_8U, 255.0/max_depth );
		//cv::imshow( "depth image", scaled_depth );

		//color data convert to OpenCV format
		const cv::Mat image_color_RGB(  cur_frame_color_.getHeight(), cur_frame_color_.getWidth(),CV_8UC3,
			( void* )cur_frame_color_.getData()  ); 


		//convert RGB format to BGR format
		cv::Mat image_color_BGR;
		cv::cvtColor( image_color_RGB, image_color_BGR, CV_RGB2BGR );
		cv::imshow( "color image", image_color_BGR );

		//save to video file
		video_writer << image_color_BGR;

		//Generate Point Cloud

		if ( cur_frame_idx<max_frame_ && cur_frame_idx==min_frame_+ith_step*scan_frame_step_size_ )
		{
			Logger << ith_step<<"th scan\n";
			Sample* new_sample = generate_point_cloud( cur_frame_idx-min_frame_ );
			SampleSet::get_instance().push_back( new_sample );
			ith_step++;
		}

		cur_frame_idx++;

		if ( cv::waitKey(1)=='q' )
		{
			break;
		}

	}
	
	video_writer.release();

	Logger<<"End scan."<<std::endl;


	stream_depth_.destroy();
	stream_color_.destroy();

	device_.close();
	//shut down device
	OpenNI::shutdown();

	emit finished_scan();

}

Sample* Scanner::generate_point_cloud( const IndexType sample_idx)
{

	Sample*	new_sample = new Sample;

	IndexType		num_point = cur_frame_depth_.getWidth() * cur_frame_depth_.getHeight();

	const DepthPixel*		depth_data_raw = (const DepthPixel*)cur_frame_depth_.getData();
	const RGB888Pixel*  image_data_raw = (const RGB888Pixel*)cur_frame_color_.getData();

	IndexType x, y;

	for ( y=0; y < cur_frame_depth_.getHeight(); y+=1 )
	{
		IndexType idx_shift = y * cur_frame_depth_.getWidth();
		for (x=0; x < cur_frame_depth_.getWidth(); x+=1 )
		{
			IndexType idx = idx_shift + x;
			ScalarType fx, fy, fz;
			
			if ( depth_data_raw[idx] <min_scan_dist_ || depth_data_raw[idx]>max_scan_dist_ )continue;

			CoordinateConverter::convertDepthToWorld( stream_depth_, x, y, depth_data_raw[idx], &fx, &fy, &fz );

			IndexType pcx, pcy;
			CoordinateConverter::convertDepthToColor( stream_depth_, stream_color_, x,y,depth_data_raw[idx],&pcx,&pcy );

			IndexType color_idx = pcy * cur_frame_depth_.getWidth() + pcx;
			RGB888Pixel color = image_data_raw[ color_idx ];

			//auto r = image_data_raw[idx].r;
			//auto g = image_data_raw[idx].g;
			//auto b = image_data_raw[idx].b;
			//ColorType cv(r/255., g/255., b/255., 1.);

			PointType v(fx, fy, fz);
			ColorType cv(color.r/255.,color.g/255.,color.b/255.,1.);
			
			NormalType nv(0.,0.,-1.);

			new_sample->add_vertex( v, nv, cv );
			
		}
	}


	new_sample->set_visble(false);
	new_sample->set_color( Color_Utility::span_color_from_table( sample_idx ) );
	new_sample->build_kdtree();


	return new_sample;


}

