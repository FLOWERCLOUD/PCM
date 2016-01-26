#ifndef _MULTI_OBJECT_SEG_H
#define _MULTI_OBJECT_SEG_H

#include <QThread>
#include "sample_set.h"
#include <vector>
#include "globals.h"
#include "basic_types.h"
#include "sample.h"
#include <algorithm>

#include"pcl/octree/octree.h"
#include "pcl/octree/octree_pointcloud_changedetector.h"
#include "pcl/point_cloud.h"

#include "pcl/octree/impl/octree2buf_base.hpp"
#include "pcl/octree/impl/octree_base.hpp"
#include"pcl/octree/impl/octree_pointcloud.hpp"

using namespace std;

typedef pcl::PointCloud<pcl::PointXYZ>::Ptr PCLXYZCloud;

class MultiObjectSeg : public QThread
{
	Q_OBJECT

public:
	MultiObjectSeg(){};
	~MultiObjectSeg(){};

public:
	void run() Q_DECL_OVERRIDE;

signals:
	void finish_compute();

public:
	void detectChangePoints(std::vector<IndexType>& changePoints);

	void convertCloudFormat(PCLXYZCloud& tgCloud, Sample& srCloud);

	void showChangePoints(vector<IndexType>& changePs, Sample& pCloud);
};

#endif