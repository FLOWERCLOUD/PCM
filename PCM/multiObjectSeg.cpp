#include "multiObjectSeg.h"

void MultiObjectSeg::run()
{
	Logger<<"Begin multi-object segment!\n";

	std::vector<IndexType> changePs;

	detectChangePoints(changePs);

	Logger<<"End multi-object segment!\n";
}

void MultiObjectSeg::detectChangePoints(std::vector<IndexType>& changePoints)
{

	SampleSet& set = SampleSet::get_instance();

	if (set.empty())
	{
		Logger<< " End Clustering(Warning: Empty Set).\n";
		emit finish_compute();
		return;
	}

	PCLXYZCloud srCloud (new pcl::PointCloud<pcl::PointXYZ>); 

	PCLXYZCloud tgCloud (new pcl::PointCloud<pcl::PointXYZ>); 

	convertCloudFormat(srCloud, set[0]);

	convertCloudFormat(tgCloud, set[5]);

	ScalarType resolution = 0.05f;

	pcl::octree::OctreePointCloudChangeDetector<pcl::PointXYZ> cgOctree (resolution);
	
	Logger<<"Tree depth = "<<cgOctree.getTreeDepth()<<endl;

	cgOctree.setInputCloud(srCloud);

	cgOctree.addPointsFromInputCloud();

	cgOctree.switchBuffers();

	cgOctree.setInputCloud(tgCloud);

	cgOctree.addPointsFromInputCloud();

	Logger<<"Tree depth = "<<cgOctree.getTreeDepth()<<endl;

	cgOctree.getPointIndicesFromNewVoxels(changePoints);

	// visualization of the change points
	showChangePoints(changePoints,set[5]);

}

void MultiObjectSeg::convertCloudFormat(PCLXYZCloud& tgCloud, Sample& srCloud)
{
	IndexType psSize = srCloud.num_vertices();

	tgCloud->width = psSize;

	tgCloud->height = 1;

	tgCloud->resize(psSize);

	for (unsigned int pIt = 0; pIt < psSize; ++ pIt )
	{
		Vertex& vtx = srCloud[pIt];

		tgCloud->points[pIt].x = vtx.x();
		tgCloud->points[pIt].y = vtx.y();
	    tgCloud->points[pIt].z = vtx.z();	
	}

}

void MultiObjectSeg::showChangePoints(vector<IndexType>& changePs, Sample& pCLoud)
{
	vector<bool> isChange(pCLoud.num_vertices(),false);

	for (auto cIter = changePs.begin(); cIter != changePs.end(); ++ cIter)
	{
		isChange[(*cIter)] = true;
	}

	IndexType i = 0;
	for (Sample::vtx_iterator vIter = pCLoud.begin();
		  vIter != pCLoud.end(); ++ vIter, ++i)
	{
		if (isChange[i])
		{
			(*vIter)->set_visble(false);
			//(*vIter)->set_label();
		}
	}
}

