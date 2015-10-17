#ifndef _MULTIWAY_PROPAGATE_H
#define _MULTIWAY_PROPAGATE_H

#include "pool_allocator.h"
#include "tracer.h"

// #include "boost/graph/adjacency_list.hpp"
// #include "boost/graph/named_graph.hpp"


#include <queue>
#include <vector>

#include "GraphMatching.h"

#include <gco-v3.0/GCoptimization.h>

#include "DeformaleRegistration.h"

#define  COUT_DEBUG 0

using namespace std;


class DualwayPropagation
{
public:

	DualwayPropagation()
	{
		downSample = NULL;
	}
	~DualwayPropagation()
	{
		if( downSample){
			delete downSample;
		}
	}

	void propagate_back_with_link(IndexType f);
	
	void propagate_back(IndexType f)
	{
		if ( components_.find(f+1)==components_.end() )
		{
			return;
		}
		CFrame& cur_frame = components_[f];
		CFrame& next_frame = components_[f+1];

		set<IndexType> new_label_set;

		for ( IndexType l=0; l<cur_frame.label_bucket.size(); ++l )
		{
			CLabel& cur_label = *(cur_frame.label_bucket[l]);
			map<IndexType, set<IndexType> > mapper;
			for (auto v = cur_label.vertex_bucket.begin();
					v != cur_label.vertex_bucket.end();
					v++)
			{
				CVertex& cur_vtx = *(v->second);
				CVertex& corr_vtx = *cur_vtx.next_corr;
				IndexType corr_label = corr_vtx.label_parent->label_id;
				auto find_if_exist = mapper.find( corr_label );
				if ( find_if_exist==mapper.end() )
				{
					mapper.insert( make_pair(corr_label, set<IndexType>()) );
				}
				mapper[corr_label].insert( corr_vtx.vtx_id );
			}

			//check whether split two part
			bool to_spilt = false;
			map<IndexType, set<IndexType> >::iterator component_to_spilt;
			if (mapper.size()==1)
			{
				component_to_spilt = mapper.begin();
				to_spilt = true;
			}
			else
			{
				bool first_encounter_no_new_label = true;
				for ( auto miter = mapper.begin(); miter != mapper.end(); miter++ )
				{
					IndexType l = miter->first;
					if( new_label_set.find(l)!=new_label_set.end() )
						continue;
					if ( first_encounter_no_new_label )
					{
						first_encounter_no_new_label = !first_encounter_no_new_label;
						component_to_spilt = miter;
						to_spilt = true;
					}
					else
					{
						// contain two old components
						to_spilt = false;
						break;
					}
				}
			}


			if( to_spilt )
			{
				 auto only = component_to_spilt;
				 IndexType map_label = only->first;
				 set<IndexType>& map_bucket = only->second;
				 map<IndexType, CVertex*>&  real_bucket = next_frame.label_bucket[map_label]->vertex_bucket;
				 if( map_bucket.size() *2 < real_bucket.size() )
				 {
					 //now split

					 //now create new component
					 IndexType new_label = next_frame.label_bucket.size();
					 new_label_set.insert(new_label);
					 next_frame.label_bucket.push_back( (CLabel*)0 );
					 CLabel* new_label_space = allocator_.allocate<CLabel>();
					 CLabel* new_label_obj = new (new_label_space)CLabel;
					 next_frame.label_bucket[new_label]	= new_label_obj;			
					 next_frame.label_bucket[new_label]->label_id = new_label;
					 next_frame.label_bucket[new_label]->frame_parent = &next_frame;

					 //delete elem from real_bucket that appear in map_bucket
					 auto to_remove = map_bucket.begin();
					 for (auto vv = real_bucket.begin(); vv!=real_bucket.end() && to_remove!=map_bucket.end(); )
					 {
						 IndexType vtx_idx = vv->first;
						 if( *to_remove == vtx_idx )
						 {
							 next_frame.label_bucket[new_label]->vertex_bucket.insert( *vv );
							 vv = real_bucket.erase( vv );
							 next_frame.label_of_vtx[ vtx_idx ] = new_label;
						 }
						 else if( *to_remove > vv->first )
							vv++;
						 else
							 to_remove++;
					 }

					 for (auto vv = next_frame.label_bucket[new_label]->vertex_bucket.begin();
						 vv !=next_frame.label_bucket[new_label]->vertex_bucket.end(); vv++ )
					 {
						 vv->second->label_parent = next_frame.label_bucket[new_label];
					 }
				 }


			}
		}
		smooth_label(f+1);
	}

	void propagate_front(IndexType f)
	{
		if ( components_.find(f-1)==components_.end() )
		{
			return;
		}
		CFrame& cur_frame = components_[f];
		CFrame& prev_frame = components_[f-1];

		set<IndexType> new_label_set;
		for ( IndexType l=0; l<cur_frame.label_bucket.size(); ++l )
		{
			CLabel& cur_label = *(cur_frame.label_bucket[l]);
			map<IndexType, set<IndexType> > mapper;
			for (auto v = cur_label.vertex_bucket.begin();
				v != cur_label.vertex_bucket.end();
				v++)
			{
				CVertex& cur_vtx = *(v->second);
				CVertex& corr_vtx = *cur_vtx.prev_corr;
				IndexType corr_label = corr_vtx.label_parent->label_id;
				auto find_if_exist = mapper.find( corr_label );
				if ( find_if_exist==mapper.end() )
				{
					mapper.insert( make_pair(corr_label, set<IndexType>()) );
				}
				mapper[corr_label].insert( corr_vtx.vtx_id );
			}

			//check whether split two part
			bool to_spilt = false;
			map<IndexType, set<IndexType> >::iterator component_to_spilt;
			if (mapper.size()==1)
			{
				component_to_spilt = mapper.begin();
				to_spilt = true;
			}
			else
			{
				bool first_encounter_no_new_label = true;
				for ( auto miter = mapper.begin(); miter != mapper.end(); miter++ )
				{
					IndexType l = miter->first;
					if( new_label_set.find(l)!=new_label_set.end() )
						continue;
					if ( first_encounter_no_new_label )
					{
						first_encounter_no_new_label = !first_encounter_no_new_label;
						component_to_spilt = miter;
						to_spilt = true;
					}
					else
					{
						// contain two old components
						to_spilt = false;
						break;
					}
				}
			}


			if( to_spilt )
			{
				auto only = component_to_spilt;
				IndexType map_label = only->first;
				set<IndexType>& map_bucket = only->second;
				map<IndexType, CVertex*>&  real_bucket = prev_frame.label_bucket[map_label]->vertex_bucket;
				if( map_bucket.size() *2 < real_bucket.size() )
				{
					//now split

					//now create new component
					IndexType new_label = prev_frame.label_bucket.size();
					new_label_set.insert(new_label);
					prev_frame.label_bucket.push_back( (CLabel*)0 );
					CLabel* new_label_space = allocator_.allocate<CLabel>();
					CLabel* new_label_obj = new (new_label_space)CLabel;
					prev_frame.label_bucket[new_label]	= new_label_obj;			
					prev_frame.label_bucket[new_label]->label_id = new_label;
					prev_frame.label_bucket[new_label]->frame_parent = &prev_frame;

					//delete elem from real_bucket that appear in map_bucket
					auto to_remove = map_bucket.begin();
					for (auto vv = real_bucket.begin(); vv!=real_bucket.end() && to_remove!=map_bucket.end(); )
					{
						IndexType vtx_idx = vv->first;
						if( *to_remove == vtx_idx )
						{
							prev_frame.label_bucket[new_label]->vertex_bucket.insert( *vv );
							vv = real_bucket.erase( vv );
							prev_frame.label_of_vtx[ vtx_idx ] = new_label;
						}
						else if( *to_remove > vv->first )
							vv++;
						else
							to_remove++;
					}
					
					for (auto vv = prev_frame.label_bucket[new_label]->vertex_bucket.begin();
						vv !=prev_frame.label_bucket[new_label]->vertex_bucket.end(); vv++ )
					{
						vv->second->label_parent = prev_frame.label_bucket[new_label];
					}
				}


			}
		}
		smooth_label(f-1);
	}

	void read_data(char *label_name,char *corr_name);


	void init_labeles_graph();

	void getEdgeVertexs( IndexType _CFrameId ,IndexType lLabelId ,IndexType rLabelId , map<IndexType, map<IndexType ,HVertex*> >& _edgepoints );
	void getEdgeVertexs2( IndexType _CFrameId , distanPriQueue& _PriQuemap, map<IndexType, map<IndexType ,HVertex*> >& _edgepoints );
	void init_labeles_graph_hier(ScalarType distThr);


	void read_label_file(char *filename)
	{
		FILE* in_file = fopen(filename, "r");
		if (in_file==NULL)
		{
			return;
		}
		IndexType frame, label, vtx_idx;
		while ( true )
		{
			int stat =  fscanf(in_file, "%d %d %d\n",&frame, &label, &vtx_idx);
			if (stat == EOF)
				break;
			if ( components_.find(frame)==components_.end() )
			{
				components_.insert(make_pair(frame, CFrame()));
				components_[frame].frame_id = frame;
			}
			if ( label >= components_[frame].label_bucket.size() )
			{
				components_[frame].label_bucket.resize( label+1 );
			}
			if (   nullptr==components_[frame].label_bucket[label] )
			{
				CLabel* new_label_space = allocator_.allocate<CLabel>();
				CLabel* new_label = new (new_label_space)CLabel;
				components_[frame].label_bucket[label] = new_label;
				components_[frame].label_bucket[label]->frame_parent = &components_[frame];
				components_[frame].label_bucket[label]->label_id = label;
			}
			CVertex* new_space = allocator_.allocate<CVertex>();
			CVertex* new_vtx = new (new_space)CVertex(vtx_idx, components_[frame].label_bucket[label]);
			components_[frame].label_bucket[label]->vertex_bucket.insert( make_pair(vtx_idx,new_vtx) );
			components_[frame].label_of_vtx.insert( make_pair(vtx_idx, label) );
			//ADD BY HUAYUN 
			SampleSet& sample_set = SampleSet::get_instance();
			Sample& sample = sample_set[frame];
			 
		 

		 	 

		}
	}

	void read_corres_file(char *filename)
	{
		FILE* in_file = fopen(filename,"r");
		if(in_file==NULL)
		{
			return;
		}
		IndexType cur_frame, cur_vtx_idx, next_frame, next_vtx_idx;
		while (true)
		{
			int stat = fscanf(in_file,"%d %d %d %d\n",&cur_frame, &cur_vtx_idx, &next_frame, &next_vtx_idx);
			if(stat==EOF)break;

			if( hier_componets_.find(cur_frame)==hier_componets_.end() || hier_componets_.find(next_frame)==hier_componets_.end())
				continue;

			IndexType label = hier_componets_[cur_frame].label_of_vtx[cur_vtx_idx];

			IndexType next_label = hier_componets_[next_frame].label_of_vtx[next_vtx_idx];

			HVertex& cur_vtx = *hier_componets_[cur_frame].hier_label_bucket[0][label]->vertex_bucket[cur_vtx_idx];

			if ( cur_frame+1 == next_frame  )
			{
				cur_vtx.next_corr = hier_componets_[next_frame].hier_label_bucket[0][ next_label]->vertex_bucket[next_vtx_idx];
			}
			else if (cur_frame-1 == next_frame)
			{
				cur_vtx.prev_corr = hier_componets_[next_frame].hier_label_bucket[0][ next_label]->vertex_bucket[next_vtx_idx];
			}
		}
	}

	void  read_label_file_hier(char *filename);



	void smooth_label(IndexType frame_idx)
	{
		Sample& orig_smp = SampleSet::get_instance()[frame_idx];
		Sample* downsmp_ptr = new Sample;
		map<IndexType, IndexType> idx_mapper;
		IndexType i=0;
		for ( auto viter = components_[frame_idx].label_of_vtx.begin();
				viter != components_[frame_idx].label_of_vtx.end();
				viter++,i++)
		{
			IndexType vtx_idx = viter->first;
			Vertex& vtx = orig_smp[vtx_idx];
			
			PointType v( vtx.x(), vtx.y(), vtx.z() );
			ColorType cv(vtx.r(), vtx.g(), vtx.b(), vtx.alpha());
			NormalType nv(vtx.nx(), vtx.ny(), vtx.nz());

			downsmp_ptr->add_vertex(v,nv,cv);

			idx_mapper.insert(make_pair(i, vtx_idx));
		}

		downsmp_ptr->build_kdtree();

		const IndexType k =60;
		IndexType neighbours[k];
		ScalarType dist[k];
		i = 0;
		map<IndexType, IndexType> old_label_map = components_[frame_idx].label_of_vtx;
		for ( i=0; i<idx_mapper.size(); i++)
		{
			vector<int> label_count( components_[frame_idx].label_bucket.size(), 0 );
			IndexType real_vtx_idx = idx_mapper[i];
			downsmp_ptr->neighbours(i, k, neighbours, dist );
			for ( IndexType neig_id = 0; neig_id<k; neig_id++ )
			{
				IndexType real_neig_id = idx_mapper[neighbours[neig_id]];
				IndexType neig_label = old_label_map[real_neig_id];
				label_count[neig_label]++;
			}
			IndexType max_freq_label = max_element(label_count.begin(), label_count.end()) - label_count.begin();
			IndexType cur_vtx_label = old_label_map[real_vtx_idx];

			if ( max_freq_label!=cur_vtx_label )
			{
				//change label
				auto vv = components_[frame_idx].label_bucket[cur_vtx_label]->vertex_bucket.find(real_vtx_idx);
				assert(vv!=components_[frame_idx].label_bucket[cur_vtx_label]->vertex_bucket.end());
				components_[frame_idx].label_bucket[max_freq_label]->vertex_bucket.insert(*vv);
				components_[frame_idx].label_bucket[cur_vtx_label]->vertex_bucket.erase(vv);
				components_[frame_idx].label_of_vtx[real_vtx_idx] = max_freq_label;
			}
		}

		delete downsmp_ptr;

	}

	void write_label_file(std::string filename)
	{

		FILE* outfile = fopen(filename.c_str(),"w");
		for ( auto frame_iter = components_.begin();
				frame_iter != components_.end();
				frame_iter++ )
		{
			CFrame& frame = frame_iter->second;
			IndexType frame_idx = frame_iter->first;
			for ( auto label_iter = frame.label_bucket.begin(); 
				label_iter!=frame.label_bucket.end(); label_iter++  )
			{
				CLabel& label = **label_iter;
				IndexType label_idx = label.label_id;
				for ( auto vtx_iter = label.vertex_bucket.begin();
					vtx_iter!=label.vertex_bucket.end();
					vtx_iter++)
				{
					CVertex& vtx = *(vtx_iter->second);
					fprintf(outfile, "%d %d %d\n", frame_idx, label_idx, vtx.vtx_id);
				}
			}
		}

		fclose(outfile);
	}

	void wirteSplitGraphLables(std::string filename);

	void compute()
	{
		for (auto citer = components_.begin(); citer!=components_.end(); citer++)
		{
			propagate_back_with_link(citer->first);
			//propagate_back(citer->first);
			
		}

		//show_corresponding(14);


		for (auto citer = components_.rbegin(); citer!=components_.rend(); citer++)
		{
			//propagate_front(citer->first);
		}


 		for (auto citer = components_.begin(); citer!=components_.end(); citer++)
 		{
 			//propagate_back(citer->first);
 		}

		//for (auto citer = components_.rbegin(); citer!=components_.rend(); citer++)
		//{
		//	propagate_front(citer->first);
		//}

		return;
	}

	void buildSmapleKDtree(CFrame* smpCloud, Sample* smp);

	void point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec);

	void buildKdTree(IndexType _cframeid);

	//0821
	void split_twoAjacent_graph_next(IndexType srFrame, IndexType tgFrame);

	void split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame);

	void split_nest_graph_prev(IndexType startFrame,IndexType srFrame, IndexType tgFrame);

	void splitAllSquenceGraph(IndexType iterN);

	void getNextCorVtx(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor);

	IndexType checkNextLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor);

	bool checkPrevLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor, IndexType& isSplit);

	void minDistBeTwoParts(IndexType cFrame,map<IndexType,HVertex*>& fPart,map<IndexType,HVertex*>& sPart, ScalarType& misDis);

	void determinateUnmarkPoints(IndexType cFrame,map<IndexType,HVertex*>& unMarkPs,vector<HLabel*> oriLabelBucket,IndexType nodeId,IndexType newLabe,IndexType tgSize);

	ScalarType p2PatchAvgDis(IndexType cFrame, PointType& pCoor,map<IndexType,HVertex*>& parthPs);//随机点的距离

	ScalarType p2PatchMinDis(IndexType cFrame, PointType& pCoor,map<IndexType,HVertex*>& parthPs);//所有点的最小距离

	void copyLabelBucket(vector<HLabel*>& leftLabels, const vector<HLabel*>& oriLabels);

	void checkPsNewLabelParentPtr(vector<HLabel*> oriLabelBucket,IndexType labParSize);

	void smoothAfterSplit();

	void smoothSingleFrame(IndexType frameId);

	void buildMainPatchMatching();

	void buildNextPatchCorr(IndexType srFrame,IndexType tgFrame);
	
	void buildPrevPatchCorr(IndexType srFrame,IndexType tgFrame);

	void buildSquenceUniqLabel(std::string filename);

	void writeSquenceLabel(map<IndexType,IndexType>& patchTrajId,std::string filename);

	map<IndexType,HFrame>* getCompents(){return &hier_componets_;}

	void read_label_file_coseg(char *label_file_name);

	void read_corr_file_coseg(char* corr_file_name);

	void buildPatchCorrespondenceByLabel();

	void mergePatchesAfterCoSeg();

	void findBestPatches(IndexType srFrame, IndexType tgFrame, set<IndexType>& srBestSet, set<IndexType>& tgBestSet); 

	void mergeSingleTinyPatches(IndexType vSize);

	//split by ordered edges 0903
	void split_twoAjacent_graph_next_order(IndexType srFrame, IndexType tgFrame);

	void split_twoAjacent_graph_prev_order(IndexType srFrame, IndexType tgFrame);

	void generateOrderededges(IndexType srFrame, IndexType tgFrame);

	void constructPCloudGraph();

	void addGraphVertex(PCloudGraph& pcGraph, IndexType frameId);

	void addGraphEdge(PCloudGraph& pcGraph, IndexType frameId);

	ScalarType p2PatchGeoDis(IndexType cFrame, HVertex& oriP, map<IndexType,HVertex*>& parthPs);

	//
	void generateOrderPrevEdges(IndexType srFrame, IndexType tgFrame);

	void show_corresponding(IndexType f);

	void mergePatchTraj();

	void generTrajNodes(vector<PatchTraj>& pNodes);// 生成轨迹节点

	void graphCuts(vector<PatchTraj>& pNodes, vector<IndexType>& labels);

	void mergeSquences(vector<IndexType>& labes);

	void setSegNeihbor(vector<PatchTraj>& pNodes, GCoptimizationGeneralGraph& segGraphC);

	void setSegDataItem(GCoptimizationGeneralGraph& segGraphC);

	void setSegSmoothItem(GCoptimizationGeneralGraph& segGraphC);

	bool isAdjInSeq(PatchTraj& nodeA, PatchTraj& nodeB);

	void getSegLabels(GCoptimizationGeneralGraph& segGraphC, vector<IndexType>& labels);

	ScalarType motionSimilarityBetw2Nodes(IndexType i, IndexType j, vector<PatchTraj>& oriTraj);// i connect j

	void calculateTrans(IndexType lab,IndexType sFrame, IndexType tFrame, Matrix34& toTrans, Matrix34& backTrans);

	void getCoorByVtxBucket(IndexType lab, IndexType frameId, Matrix3X& vtxCoor);


public:

	Sample* downSample;

private:
	map<IndexType, CFrame> components_;
	map<IndexType,IndexType> frame_index_map_label_; 
	static  PoolAllocator allocator_;

	map<IndexType, HFrame> hier_componets_;

	map<IndexType,IndexType> o2d_map; // origin to down
	map<IndexType,IndexType> d2o_map;

public:
		OrderEdgeQueue orderedEdgeQ;
};


#endif