#ifndef _MULTIWAY_PROPAGATE__H
#define _MULTIWAY_PROPAGATE__H

#include "co_segmentation.h"
#include "pool_allocator.h"
#include "tracer.h"

#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/named_graph.hpp"

#define  INF 1000000



#define frame_index_to_key(f,i) ((f<<16)|i)
#define frame_label_to_key(f,l) ((f<<16)|l)
#define get_index_from_key(k) (k&0xffff)
#define get_frame_from_key(k) (k>>16)

struct CVertex;
struct CLabel;
struct CFrame; 


struct HVertex;
struct HLabel;
struct HFrame; 

struct GraphVertexProperty;
struct GraphEdgeProperty;

typedef boost::adjacency_list<boost::listS, boost::vecS,boost::undirectedS,GraphVertexProperty,GraphEdgeProperty> LabelsGraph;
// 节点描述符
typedef boost::graph_traits<LabelsGraph>::vertex_descriptor   VertexDescriptor; 
// 边描述符
typedef boost::graph_traits<LabelsGraph>::edge_descriptor     EdgeDescriptor;
// 下面是邻接链表的一些遍历器
typedef boost::graph_traits<LabelsGraph>::vertex_iterator   VertexIterator;
typedef boost::graph_traits<LabelsGraph>::edge_iterator   EdgeIterator;
typedef boost::graph_traits<LabelsGraph>::adjacency_iterator   AdjacencyIterator;
typedef boost::graph_traits<LabelsGraph>::out_edge_iterator   OutEdgeIterator;



struct CVertex
{

	CVertex(IndexType vi, CLabel* lab):vtx_id(vi), label_parent(lab)
		,prev_corr(NULL),next_corr(NULL){}

	IndexType vtx_id;
	CLabel* label_parent;
	CVertex* prev_corr;
	CVertex* next_corr;
};

struct CLabel
{
	IndexType label_id;
	CFrame* frame_parent;
	map<IndexType,CVertex*> vertex_bucket;

	CLabel* prev_corr;
	CLabel* next_corr; //path correspondence

};

struct CFrame
{
	IndexType frame_id;
	vector<CLabel*> label_bucket;
	map<IndexType, IndexType> label_of_vtx;

	LabelsGraph* labelGraph; //保存块之间的图结构,一块对应一个节点

	vector<LabelsGraph*> hier_graph;
	//CFrame(){}
	//~CFrame(){}
	//CFrame(const CFrame& frame):frame_id(frame.frame_id),label_of_vtx(frame.label_of_vtx),labelGraph(new LabelsGraph(*frame.labelGraph) )
	//{
	//	if (! frame.label_bucket.empty())
	//	{
	//		
	//		for (auto iter = frame.label_bucket.begin(); iter != frame.label_bucket.end(); iter ++)
	//		{
	//			
	//			CLabel * temp = new CLabel;
	//			memcpy(temp,*iter,sizeof(CLabel) );
	//			label_bucket.push_back(temp);
	//		}
	//	}
	//}
};

struct HVertex
{

	HVertex(IndexType vi, HLabel* lab):vtx_id(vi), label_parent(lab)
		,prev_corr(NULL),next_corr(NULL){}

	IndexType vtx_id;
	HLabel* label_parent;
	HVertex* prev_corr;
	HVertex* next_corr;
};

struct HLabel
{
	IndexType label_id;
	HFrame* frame_parent;
	map<IndexType,HVertex*> vertex_bucket;

	HLabel* prev_corr;
	HLabel* next_corr; //path correspondence

};


struct HFrame
{
	IndexType frame_id;
	vector< vector<HLabel*> > hier_label_bucket; //帧的层次块结构
	map<IndexType,IndexType> label_of_vtx;
	vector<LabelsGraph*> hier_graph;//帧的层次图
};

struct GraphVertexProperty
{
	GraphVertexProperty(){}
	GraphVertexProperty(IndexType id, IndexType _prev, IndexType _next):
		index(id),prev(_prev),next(_next){}
		                
	IndexType index;
	IndexType prev;
	IndexType next;
};

struct GraphEdgeProperty
{
    IndexType index;
	IndexType start_;
	IndexType end_;

	map<IndexType, map<IndexType ,HVertex*> > edgePoints;

	GraphEdgeProperty(){}

	bool operator<(const GraphEdgeProperty& first) const //后面必须要const,set插入原始必须要有<运算符.
	{
		return this->index < first.index;
	}

};




class DualwayPropagation
{
public:
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

			//check whether split two parts
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

	void read_data(char *label_name,char *corr_name)
	{
		read_label_file(label_name);

		//read_label_file_hier(label_name);

		read_corres_file(corr_name);

		// test patches distances

// 		set<IndexType> oriPatech;
// 		set<IndexType> tarPatch;
// 		oriPatech.insert(0);
// 		oriPatech.insert(2);
// 		tarPatch.insert(0);
// 
// 		//ScalarType error = distance_between_paches(1,oriPatech,2,tarPatch);
// 
// 		split_patches(1,oriPatech,2,tarPatch);
		
	}

	void init_labeles_graph()
	{
		IndexType nodeSize = 0;

		for (auto citer = components_.begin(); citer!=components_.end(); citer++)
		{
		    nodeSize = citer->second.label_bucket.size();

			LabelsGraph* new_labelGraph_space = allocator_.allocate<LabelsGraph>();

			LabelsGraph* new_labelGraph = new (new_labelGraph_space)LabelsGraph(nodeSize);


			//new_labelGraph->added_vertex(nodeSize);

			for (IndexType i = 0; i < nodeSize; i ++)
			{
				for (IndexType j = i + 1; j <= nodeSize; j++)
				{

					CLabel* fir = citer->second.label_bucket[i];

					CLabel* sec = citer->second.label_bucket[j];

					map<IndexType,CVertex*> fir_vtx = fir->vertex_bucket;

					map<IndexType,CVertex*> sec_vtx = sec->vertex_bucket;

					//计算两个块之间的最短距离
					ScalarType minDis = 0.0;

					if (minDis < 1.)
					{

			           boost::add_edge(i,j,*new_labelGraph);

					}

				}
			}

			citer->second.labelGraph = new_labelGraph;

		}




	
	}



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

			if( components_.find(cur_frame)==components_.end() || components_.find(next_frame)==components_.end())
				continue;

			IndexType label = components_[cur_frame].label_of_vtx[cur_vtx_idx];

			IndexType next_label = components_[next_frame].label_of_vtx[next_vtx_idx];
			CVertex& cur_vtx = *components_[cur_frame].label_bucket[label]->vertex_bucket[cur_vtx_idx];
			if ( cur_frame+1 == next_frame  )
			{
				cur_vtx.next_corr = components_[next_frame].label_bucket[ next_label]->vertex_bucket[next_vtx_idx];
			}
			else if (cur_frame-1 == next_frame)
			{
				cur_vtx.prev_corr = components_[next_frame].label_bucket[ next_label]->vertex_bucket[next_vtx_idx];
			}
		}
	}

	void  read_label_file_hier(char *filename);

	void show_corresponding(int f)
	{
		for ( IndexType l = 0; l<components_[f].label_bucket.size(); l++ )
		{
			CLabel& label = *components_[f].label_bucket[l];
			for ( auto viter = label.vertex_bucket.begin();
					viter!=label.vertex_bucket.end();
					viter++)
			{
				CVertex& vtx = *(viter->second);
				if ( vtx.next_corr )
				{
					Tracer::get_instance().add_record(f, vtx.vtx_id, f+1, vtx.next_corr->vtx_id);
				}
				if (vtx.prev_corr)
				{
					Tracer::get_instance().add_record(f, vtx.vtx_id, f-1, vtx.prev_corr->vtx_id);
				}
			}
		}
	}

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

	void compute()
	{
		for (auto citer = components_.begin(); citer!=components_.end(); citer++)
		{
			propagate_back(citer->first);
			
		}

		//show_corresponding(14);


		for (auto citer = components_.rbegin(); citer!=components_.rend(); citer++)
		{
			propagate_front(citer->first);
		}


 		for (auto citer = components_.begin(); citer!=components_.end(); citer++)
 		{
 			propagate_back(citer->first);
 		}

		//for (auto citer = components_.rbegin(); citer!=components_.rend(); citer++)
		//{
		//	propagate_front(citer->first);
		//}

		return;
	}

	ScalarType distance_between_paches(IndexType sFrame, set<IndexType>& sPatches,
		                               IndexType tFrame, set<IndexType>& tPatches)
	{
		assert(sPatches.size() > 0 && tPatches.size() > 0);

		map<IndexType,CVertex*> srVertex;

		set<IndexType>::iterator sIter = sPatches.begin();

		for (; sIter != sPatches.end(); sIter ++)
		{
			CLabel* temp = components_[sFrame].label_bucket[*sIter];
			srVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
		}

		map<IndexType,CVertex*> tgVertex;

		set<IndexType>::iterator tIter = tPatches.begin();

		for (; tIter != tPatches.end(); tIter++)
		{
			CLabel* temp = components_[tFrame].label_bucket[*tIter];
			tgVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
		}

		map<IndexType,CVertex*> toComVtx; // 记录那些点它们的对应会落在对应块上

		map<IndexType,CVertex*> backComVtx; // 相互记录

		bool isTo = true;
		computer_common(srVertex,tgVertex,toComVtx,isTo);

		isTo = false;
		computer_common(srVertex, tgVertex,backComVtx,isTo);


		// the distances should includes two items Jarcard and deformation distortion

		ScalarType srV_size = srVertex.size();
		ScalarType tgV_size = tgVertex.size();

		IndexType toCom_size = toComVtx.size();
		IndexType backCom_size = backComVtx.size();

		ScalarType corValue =  0.5 * (toCom_size/srV_size + backCom_size/tgV_size );

		ScalarType scaleShape =  abs( srV_size - tgV_size) / max(srV_size, tgV_size );
		
		//if ( scaleShape > 0.5 || corValue < 0.5)
		//{
		//	return INF;
		//}else
		//{
		//	return  1 - corValue;
		//}
		

		// calculate deformation distortion

		ScalarType totError = 0.0;
		const IndexType k = 300;
		IndexType neighbours[k];
		ScalarType dist[k];

		map<IndexType,CVertex*> valid_neig;

		map<IndexType,CVertex*>::iterator srComIter = toComVtx.begin();
			
		for (; srComIter != toComVtx.end(); srComIter ++)//---->
		{
			valid_neig.clear();
			Sample& curF = SampleSet::get_instance()[sFrame];
			curF.neighbours(srComIter->first, k, neighbours, dist);

			for (IndexType kk = 0; kk < k; kk ++)
			{
				 auto fiter = srVertex.find(neighbours[kk] );

				if ( fiter != srVertex.end() )
				{
					valid_neig.insert( *fiter );
				}
			}

			IndexType neig_siz = valid_neig.size();

			if(neig_siz<4)
			{
				Logger<<" neighbor size is small!.\n";
				continue; 
			}

			Matrix3X s_coord, t_coord;
			s_coord.setZero(3, neig_siz);
			t_coord.setZero(3, neig_siz);

			Sample& nextF = SampleSet::get_instance()[ tFrame ];

			IndexType vi = 0;
			IndexType kk = 0;

			for ( map<IndexType,CVertex*>::iterator neig_iter = valid_neig.begin();
				neig_iter != valid_neig.end(); ++neig_iter )
			{
				s_coord.col(kk) = curF.vertices_matrix().col(neig_iter->first );
				t_coord.col(kk) = nextF.vertices_matrix().col(neig_iter->second->next_corr->vtx_id);
			}
			
			Matrix33 rot_mat;
			MatrixXX tran_vec;


			point2point(s_coord, t_coord, rot_mat, tran_vec);

			auto bias = rot_mat * curF.vertices_matrix().col(srComIter->first) + tran_vec
				        - nextF.vertices_matrix().col(srComIter->second->next_corr->vtx_id);

			totError += bias.norm();

		}

		//<---- back deformation error

		map<IndexType,CVertex*>::iterator tgComIter = backComVtx.begin();

		for (; tgComIter != backComVtx.end(); tgComIter ++)
		{
			valid_neig.clear();

			Sample& curF = SampleSet::get_instance()[tFrame];

			curF.neighbours(tgComIter->first, k, neighbours, dist);

			for (IndexType kk = 0; kk < k; kk ++)
			{
				auto fiter = srVertex.find(neighbours[kk] );

				if ( fiter != srVertex.end() )
				{
					valid_neig.insert( *fiter );
				}
			}

			IndexType neig_siz = valid_neig.size();

			if(neig_siz<4)
			{
				Logger<<" neighbor size is small!.\n";
				continue; 
			}

			Matrix3X s_coord, t_coord;
			s_coord.setZero(3, neig_siz);
			t_coord.setZero(3, neig_siz);

			Sample& prevF = SampleSet::get_instance()[ sFrame ];

			IndexType vi = 0;
			IndexType kk = 0;

			for ( map<IndexType,CVertex*>::iterator neig_iter = valid_neig.begin();
				neig_iter != valid_neig.end(); ++neig_iter )
			{
				s_coord.col(kk) = curF.vertices_matrix().col(neig_iter->first );
				t_coord.col(kk) = prevF.vertices_matrix().col(neig_iter->second->prev_corr->vtx_id);
			}

			Matrix33 rot_mat;
			MatrixXX tran_vec;

			point2point(s_coord, t_coord, rot_mat, tran_vec);

			auto bias = rot_mat * curF.vertices_matrix().col(tgComIter->first) + tran_vec
				- prevF.vertices_matrix().col(tgComIter->second->prev_corr->vtx_id);
		}

		return totError;

	}

	// 分裂其中的对应块

	void split_patches(IndexType sFrame, set<IndexType>& sPatches,
		               IndexType tFrame, set<IndexType>& tPatches )
	{
		assert(sPatches.size() > 0 && tPatches.size() > 0);

		if ( sPatches.size() == 1 && tPatches.size() == 1)
		{
			Logger<<"不需要分裂,直接建立两块之间的对应!.\n";
			return;
		}

		map<IndexType,CVertex*> srVertex;

		set<IndexType>::iterator sIter = sPatches.begin();

		for (; sIter != sPatches.end(); sIter ++)
		{
			CLabel* temp = components_[sFrame].label_bucket[*sIter];
			srVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
		}

		map<IndexType,CVertex*> tgVertex;

		set<IndexType>::iterator tIter = tPatches.begin();

		for (; tIter != tPatches.end(); tIter++)
		{
			CLabel* temp = components_[tFrame].label_bucket[*tIter];
			tgVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
		}
		map<IndexType,CVertex*> toComVtx; // 记录那些点它们的对应会落在对应块上

		map<IndexType,CVertex*> backComVtx; // 相互记录

		bool isTo = true;
		computer_common(srVertex,tgVertex,toComVtx,isTo);

		isTo = false;
		computer_common(srVertex, tgVertex,backComVtx,isTo);


		//split 策略
		if ( sPatches.size() < tPatches.size()) // 先分裂source  
		{      
			// 保存一份在分割之前的label_bucket数据,先申请内存在赋值
	        //包含前传和后传两步骤

			CFrame beforSplit(components_[sFrame] );

			for (auto pIter = tPatches.begin(); pIter != tPatches.end(); pIter ++)   //A <---------------B
			{
				CLabel* tPatch = components_[tFrame].label_bucket[*pIter];

				map<IndexType,CVertex*> tPatchVtx(tPatch->vertex_bucket.begin(), tPatch->vertex_bucket.end() );//获得patch上的点集

			    map<IndexType,CVertex*> backComPatVtx;//该块中点在上帧的对应落在对方块上
				
				backComPatVtx.clear();

				computer_common(srVertex,tPatchVtx,backComPatVtx,false);

				if ( backComPatVtx.size() < 15)
				{
					Logger<<"顶点个数太少,不分裂.\n";
					continue;
				}else
				{
					//now split
					map<IndexType,CVertex*> prev_map; //在上一帧中对应的结点
					prev_map.clear();

					for (auto comIter = backComPatVtx.begin(); comIter != backComPatVtx.end(); comIter ++)
					{
						CVertex& prevV = *comIter->second->prev_corr;

						prev_map.insert(make_pair(prevV.vtx_id, &prevV) );

					}

					//now create new component
					IndexType new_label = components_[sFrame].label_bucket.size();
					components_[sFrame].label_bucket.push_back( (CLabel*)0 );
					CLabel* new_label_space = allocator_.allocate<CLabel>();
					CLabel* new_label_obj = new (new_label_space)CLabel;
					components_[sFrame].label_bucket[new_label]    = new_label_obj;          
					components_[sFrame].label_bucket[new_label]->label_id = new_label;
					components_[sFrame].label_bucket[new_label]->frame_parent = &components_[sFrame];

					//遍历srVertex所有点,以便判断是保持label不变还是move到新的label中

					auto to_remove = prev_map.begin();
					for(auto vIter = srVertex.begin(); vIter != srVertex.end() && to_remove != prev_map.end();)
					{
						IndexType vtx_idx = (*vIter).first;
						if (prev_map.find(vtx_idx) != prev_map.end() )
						{
							components_[sFrame].label_bucket[new_label]->vertex_bucket.insert( * vIter);
							vIter = srVertex.erase( vIter);
							components_[sFrame].label_of_vtx[vtx_idx] = new_label;
						}else if( vIter->second->next_corr->label_parent->label_id  == (*pIter) ) //对应点落在下帧的块中
						{
							components_[sFrame].label_bucket[new_label]->vertex_bucket.insert( * vIter);
							vIter = srVertex.erase( vIter);
							components_[sFrame].label_of_vtx[vtx_idx] = new_label;
						}else
						{
							vIter ++;
						}
					}

					for (auto vv = components_[sFrame].label_bucket[new_label]->vertex_bucket.begin();
						 vv != components_[sFrame].label_bucket[new_label]->vertex_bucket.end(); vv++) 
					{
						vv->second->label_parent = components_[sFrame].label_bucket[new_label];
					}

				}

			}// 遍历target 上组合块上的所有块,决定sr上的分裂情况

			// 用分裂前的label_burcket来指导target的分裂. ------------------>>
			for (auto sPIter = sPatches.begin(); sPIter != sPatches.end(); sPIter ++)
			{
				CLabel* temp = beforSplit.label_bucket[*sPIter];
				
				map<IndexType,CVertex*> sPatchVtx(temp->vertex_bucket.begin(),temp->vertex_bucket.end() );

				map<IndexType,CVertex*> toComPatVtx;

				toComPatVtx.clear();

				computer_common(sPatchVtx,tgVertex,toComPatVtx,true);

				if (toComPatVtx.size() < 15)
				{
					Logger<<"顶点个数太少,不分裂.\n";
					continue;
				}else
				{
					//start split
					map<IndexType,CVertex*> next_map;// 在下一帧中对应的结点
					next_map.clear();

					for (auto comIter = toComPatVtx.begin(); comIter != toComPatVtx.end(); comIter ++)
					{
						CVertex* prevV = comIter->second->next_corr;

						next_map.insert(make_pair(prevV->vtx_id,prevV) );
					}
                    
					//now create new component
					IndexType new_label = components_[tFrame].label_bucket.size();
					components_[tFrame].label_bucket.push_back( (CLabel*)0 );
					CLabel* new_label_space = allocator_.allocate<CLabel>();
					CLabel* new_label_obj = new (new_label_space)CLabel;
					components_[tFrame].label_bucket[new_label]    = new_label_obj;          
					components_[tFrame].label_bucket[new_label]->label_id = new_label;
					components_[tFrame].label_bucket[new_label]->frame_parent = &components_[tFrame];


					for(auto vIter = tgVertex.begin(); vIter != tgVertex.end();)
					{
						IndexType vtx_idx = (*vIter).first;
						if (next_map.find(vtx_idx) != next_map.end() )
						{
							components_[tFrame].label_bucket[new_label]->vertex_bucket.insert( * vIter);
							vIter = tgVertex.erase( vIter);
							components_[tFrame].label_of_vtx[vtx_idx] = new_label;
						}else if( vIter->second->prev_corr->label_parent->label_id  == (*sPIter) ) //对应点落在下帧的块中
						{
							components_[tFrame].label_bucket[new_label]->vertex_bucket.insert( * vIter);
							vIter = tgVertex.erase( vIter);
							components_[tFrame].label_of_vtx[vtx_idx] = new_label;
						}else
						{
							vIter ++;
						}
					}

					for (auto vv = components_[tFrame].label_bucket[new_label]->vertex_bucket.begin();
						vv != components_[tFrame].label_bucket[new_label]->vertex_bucket.end(); vv++) 
					{
						vv->second->label_parent = components_[tFrame].label_bucket[new_label];
					}
				}

			}

		}else //先分裂target   A----->B
		{
			CFrame beforSplit(components_[tFrame] );

			for (auto sPIter = sPatches.begin(); sPIter != sPatches.end(); sPIter ++)
			{
				CLabel* temp = components_[sFrame].label_bucket[*sPIter];

				map<IndexType,CVertex*> sPatchVtx(temp->vertex_bucket.begin(),temp->vertex_bucket.end() );

				map<IndexType,CVertex*> toComPatVtx;

				toComPatVtx.clear();

				computer_common(sPatchVtx,tgVertex,toComPatVtx,true);

				if (toComPatVtx.size() < 15)
				{
					Logger<<"顶点个数太少,不分裂.\n";
					continue;
				}else
				{
					//start split
					map<IndexType,CVertex*> next_map;// 在下一帧中对应的结点
					next_map.clear();

					for (auto comIter = toComPatVtx.begin(); comIter != toComPatVtx.end(); comIter ++)
					{
						CVertex* prevV = comIter->second->next_corr;

						next_map.insert(make_pair(prevV->vtx_id,prevV) );
					}

					//now create new component
					IndexType new_label = components_[tFrame].label_bucket.size();
					components_[tFrame].label_bucket.push_back( (CLabel*)0 );
					CLabel* new_label_space = allocator_.allocate<CLabel>();
					CLabel* new_label_obj = new (new_label_space)CLabel;
					components_[tFrame].label_bucket[new_label]    = new_label_obj;          
					components_[tFrame].label_bucket[new_label]->label_id = new_label;
					components_[tFrame].label_bucket[new_label]->frame_parent = &components_[tFrame];


					for(auto vIter = tgVertex.begin(); vIter != tgVertex.end();)
					{
						IndexType vtx_idx = (*vIter).first;
						if (next_map.find(vtx_idx) != next_map.end() )
						{
							components_[tFrame].label_bucket[new_label]->vertex_bucket.insert( * vIter);
							vIter = tgVertex.erase( vIter);
							components_[tFrame].label_of_vtx[vtx_idx] = new_label;
						}else if( vIter->second->prev_corr->label_parent->label_id  == (*sPIter) ) //对应点落在下帧的块中
						{
							components_[tFrame].label_bucket[new_label]->vertex_bucket.insert( * vIter);
							vIter = tgVertex.erase( vIter);
							components_[tFrame].label_of_vtx[vtx_idx] = new_label;
						}else
						{
							vIter ++;
						}
					}

					for (auto vv = components_[tFrame].label_bucket[new_label]->vertex_bucket.begin();
						vv != components_[tFrame].label_bucket[new_label]->vertex_bucket.end(); vv++) 
					{
						vv->second->label_parent = components_[tFrame].label_bucket[new_label];
					}
				}

			}
		}


	}

	void computer_common(map<IndexType,CVertex*>& srPatchesVtx, map<IndexType,CVertex*>& tgPatchesVtx, 
		                 map<IndexType,CVertex*>& comVtx,bool isTo)
	{
       

		if (isTo) // -->
		{
            map<IndexType,CVertex*>::iterator iter = srPatchesVtx.begin();

			for (; iter != srPatchesVtx.end(); iter ++)
			{
                if( tgPatchesVtx.find( iter->second->next_corr->vtx_id) != tgPatchesVtx.end())
				{
					comVtx.insert(*iter);
				}

			}

		}else// <--
		{
			map<IndexType,CVertex*>::iterator iter = tgPatchesVtx.begin();

			for (; iter != tgPatchesVtx.end(); iter++)
			{
				if ( srPatchesVtx.find(iter->second->prev_corr->vtx_id) != srPatchesVtx.end())
				{
					comVtx.insert(*iter);
				}
			}
		}

		return;
	}

	void buildSmapleKDtree(CFrame* smpCloud, Sample* smp)
	{
/*		smp = new Sample;*/

		Sample& ori = SampleSet::get_instance()[smpCloud->frame_id];

		map<IndexType,IndexType>::iterator vtxIter = smpCloud->label_of_vtx.begin();

		for (; vtxIter != smpCloud->label_of_vtx.end(); vtxIter ++)
		{
			Vertex& vtx = ori[vtxIter->first];

			PointType v( vtx.x(), vtx.y(), vtx.z() );

			ColorType cv(vtx.r(), vtx.g(), vtx.b(), vtx.alpha());

			NormalType nv(vtx.nx(), vtx.ny(), vtx.nz());

		    smp->add_vertex(v,nv,cv);
		}

		smp->build_kdtree();
		
		return; 
	}

	void point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec)
	{
		Eigen::Vector3f X_mean, Y_mean;

		for(int i=0; i<3; ++i) //计算两点云的均值
		{
			X_mean(i) = srCloud.row(i).sum()/srCloud.cols();
			Y_mean(i) = tgCloud.row(i).sum()/tgCloud.cols();
		}

		srCloud.colwise() -= X_mean;
		tgCloud.colwise() -= Y_mean;

		/// Compute transformation
		Eigen::Affine3f transformation;
		Eigen::Matrix3f sigma = srCloud * tgCloud.transpose();
		Eigen::JacobiSVD<Eigen::Matrix3f> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);
		if(svd.matrixU().determinant()*svd.matrixV().determinant() < 0.0)//contains reflection
		{
			Eigen::Vector3f S = Eigen::Vector3f::Ones(); S(2) = -1.0;
			transformation.linear().noalias() = svd.matrixV()*S.asDiagonal()*svd.matrixU().transpose();
		} else 
		{
			transformation.linear().noalias() = svd.matrixV()*svd.matrixU().transpose();//计算旋转矩阵
		}

		transVec = Y_mean - transformation.linear()*X_mean;
		rotMat = transformation.linear() ;

		srCloud.colwise() += X_mean;
		tgCloud.colwise() += Y_mean;

	}


	//0821
	void split_twoAjacent_graph_next(IndexType srFrame, IndexType tgFrame);

	void split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame);

	void split_nest_graph_prev(IndexType srFrame, IndexType tgFrame);

	void splitAllSquenceGraph();

	void getNextCorVtx(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor);

	IndexType checkNextLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor);

private:
	map<IndexType, CFrame> components_;
	map<IndexType,IndexType> frame_index_map_label_; 
	PoolAllocator allocator_;

	map<IndexType, HFrame> hier_componets_;
};


#endif