#ifndef _CO_SEGMENT_H
#define _CO_SEGMENT_H

#include <QThread>
#include <basic_types.h>
#include <unordered_map>
#include <queue>
#include <vector>
#include <set>
#include "pool_allocator.h"
#include "sample_set.h"
#include "graph_cut_node.h"

#include "GraphMatching.h"

//  struct HVertex;
//  struct HLabel;
//  struct HFrame;
//  struct GraphVertexProperty;
//  struct GraphEdgeProperty;
//  struct GraphVertexProperty;
//  struct GraphEdgeProperty;


using namespace std;

// #define frame_index_to_key(f,i) ((f<<16)|i)
// #define frame_label_to_key(f,l) ((f<<16)|l)
// #define get_index_from_key(k) (k&0xffff)
// #define get_frame_from_key(k) (k>>16)

struct Component{
	IndexType frame;
	IndexType label;
	map<IndexType, IndexType> vtx_corr_next_frame;
	map<IndexType,IndexType> vtx_corr_prev_frame;
	Component(IndexType f, IndexType l):frame(f),label(l){}
};

struct MinTree
{
	void init(IndexType n)
	{
		parent.clear();
		rank.clear();
		parent.resize(n);
		rank.resize(n);
		for(IndexType i=0; i<n; ++i)
		{
			parent[i] = i;
			rank[i] = 0;
			union_find_set.insert( make_pair(i, set<IndexType>()) );
			union_find_set[i].insert(i);
		}
	}

	IndexType find(IndexType x)
	{
		while( x!=parent[x] )
			x = parent[x];

		return x;
	}

	void unite(IndexType x, IndexType y)
	{
		x = find(x);
		y = find(y);
		if(x==y)return;
		if(rank[x]<rank[y])
		{
			parent[x] = y;
			union_find_set[y].insert( union_find_set[x].begin(), union_find_set[x].end() );
			union_find_set.erase(x);

		}
		else
		{
			parent[y] = x;
			union_find_set[x].insert( union_find_set[y].begin(), union_find_set[y].end() );
			union_find_set.erase(y);
			if(rank[x]==rank[y])rank[x]++;
		}

	}

	bool same(IndexType x, IndexType y)
	{
		return find(x) == find(y);
	}

	vector<IndexType> parent;
	vector<IndexType> rank;
	map<IndexType, set<IndexType> > union_find_set;
};

class CoSegmentation
{

	struct compo_dist_node
	{
		compo_dist_node(IndexType _i, IndexType _j, ScalarType _d):i(_i),j(_j),dist(_d){}
		IndexType i,j;
		ScalarType dist;
	};

	struct compo_dist_greater_than
	{
		inline bool operator()( const compo_dist_node& lhs, const compo_dist_node& rhs )
		{
			return lhs.dist < rhs.dist;
		}
	};

public:
	//CoSegmentation(){}

	CoSegmentation(SampleSet& _smpSet, map<IndexType,HFrame>* oriComponents): smpSet(_smpSet),hier_componets(oriComponents)
	{
	}

	~CoSegmentation(){}

public:
	void read_data(char *label_name,char *corr_name);
	void read_label_file(char *filename);
	void read_corres_file(char *filename);
	void write_label_file(char *filename);
	void visualize_seg();
	void build_mini_seg();
	void compute();
	void calc_similarity();
	ScalarType similarity_between_component(IndexType i, IndexType j);
	void cluster_component();
	void co_segment();  
	void propagate_component();

	void calc_similarity_with_graph();

    //20151010
	void hierComponets2Components();
	void components2HierComponets();

public:
	ScalarType similarity_between_componet_with_graph(IndexType srFrame,IndexType tgFrame);
	void matchAndMergePatches();

public:
	IndexType compute_common( Component & c1, Component & c2, map<IndexType,IndexType> & common_part );
	void point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec);
	IndexType compute_common_back(Component & c1, Component & c2, map<IndexType,IndexType> & common_part );
	ScalarType similarity_between_component_back(IndexType i, IndexType j);
private:
		vector<Component> components_;
		typedef priority_queue<compo_dist_node, vector<compo_dist_node>, compo_dist_greater_than> Cqueue;
		Cqueue cq_;
		MinTree mint_;
		map<IndexType,IndexType> compo_fast_frame_label;
		map<IndexType,IndexType> compo_fast_frame_vtx;
		vector<IndexType> mini_seg;
		IndexType *dp_common_part;

public:
	SampleSet& smpSet;
	map<IndexType,HFrame>* hier_componets;
	
	static PoolAllocator allocator_;

};

// #undef frame_index_to_key
// #undef frame_label_to_key
// #undef get_index_from_key
// #undef get_frame_from_key

#endif