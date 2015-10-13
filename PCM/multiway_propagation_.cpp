#include"multiway_propagation_.h"


void DualwayPropagation::split_twoAjacent_graph_next(IndexType srFrame, IndexType tgFrame)
{
	//向前分裂
	
	//get the new graph of tgGrame--需要深拷贝

	IndexType tgGraSize = hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* oriGra = hier_componets_[tgFrame].hier_graph[tgGraSize - 1];
	LabelsGraph* new_graph = new LabelsGraph(*oriGra);
	//new_graph = oriGra;

	// 
	vector<HLabel* > new_label_bucket = hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1];
	//
	IndexType srGraSize = hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* srGraLat = hier_componets_[tgFrame].hier_graph[srGraSize - 1];



	for (;;)//"遍历srGraLat所有的边"Edge
	{
		//取出Edge的边界点属性,map<IndexType,CVertex*> edgeMapVtx;
		//根据edgeMapVtx的点在下一帧的对应点,找到它们在new_graph对应的节点node;
		//保存所有与node相关联的边的边属性, 用map<int, map<IndexType,CVertex*> >记录,其中的key保存与node相连的节点编号;
		//删除所有与node相关联的边;
		//增加一个新顶点new_node
		//把CLabel[node]中的点分类两类,根据其中的点对应到上一帧的标签情况,相同标签的点放在一起(***);
		//对两类赋标签,一组还是node,另外一组为新增加的标签;
		//对node 和new_node建立新的边关系(*****), 对新增加的边添加属性;      
	}

	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*srGraLat);

	for (EdgeIterator eit = ei.first; eit != ei.second; eit++)
	{
		EdgeDescriptor ed = *eit;

		GraphEdgeProperty& ep = (*srGraLat)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 3)
		{
			Logger<<"边上的顶点数太少,无法分裂.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//获得边界点在下一帧对应的块和对应点

		//遍历nodeId 上相连接的边

		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);

		VertexIterator nodeIter = (vi.first + nodeId);
        
		VertexDescriptor nodeDesc = *nodeIter;

        //节点对应的所有出边
		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(nodeDesc,*new_graph);

		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;

		set<GraphEdgeProperty> collapseEdges;

		for (OutEdgeIterator eit = nodeEiter.first; eit != nodeEiter.second; eit ++)
		{
			EdgeDescriptor nextEdgeD = *eit;

			 GraphEdgeProperty& nextEP = (*new_graph)[nextEdgeD];

			collapseEdges.insert(nextEP);

// 			map<IndexType,HVertex*> edgeProPoints;
// 
// 			auto ePsIt = nextEP.edgePoints.begin();
// 			edgeProPoints.insert(ePsIt->second.begin(),ePsIt->second.end() );
// 
// 
// 			VertexDescriptor srVtxD = boost::source(nextEdgeD,*new_graph);
// 
// 			VertexDescriptor tgVtxD = boost::target(nextEdgeD,*new_graph);
// 
// 			GraphVertexProperty& srVtxP = (*new_graph)[srVtxD];
// 
// 			GraphVertexProperty& tgVtxP = (*new_graph)[tgVtxD];
// 
// 			IndexType srId = srVtxP.index;
// 
// 			IndexType tgId = tgVtxP.index;
// 
// 			IndexType seIdKey = frame_index_to_key(srId,tgId);
// 			//记录边的起点和终点
// 
// 			recordColapseEdges[seIdKey] = edgeProPoints;


			//删除这条边
			boost::remove_edge(nextEdgeD,*new_graph);

		}

		//增加一个节点

		IndexType nSize = boost::num_vertices(*new_graph);

        GraphVertexProperty vp(nSize,-1,-1);

		boost::add_vertex(vp,*new_graph);


		//更新分割块信息,新增加的Label标号为nSize. 被分裂的点为nodeId
 
		IndexType new_label = nSize;
		
		new_label_bucket.push_back((HLabel*)0 );
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label_obj = new (new_label_space)HLabel;
		new_label_bucket[new_label] = new_label_obj;
		new_label_bucket[new_label]->label_id = new_label;
		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];
	
		HLabel* splitedLabel = new_label_bucket[nodeId];

		//对应回去,标签是起始点,则保持不变
		IndexType eS = ep.start_;
		IndexType eE = ep.end_;
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end();)
		{
			IndexType prev_id = iter->second->prev_corr->label_parent->label_id;
			IndexType vtx_id = iter->first;
			if (prev_id == eS)
			{
				Logger<<"不改变标签.\n";
				iter ++;
			}else if (prev_id == eE)
			{
				new_label_bucket[new_label]->vertex_bucket.insert(*iter);
				iter = splitedLabel->vertex_bucket.erase(iter);
				hier_componets_[tgFrame].label_of_vtx[vtx_id ] = new_label;
			}else
			{
				Logger<<"一些待确定label的点.\n";
			}

		}

		//对这两个点进行加边操作,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = nodeId;
		newEP.end_ = nSize;

		//还没对index赋值

		auto sBIter = newEP.edgePoints.begin();

		(*sBIter).second.insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());

		boost::add_edge(nodeId,nSize,newEP,*new_graph);

		//断定查找两个节点与其它节点进行连边操作只会出现 recordColapseEdges.size()次数.
		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
		{

			GraphEdgeProperty glueEdge;

			glueEdge = (*iter);
			IndexType  sVtxId = glueEdge.start_;

 			if (sVtxId < nodeId) //nodeid为终点
 			{
 				// 计算svtxId 和nodeid,nSize之间距离,判断终点节点是哪个.以便加边
 				glueEdge.end_ = nodeId;//或者 = nSize;
 
 			}else//nodeid为起点
 			{
 				//判断 eVtxId和nodeid,nSize之间距离,判断终点节点是哪个.以便加边
 			}
		}


	}//遍历每条边,每条边都会使得new_graph增加一个新的节点


	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//保存最新的graph

}

void DualwayPropagation::split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame)
{
	
	//首先更新srFrame的grpah.

}

void DualwayPropagation::split_nest_graph_prev(IndexType srFrame, IndexType tgFrame)
{
	if(srFrame == 1)
	{
		split_twoAjacent_graph_prev(srFrame,tgFrame);
		return;
	}

	while (srFrame > 0)
	{
		split_twoAjacent_graph_prev(srFrame,tgFrame);

		srFrame --;

		tgFrame --;
	}

	return;
}

void DualwayPropagation::splitAllSquenceGraph()
{
	for (auto cIter = components_.begin(); cIter != components_.end(); cIter++)
	{
        //最后一个节点不可以运算
		IndexType srFrame = cIter->first;
		IndexType tgFrame = srFrame + 1;

		split_twoAjacent_graph_next(srFrame,tgFrame );

		split_nest_graph_prev(srFrame,tgFrame);
	}

}

void DualwayPropagation::read_label_file_hier(char *filename)
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


		if ( hier_componets_.find(frame)== hier_componets_.end() )
		{
			hier_componets_.insert(make_pair(frame, HFrame()));
			hier_componets_[frame].frame_id = frame;	
			hier_componets_[frame].hier_label_bucket.resize(1);
		}
		if ( label >= hier_componets_[frame].hier_label_bucket[0].size() )
		{
			hier_componets_[frame].hier_label_bucket[0].resize( label+1 );
		}
		if (   nullptr==hier_componets_[frame].hier_label_bucket[0][label] )
		{
			HLabel* new_label_space = allocator_.allocate<HLabel>();
			HLabel* new_label = new (new_label_space)HLabel;
			hier_componets_[frame].hier_label_bucket[0][label] = new_label;
			hier_componets_[frame].hier_label_bucket[0][label]->frame_parent = &hier_componets_[frame];
			hier_componets_[frame].hier_label_bucket[0][label]->label_id = label;
		}
		HVertex* new_space = allocator_.allocate<HVertex>();
		HVertex* new_vtx = new (new_space)HVertex(vtx_idx, hier_componets_[frame].hier_label_bucket[0][label]);
		hier_componets_[frame].hier_label_bucket[0][label]->vertex_bucket.insert( make_pair(vtx_idx,new_vtx) );
		hier_componets_[frame].label_of_vtx.insert( make_pair(vtx_idx, label) );

	}
}

IndexType DualwayPropagation::checkNextLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor)
{
	assert(edgePs.size() > 0);

	set<IndexType > labelS;
	labelS.clear();

	for (auto iter = edgePs.begin(); iter != edgePs.end(); iter ++)
	{
		IndexType nextVtxLabid = iter->second->next_corr->label_parent->label_id;
		labelS.insert(nextVtxLabid);

		HVertex* nextP = iter->second->next_corr;
		IndexType corPId = iter->second->next_corr->vtx_id;
		edgePsCor[corPId] = nextP;

	}

	if (labelS.size() > 1)
	{
		Logger<<"边界点指向了下个块的两个Label!.\n";

	}

	return (*labelS.begin() );
}
void DualwayPropagation::getNextCorVtx(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor)
{
// 	for (auto iter = edgePs.begin(); iter != edgePs.end(); iter ++ )
// 	{
// 		HVertex& nextP = iter->sec
// 	}
}