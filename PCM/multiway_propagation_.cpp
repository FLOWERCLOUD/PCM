#include"multiway_propagation_.h"


void DualwayPropagation::split_twoAjacent_graph_next(IndexType srFrame, IndexType tgFrame)
{
	//��ǰ����
	
	//get the new graph of tgGrame--��Ҫ���

	IndexType tgGraSize = hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* oriGra = hier_componets_[tgFrame].hier_graph[tgGraSize - 1];
	LabelsGraph* new_graph = new LabelsGraph(*oriGra);
	//new_graph = oriGra;

	// 
	vector<HLabel* > new_label_bucket = hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1];
	//
	IndexType srGraSize = hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* srGraLat = hier_componets_[tgFrame].hier_graph[srGraSize - 1];



	for (;;)//"����srGraLat���еı�"Edge
	{
		//ȡ��Edge�ı߽������,map<IndexType,CVertex*> edgeMapVtx;
		//����edgeMapVtx�ĵ�����һ֡�Ķ�Ӧ��,�ҵ�������new_graph��Ӧ�Ľڵ�node;
		//����������node������ıߵı�����, ��map<int, map<IndexType,CVertex*> >��¼,���е�key������node�����Ľڵ���;
		//ɾ��������node������ı�;
		//����һ���¶���new_node
		//��CLabel[node]�еĵ��������,�������еĵ��Ӧ����һ֡�ı�ǩ���,��ͬ��ǩ�ĵ����һ��(***);
		//�����ำ��ǩ,һ�黹��node,����һ��Ϊ�����ӵı�ǩ;
		//��node ��new_node�����µı߹�ϵ(*****), �������ӵı��������;      
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
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��

		//����nodeId �������ӵı�

		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);

		VertexIterator nodeIter = (vi.first + nodeId);
        
		VertexDescriptor nodeDesc = *nodeIter;

        //�ڵ��Ӧ�����г���
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
// 			//��¼�ߵ������յ�
// 
// 			recordColapseEdges[seIdKey] = edgeProPoints;


			//ɾ��������
			boost::remove_edge(nextEdgeD,*new_graph);

		}

		//����һ���ڵ�

		IndexType nSize = boost::num_vertices(*new_graph);

        GraphVertexProperty vp(nSize,-1,-1);

		boost::add_vertex(vp,*new_graph);


		//���·ָ����Ϣ,�����ӵ�Label���ΪnSize. �����ѵĵ�ΪnodeId
 
		IndexType new_label = nSize;
		
		new_label_bucket.push_back((HLabel*)0 );
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label_obj = new (new_label_space)HLabel;
		new_label_bucket[new_label] = new_label_obj;
		new_label_bucket[new_label]->label_id = new_label;
		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];
	
		HLabel* splitedLabel = new_label_bucket[nodeId];

		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
		IndexType eS = ep.start_;
		IndexType eE = ep.end_;
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end();)
		{
			IndexType prev_id = iter->second->prev_corr->label_parent->label_id;
			IndexType vtx_id = iter->first;
			if (prev_id == eS)
			{
				Logger<<"���ı��ǩ.\n";
				iter ++;
			}else if (prev_id == eE)
			{
				new_label_bucket[new_label]->vertex_bucket.insert(*iter);
				iter = splitedLabel->vertex_bucket.erase(iter);
				hier_componets_[tgFrame].label_of_vtx[vtx_id ] = new_label;
			}else
			{
				Logger<<"һЩ��ȷ��label�ĵ�.\n";
			}

		}

		//������������мӱ߲���,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = nodeId;
		newEP.end_ = nSize;

		//��û��index��ֵ

		auto sBIter = newEP.edgePoints.begin();

		(*sBIter).second.insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());

		boost::add_edge(nodeId,nSize,newEP,*new_graph);

		//�϶����������ڵ��������ڵ�������߲���ֻ����� recordColapseEdges.size()����.
		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
		{

			GraphEdgeProperty glueEdge;

			glueEdge = (*iter);
			IndexType  sVtxId = glueEdge.start_;

 			if (sVtxId < nodeId) //nodeidΪ�յ�
 			{
 				// ����svtxId ��nodeid,nSize֮�����,�ж��յ�ڵ����ĸ�.�Ա�ӱ�
 				glueEdge.end_ = nodeId;//���� = nSize;
 
 			}else//nodeidΪ���
 			{
 				//�ж� eVtxId��nodeid,nSize֮�����,�ж��յ�ڵ����ĸ�.�Ա�ӱ�
 			}
		}


	}//����ÿ����,ÿ���߶���ʹ��new_graph����һ���µĽڵ�


	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//�������µ�graph

}

void DualwayPropagation::split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame)
{
	
	//���ȸ���srFrame��grpah.

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
        //���һ���ڵ㲻��������
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
		Logger<<"�߽��ָ�����¸��������Label!.\n";

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