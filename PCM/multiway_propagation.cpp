#include"multiway_propagation.h"


#define  INF 1000000
#define  COUT_DEBUG 0

#define frame_index_to_key(f,i) ((f<<16)|i)
#define frame_label_to_key(f,l) ((f<<16)|l)
#define get_index_from_key(k) (k&0xffff)
#define get_frame_from_key(k) (k>>16)




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


void DualwayPropagation::read_label_file_coseg(char *label_file_name)
{
	FILE* in_file = fopen(label_file_name, "r");
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


void DualwayPropagation::read_corr_file_coseg(char* corr_file_name)
{
	FILE* in_file = fopen(corr_file_name,"r");
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


	//ɾ��label_buckert��Ϊ�յ�Ԫ��
	for (auto fiter = hier_componets_.begin(); fiter != hier_componets_.end(); ++ fiter)
	{
		IndexType frame_id = fiter->first;
		vector<HLabel*> fHierBucket = hier_componets_[frame_id].hier_label_bucket[0];

		vector<HLabel*>::iterator  hlabelIter = fHierBucket.begin();

		map<IndexType,IndexType> vtxIndex;
		IndexType vtxLevel = 0;
		for (; hlabelIter != fHierBucket.end(); )
		{
			if (*hlabelIter == NULL)
			{
				hlabelIter = fHierBucket.erase(hlabelIter);
			}else
			{
				IndexType labelId = (*hlabelIter)->label_id;
				vtxIndex[labelId] = vtxLevel;

				++ vtxLevel;
				++ hlabelIter;
			}
		}

		hier_componets_[frame_id].hier_label_bucket[0] = fHierBucket;

		hier_componets_[frame_id].hier_label_vtxBucket_index.resize(1);

		hier_componets_[frame_id].hier_label_vtxBucket_index[0] = vtxIndex;

	}

}

void DualwayPropagation::getEdgeVertexs( IndexType _CFrameId ,IndexType lLabelId ,IndexType rLabelId , map<IndexType, map<IndexType ,HVertex*> >& _edgepoints )
{
	map<IndexType,HVertex*> twoLabelsVtx;

	IndexType sr_labelLevel = hier_componets_[_CFrameId].hier_label_vtxBucket_index[0][lLabelId];
	IndexType tg_labelLevel = hier_componets_[_CFrameId].hier_label_vtxBucket_index[0][rLabelId];


	twoLabelsVtx.insert(hier_componets_[_CFrameId].hier_label_bucket[0][sr_labelLevel]->vertex_bucket.begin(),
		hier_componets_[_CFrameId].hier_label_bucket[0][sr_labelLevel]->vertex_bucket.end() );

	twoLabelsVtx.insert(hier_componets_[_CFrameId].hier_label_bucket[0][tg_labelLevel]->vertex_bucket.begin(),
		hier_componets_[_CFrameId].hier_label_bucket[0][tg_labelLevel]->vertex_bucket.end() );

// 	twoLabelsVtx.insert(hier_componets_[_CFrameId].hier_label_bucket[0][lLabelId]->vertex_bucket.begin(),
// 		                hier_componets_[_CFrameId].hier_label_bucket[0][lLabelId]->vertex_bucket.end() );
// 
// 	twoLabelsVtx.insert(hier_componets_[_CFrameId].hier_label_bucket[0][rLabelId]->vertex_bucket.begin(),
// 		                hier_componets_[_CFrameId].hier_label_bucket[0][rLabelId]->vertex_bucket.end() );

	const IndexType k = 3;

	IndexType neighbours[k];

	ScalarType dist[k];

	IndexType i = 0;
	map<IndexType,HVertex*> epoints;
	for (auto iter = twoLabelsVtx.begin(); iter != twoLabelsVtx.end(); iter ++)
	{
		IndexType ori_idx = iter->first;
		IndexType map_id = o2d_map[ori_idx];

		downSample->neighbours(map_id,k,neighbours,dist);
		set<IndexType> neilabel_size;
		neilabel_size.clear();
		for (IndexType kk = 0; kk < k; kk++)
		{
			IndexType neig_downId = neighbours[kk];
			IndexType oriNeig_id = d2o_map[neig_downId];
			IndexType neig_label = hier_componets_[_CFrameId].label_of_vtx[oriNeig_id];
			neilabel_size.insert(neig_label);
		}
		
		if( neilabel_size.size() > 1)
		{
			//map<IndexType, map<IndexType ,HVertex*> >& _edgepoints
			epoints[ori_idx] = (*iter).second;
			IndexType seKey = frame_label_to_key(lLabelId,rLabelId);
			_edgepoints[seKey] = epoints;
			
		} ;//�����Ǳ߽��
	}

	return;
}



void DualwayPropagation::getEdgeVertexs2( IndexType _CFrameId , distanPriQueue& _PriQuemap, map<IndexType, map<IndexType ,HVertex*> >& _edgepoints )
{
	map<IndexType,HVertex*> edgedps;
	pointdistance pd;


	IndexType sr_labelLevel,tg_labelLevel;
	for( IndexType i  = 10 ; i >0 ; --i)
	{
		if (_PriQuemap.empty())
		{
			Logger<<"Queue empty.\n";
			break;
		}
		pd = _PriQuemap.top();
		_PriQuemap.pop();

	    sr_labelLevel = hier_componets_[_CFrameId].hier_label_vtxBucket_index[0][pd.labelofvtx1_];
		tg_labelLevel = hier_componets_[_CFrameId].hier_label_vtxBucket_index[0][pd.labelofvtx2_];

		HVertex* tmpVtx;


		tmpVtx = hier_componets_[_CFrameId].hier_label_bucket[0][sr_labelLevel]->vertex_bucket[pd.vtx1Id_];
		edgedps[ pd.vtx1Id_] = tmpVtx;

// 		tmpVtx = hier_componets_[_CFrameId].hier_label_bucket[0][pd.labelofvtx2_]->vertex_bucket[pd.vtx2Id_];
// 		edgedps[ pd.vtx2Id_] = tmpVtx;



	}

	IndexType seKey = frame_label_to_key(sr_labelLevel, tg_labelLevel);
	//IndexType seKey = frame_label_to_key( pd.labelofvtx1_, pd.labelofvtx2_);

	_edgepoints[ seKey] = edgedps;

}

void DualwayPropagation::buildKdTree( IndexType _cframeId)
{
	if (downSample != NULL)
	{
		delete downSample;
	}

	downSample =new Sample();
	Sample& smp = SampleSet::get_instance()[_cframeId];
	auto iter = hier_componets_[_cframeId].label_of_vtx.begin();
	IndexType i = 0;
	for(; iter != hier_componets_[_cframeId].label_of_vtx.end(); ++iter,i++ )
	{
		IndexType vtx_id = (*iter).first;
		Vertex& vtx = smp[vtx_id];

		PointType v( vtx.x(), vtx.y(), vtx.z() );
		ColorType cv(vtx.r(), vtx.g(), vtx.b(), vtx.alpha());
		NormalType nv(vtx.nx(), vtx.ny(), vtx.nz());

		downSample->add_vertex(v,nv,cv);
		o2d_map[vtx_id] = i;
		d2o_map[i] = vtx_id;

	}

	downSample->build_kdtree();

}


void DualwayPropagation:: init_labeles_graph_hier()
{

	Logger<<"  Begin initialize graphs.\n";

	for (auto citer = hier_componets_.begin(); citer!=hier_componets_.end(); citer++)
	{

		IndexType lbsize = citer->second.hier_label_bucket.size();
		assert( lbsize > 0);
		IndexType nodeSize = citer->second.hier_label_bucket[lbsize -1].size();//.size();

		LabelsGraph* new_labelGraph_space = allocator_.allocate<LabelsGraph>();

		LabelsGraph* new_labelGraph = new (new_labelGraph_space)LabelsGraph;

		GraphVertexProperty gvp;

		map<IndexType,IndexType> labelLevel = citer->second.hier_label_vtxBucket_index[0];

		auto labIndex = labelLevel.begin();

		IndexType gep_count  = 0;

		for( IndexType i = 0 ; i< nodeSize && labIndex != labelLevel.end(); ++i,++ labIndex)
		{
			gvp.index = i ;
			gvp.prev =  -1;
			gvp.next = -1;
			gvp.label_id = labIndex->first;
			add_vertex( gvp , *new_labelGraph);
		}

		for (IndexType i = 0; i < nodeSize-1; i ++)
		{
			for (IndexType j = i + 1; j < nodeSize; j++)
			{
				GraphEdgeProperty   gep;

				distanPriQueue PriQue;   //
				while ( ! PriQue.empty() ) PriQue.pop(); 

				HLabel* fir = citer->second.hier_label_bucket[lbsize -1][i];

				HLabel* sec = citer->second.hier_label_bucket[lbsize -1][j];

				map<IndexType,HVertex*> fir_vtx = fir->vertex_bucket;

				map<IndexType,HVertex*> sec_vtx = sec->vertex_bucket;

				//����������֮�����̾���
				ScalarType minDis = 1e5;
				SampleSet& sample_set = SampleSet::get_instance();
				map<IndexType,HVertex*>::iterator biter1 , eiter1 ,biter2,eiter2;
				eiter1 = fir_vtx.end();
				eiter2 = sec_vtx.end();
				ScalarType dia ;

				for( biter1 = fir_vtx.begin() ;biter1 != eiter1 ;++biter1)
				{
					for( biter2 = sec_vtx.begin() ;biter2 != eiter2 ;++biter2)
					{
						IndexType index1 = biter1->first;
						Sample& s = sample_set[citer->first];
						dia = s.getBox().diag();

						PointType point1 =  s.vertices_matrix().col(index1);

						IndexType index2 = biter2->first;

						PointType point2 =  s.vertices_matrix().col(index2);

						ScalarType distance = (point1 - point2).norm();
						if(distance < minDis)minDis = distance;
						IndexType label1 = hier_componets_[citer->first].label_of_vtx[ index1];
						IndexType label2 = hier_componets_[citer->first].label_of_vtx[ index2];
						if( distance < dia*0.1){
							pointdistance tpd( index1 ,  label1 , index2 ,label2 ,distance );
							PriQue.push( tpd);
						}

					}

				}

				if (minDis <  dia * 0.1 )//0.0588
				{
					getEdgeVertexs2( citer->first ,PriQue ,gep.edgePoints);

					//�߽�㼯����̫��Ҳ����Ҫ���,������������.

					gep.index = gep_count;
					if( i < j)
					{
						gep.start_ = i;
						gep.end_ = j;
						boost::add_edge(i,j,gep ,*new_labelGraph);
						++gep_count;
					}else if( i> j)
					{
						gep.start_ = j;
						gep.end_ = i;
						boost::add_edge( j, i ,gep ,*new_labelGraph);
						++gep_count;
					}else
					{
						// i == j ʱ������
					}

				}

			}

		}

		//citer->second.hier_graph.clear();

		citer->second.hier_graph.push_back( new_labelGraph);

	}//֡��������

	Logger<<"  End initialize graphs.\n";
}

void DualwayPropagation::split_twoAjacent_graph_next(IndexType srFrame, IndexType tgFrame)
{
	//��ǰ����

	//get the new graph of tgGrame--��Ҫ���


	Logger<<" .......\n";
	Logger<<"  Start next split.\n";
	IndexType tgGraSize = hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* oriGra = hier_componets_[tgFrame].hier_graph[tgGraSize - 1];
	LabelsGraph* new_graph = new LabelsGraph(*oriGra);

	//new_graph = oriGra;

	//vector<HLabel* > new_label_bucket =  hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1];
	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1] );


	//
	IndexType gLevel = 0;

	IndexType srGraSize = hier_componets_[srFrame].hier_graph.size();

	gLevel  = srGraSize - 1;//��ȡ���µĲ�
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

    IndexType labParentsize = tgGraSize + 1; //���ɵĲ���

	Logger<<srFrame<<"֡�ĵ�"<<gLevel<<"��߽�ָ�"<<tgFrame<<"��"<<tgGraSize - 1 <<"��"<<endl;

	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*srGraLat);

	for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
	{
		EdgeDescriptor ed = *eit;

		GraphEdgeProperty& ep = (*srGraLat)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 1)
		{
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType newGraphEdgeSize = new_graph->m_edges.size();

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��
		//IndexType nodeId = edgeCorrNextVtx.size();


		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		Logger<<"  �ߵ����Ϊ"<<eS<<"�յ�Ϊ"<<eE<<endl;

		IndexType recordS = 0;       
		IndexType recordE = 0;

		IndexType vtxBSzie = splitedLabel->vertex_bucket.size();
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
		{
			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //

			if (prev_id == eS)
			{
				recordS ++;

			}else if(prev_id == eE)
			{
				recordE ++;
			}
		}

		ScalarType ration = (ScalarType)(recordS + recordE)/vtxBSzie;

		//�����ѳ����ĵ������һ�����ݺ���,��ñ߲����ѱ�

		if ( recordE < 5 || recordS < 5 )
		{
			Logger<<"�߽��̫����,����Ҫ����.\n";
			continue;
		}

		if (ration < 0.2)
		{
			Logger<<"Unmark���ֵ̫��,��ʱ������.\n";
			continue;
		}

		//����nodeId �������ӵı�

		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);

		VertexIterator nodeIter = (vi.first + nodeId);

		VertexDescriptor nodeDesc = *nodeIter;

		//�ڵ��Ӧ�����г���
		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(nodeDesc,*new_graph);

		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;

		set<GraphEdgeProperty> collapseEdges;

		OutEdgeIterator oit,nextIt;

		oit = nodeEiter.first;

		for (nextIt = oit; oit != nodeEiter.second; oit = nextIt )
		{
			++nextIt;

			EdgeDescriptor nextEdgeD = *oit;

			GraphEdgeProperty& nextEP = (*new_graph)[nextEdgeD];

			collapseEdges.insert(nextEP);
			
			boost::remove_edge(*oit,*new_graph);//ɾ��������
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
 
// 		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���

		 map<IndexType,HVertex*> unMakePs;

		// Sample& smp = SampleSet::get_instance()[tgFrame];
		 //ScalarType egThreshold = smp.getBox().diag();

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //�ò������µ�label_parent��ַ.
			IndexType vtx_id = iter->first;

// 			HVertex* prevCorVtx = iter->second->prev_corr->next_corr;
// 			IndexType prevCorId = prevCorVtx->vtx_id;
// 			PointType oriCoor = smp.vertices_matrix().col(vtx_id);
// 			PointType prevnextCoor = smp.vertices_matrix().col(prevCorId); 
// 			ScalarType dist = (oriCoor - prevnextCoor).norm();

			if (prev_id == eS)
			{

					iter->second->label_parent.resize(tgGraSize + 1);
					iter->second->label_parent[tgGraSize] =  new_label_bucket[nodeId];
					++iter;

			}else if(prev_id == eE)
			{

					iter->second->label_parent.resize(tgGraSize + 1);
					iter->second->label_parent[tgGraSize] = new_label_bucket[new_label] ;
					new_label_bucket[new_label]->vertex_bucket.insert(*iter);
					hier_componets_[tgFrame].label_of_vtx[vtx_id ] = new_label;
					iter = splitedLabel->vertex_bucket.erase(iter);


 			}else//һЩ��ȷ��label�ĵ�
 			{
                   unMakePs.insert(*iter);
                   iter = splitedLabel->vertex_bucket.erase(iter);
 			}

		}

		//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
		determinateUnmarkPoints(tgFrame,unMakePs,new_label_bucket,nodeId,new_label,tgGraSize);


		//������������мӱ߲���,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = nodeId;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;

		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		//��û��index��ֵ

		newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());


		boost::add_edge(nodeId,nSize,newEP,*new_graph);

		//�϶����������ڵ��������ڵ�������߲���ֻ����� recordColapseEdges.size()����.
		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
		{

			GraphEdgeProperty glueEdge;

			glueEdge = (*iter);
			IndexType  sEdgeId = glueEdge.start_;
			IndexType  eEdgeId = glueEdge.end_;

			map<IndexType,HVertex*>  edgePoints;

			auto bIter = glueEdge.edgePoints.begin();
			edgePoints.insert(bIter->second.begin(),bIter->second.end() );

			map<IndexType,HVertex*> startVtx = new_label_bucket[sEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> endVtx = new_label_bucket[eEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> nodeVtx = new_label_bucket[nodeId]->vertex_bucket;
			map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;

			ScalarType minNode = 0.0;
			ScalarType minSize = 0.0;


			if (sEdgeId < nodeId) //nodeidΪ�յ�
			{
				minDistBeTwoParts(tgFrame,startVtx,nodeVtx,minNode);
				minDistBeTwoParts(tgFrame,startVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.end_ = nodeId;
				}else
				{
					glueEdge.end_ = nSize;
				}
							
			}else//nodeidΪ���
			{
				minDistBeTwoParts(tgFrame,endVtx,nodeVtx,minNode);
				minDistBeTwoParts(tgFrame,endVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.start_ = nodeId;
					glueEdge.end_ = eEdgeId;

				}else
				{
					glueEdge.start_ = eEdgeId;
					glueEdge.end_ = nSize;
				}

			}

			IndexType eKey = frame_index_to_key(glueEdge.start_,glueEdge.end_);

			glueEdge.edgePoints[eKey] = edgePoints;

			boost::add_edge(glueEdge.start_,glueEdge.end_,glueEdge,*new_graph);

		}//����collapse�ı�


	}//����ÿ����,ÿ���߶���ʹ��new_graph����һ���µĽڵ�


	

	checkPsNewLabelParentPtr(new_label_bucket,labParentsize);//next dirction


	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//�������µ�graph

	hier_componets_[tgFrame].hier_label_vtxBucket_index.push_back(labelIndex);

	Logger<<"  Start next split.\n";
	Logger<<" .......\n";
}

void DualwayPropagation::split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame)
{

	//��ȡ��Ҫ���µ�graph
	Logger<<" .......\n";
	Logger<<"  Start prev split.\n";

	IndexType srGraphSize = hier_componets_[srFrame].hier_graph.size();
	IndexType tgGraphSize = hier_componets_[tgFrame].hier_graph.size();

	assert(srGraphSize > 0 && tgGraphSize > 0);

	LabelsGraph* oriSpGra = hier_componets_[srFrame].hier_graph[srGraphSize - 1];
	LabelsGraph* shouldSplitGraph = new LabelsGraph(*oriSpGra);

	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[srFrame].hier_label_bucket[srGraphSize - 1] );

	//��ȡָ���ָ��graph
	LabelsGraph* guideSplitGraph;
	IndexType graphLevel = 0;

	assert(tgGraphSize > 1);

	if (tgGraphSize > 2)
	{
		graphLevel = tgGraphSize - 1;
	}else
	{
		graphLevel = tgGraphSize - 2;
	}
	guideSplitGraph = hier_componets_[tgFrame].hier_graph[graphLevel];

	//�Ա߽������ȼ�����




	// guideSplitGraph��ÿ��������һ�ηָ�
	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*guideSplitGraph);

	for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
	{
		EdgeDescriptor ed = *eit;

		GraphEdgeProperty& ep = (*guideSplitGraph)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 1)
		{
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		Logger<<tgFrame<<"֡�ı����"<<ep.start_<<"���յ�"<<ep.end_<<"��ʼ����"<<endl;

		map<IndexType,HVertex*> edgeCorrPrevVtx;

		IndexType nSize = boost::num_vertices(*shouldSplitGraph);

		IndexType newGraphEdgeSize = shouldSplitGraph->m_edges.size(); //Ϊ�˸������ӵı�������

		bool isSplit = true;

		IndexType edgePsCorNode = 0;

	    isSplit = checkPrevLabelBucket(edgePoints,edgeCorrPrevVtx,edgePsCorNode);

		if ( edgePsCorNode < 0 || edgePsCorNode > nSize - 1)
		{
			Logger<<"�߽��ҵ��Ŀ�Խ��!.\n";
			break;
		}

		if (!isSplit)
		{
			Logger<<"�߽��ָ����������,���ָ�.\n";
			continue;
		}

		HLabel* splitedLabel = new_label_bucket[edgePsCorNode]; //������Ҫ���ѵĿ�

		IndexType curEdgeStart = ep.start_;
		IndexType curEdgeEnd   = ep.end_;
		IndexType strCorPsSzie = 0;
		IndexType endCorPsSize = 0;

		//����ǰ��һ���򵥵�Ԥ�ж�
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
		{
			if (iter->second->next_corr->label_parent[graphLevel] != NULL)
			{

				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; 

				if (nextVtx_label == curEdgeStart)
				{
					strCorPsSzie ++;

				}else if(nextVtx_label == curEdgeEnd)
				{
					endCorPsSize ++;
				}
			}

		}

		IndexType  vtxBSize = splitedLabel->vertex_bucket.size();

		ScalarType ration = (ScalarType)(strCorPsSzie + endCorPsSize)/vtxBSize;

		if ( strCorPsSzie < 10 || endCorPsSize < 10)
		{
			Logger<<"�߽�㿿���߽�(��Ӧ�����),����Ҫ����.\n";
			continue;
		}

		if ( ration < 0.2)
		{
			Logger<<"������ǧ��?���˰�!.\n";
			continue;
		}

		//��eit �߿�ʼ����
		pair<VertexIterator, VertexIterator> vi = boost::vertices(*shouldSplitGraph);

		VertexIterator curNodeIter = (vi.first + edgePsCorNode);

		VertexDescriptor curNodeDesc = *curNodeIter;

		//�ڵ��Ӧ�����г���
		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(curNodeDesc,*shouldSplitGraph);

		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;

		set<GraphEdgeProperty> collapseEdges;

		OutEdgeIterator oit,nextIt;

		oit = nodeEiter.first;

		for (nextIt = oit; oit != nodeEiter.second; oit = nextIt )
		{
			++nextIt;

			EdgeDescriptor nextEdgeD = *oit;

			GraphEdgeProperty& nextEP = (*shouldSplitGraph)[nextEdgeD];

			collapseEdges.insert(nextEP);

			boost::remove_edge(*oit,*shouldSplitGraph);//ɾ��������
		}

		//����һ���ڵ�

		GraphVertexProperty vp(nSize,-1,-1);
		boost::add_vertex(vp,*shouldSplitGraph);

		//���·ָ����Ϣ,�����ӵ�Label���ΪnSize. �����ѵĵ�ΪedgePsCorNode
		IndexType new_label = nSize;
		new_label_bucket.push_back((HLabel*)0 );
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label_obj = new (new_label_space)HLabel;
		new_label_bucket[new_label] = new_label_obj;
		new_label_bucket[new_label]->label_id = new_label;
		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];

		//��ʼ����
		map<IndexType,HVertex*> unMakePs;

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			if (iter->second->next_corr->label_parent[graphLevel])
			{
				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; //�ò������µ�label_parent��ַ.

				IndexType vtx_id = iter->first;

				if (nextVtx_label == curEdgeStart)
				{
					iter->second->label_parent.resize(srGraphSize + 1);
					iter->second->label_parent[srGraphSize] = new_label_bucket[edgePsCorNode];
					++iter;

				}else if(nextVtx_label == curEdgeEnd)
				{
					iter->second->label_parent.resize(srGraphSize + 1);

					iter->second->label_parent[srGraphSize] = new_label_bucket[new_label];

					new_label_bucket[new_label]->vertex_bucket.insert(*iter);

					hier_componets_[srFrame].label_of_vtx[vtx_id ] = new_label;

					iter = splitedLabel->vertex_bucket.erase(iter);

				}else//һЩ��ȷ��label�ĵ�
				{
					unMakePs.insert(*iter);
					iter = splitedLabel->vertex_bucket.erase(iter);
				}
			}

		} //������Ҫ���ѿ��ÿ����

			//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
		determinateUnmarkPoints(srFrame,unMakePs,new_label_bucket,edgePsCorNode,new_label,srGraphSize);

		//������������мӱ߲���,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = edgePsCorNode;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;
		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		newEP.edgePoints[edgeKey].insert(edgeCorrPrevVtx.begin(),edgeCorrPrevVtx.end());
		boost::add_edge(edgePsCorNode,nSize,newEP,*shouldSplitGraph);


		//�϶����������ڵ��������ڵ�������߲���ֻ����� recordColapseEdges.size()����.
		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
		{
			GraphEdgeProperty glueEdge;
			glueEdge = (*iter);
			IndexType  sEdgeId = glueEdge.start_;
			IndexType  eEdgeId = glueEdge.end_;

			map<IndexType,HVertex*>  edgePoints;

			auto bIter = glueEdge.edgePoints.begin();
			edgePoints.insert(bIter->second.begin(),bIter->second.end() );

			map<IndexType,HVertex*> startVtx = new_label_bucket[sEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> endVtx = new_label_bucket[eEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> nodeVtx = new_label_bucket[edgePsCorNode]->vertex_bucket;
			map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;

			ScalarType minNode = 0.0;
			ScalarType minSize = 0.0;


			if (sEdgeId < edgePsCorNode) //nodeidΪ�յ�
			{
				minDistBeTwoParts(srFrame,startVtx,nodeVtx,minNode);
				minDistBeTwoParts(srFrame,startVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.end_ = edgePsCorNode;
				}else
				{
					glueEdge.end_ = nSize;
				}

			}else//nodeidΪ���
			{
				minDistBeTwoParts(srFrame,endVtx,nodeVtx,minNode);
				minDistBeTwoParts(srFrame,endVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.start_ = edgePsCorNode;
					glueEdge.end_ = eEdgeId;

				}else
				{
					glueEdge.start_ = eEdgeId;
					glueEdge.end_ = nSize;
				}

			}

			IndexType eKey = frame_index_to_key(glueEdge.start_,glueEdge.end_);

			glueEdge.edgePoints[eKey] = edgePoints;

			boost::add_edge(glueEdge.start_,glueEdge.end_,glueEdge,*shouldSplitGraph);

		}//����collapse�ı�

	} //���������ָ�ͼ��ÿ����

	checkPsNewLabelParentPtr(new_label_bucket,srGraphSize + 1);

	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[srFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[srFrame].hier_graph.push_back(shouldSplitGraph);//�������µ�graph

	hier_componets_[srFrame].hier_label_vtxBucket_index.push_back(labelIndex);

    Logger<<"  End prev split.\n";
	Logger<<" .......\n";
}

void DualwayPropagation::split_twoAjacent_graph_prev_order(IndexType srFrame, IndexType tgFrame)
{
// 	//��ȡ��Ҫ���µ�graph
// 	Logger<<" .......\n";
// 	Logger<<"  Start prev split.\n";
// 
// 	IndexType srGraphSize = hier_componets_[srFrame].hier_graph.size();
// 	IndexType tgGraphSize = hier_componets_[tgFrame].hier_graph.size();
// 
// 	assert(srGraphSize > 0 && tgGraphSize > 0);
// 
// 	LabelsGraph* oriSpGra = hier_componets_[srFrame].hier_graph[srGraphSize - 1];
// 	LabelsGraph* shouldSplitGraph = new LabelsGraph(*oriSpGra);
// 
// 	vector<HLabel* > new_label_bucket;
// 	copyLabelBucket(new_label_bucket,hier_componets_[srFrame].hier_label_bucket[srGraphSize - 1] );
// 
// 	//��ȡָ���ָ��graph
// 	LabelsGraph* guideSplitGraph;
// 	IndexType graphLevel = 0;
// 
// 	assert(tgGraphSize > 1);
// 
// 	if (tgGraphSize > 2)
// 	{
// 		graphLevel = tgGraphSize - 1;
// 	}else
// 	{
// 		graphLevel = tgGraphSize - 2;
// 	}
// 
// 	guideSplitGraph = hier_componets_[tgFrame].hier_graph[graphLevel];
// 
// 	// guideSplitGraph��ÿ��������һ�ηָ�
// 	//pair<EdgeIterator,EdgeIterator> ei = boost::edges(*guideSplitGraph);
// 
// 	generateOrderPrevEdges(srFrame,tgFrame);
// 
// 	while (!orderedEdgeQ.empty())
// 	{
// 		EdgeSplitOrder oEdge = orderedEdgeQ.top();
// 		orderedEdgeQ.pop();
// 
// 		EdgeDescriptor ed = oEdge.EdgeDec;
// 
// 		GraphEdgeProperty& ep = (*guideSplitGraph)[ed];
// 
// 		map<IndexType,HVertex*> edgePoints;
// 
// 		auto ePsIt = ep.edgePoints.begin();
// 
// 		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  
// 
// 		if (edgePoints.size() < 3)
// 		{
// 			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
// 			continue;
// 		}
// 
// 		map<IndexType,HVertex*> edgeCorrNextVtx;
// 
// 		IndexType newGraphEdgeSize = shouldSplitGraph->m_edges.size();
// 
// 		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��
// 		//IndexType nodeId = edgeCorrNextVtx.size();
// 
// 
// 		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
// 		HLabel* splitedLabel = new_label_bucket[nodeId];
// 
// 		IndexType eS = ep.start_;
// 		IndexType eE = ep.end_;
// 
// 		Logger<<"  �ߵ����Ϊ"<<eS<<"�յ�Ϊ"<<eE<<endl;
// 
// 		if ( oEdge.srCorNum < 10 || oEdge.tgCorNum < 10 )
// 		{
// 			Logger<<"�߽��̫����,����Ҫ����.\n";
// 			continue;
// 		}
// 
// 		if ( oEdge.unMarkedRation < 0.2)
// 		{
// 			Logger<<"Unmark���ֵ̫��,��ʱ������.\n";
// 			continue;
// 		}
// 
// 		//��eit �߿�ʼ����
// 		pair<VertexIterator, VertexIterator> vi = boost::vertices(*shouldSplitGraph);
// 
// 		VertexIterator curNodeIter = (vi.first + edgePsCorNode);
// 
// 		VertexDescriptor curNodeDesc = *curNodeIter;
// 
// 		//�ڵ��Ӧ�����г���
// 		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(curNodeDesc,*shouldSplitGraph);
// 
// 		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;
// 
// 		set<GraphEdgeProperty> collapseEdges;
// 
// 		OutEdgeIterator oit,nextIt;
// 
// 		oit = nodeEiter.first;
// 
// 		for (nextIt = oit; oit != nodeEiter.second; oit = nextIt )
// 		{
// 			++nextIt;
// 
// 			EdgeDescriptor nextEdgeD = *oit;
// 
// 			GraphEdgeProperty& nextEP = (*shouldSplitGraph)[nextEdgeD];
// 
// 			collapseEdges.insert(nextEP);
// 
// 			boost::remove_edge(*oit,*shouldSplitGraph);//ɾ��������
// 		}
// 
// 		//����һ���ڵ�
// 
// 		GraphVertexProperty vp(newGraphEdgeSize,-1,-1);
// 		boost::add_vertex(vp,*shouldSplitGraph);
// 
// 		//���·ָ����Ϣ,�����ӵ�Label���ΪnSize. �����ѵĵ�ΪedgePsCorNode
// 		IndexType new_label = newGraphEdgeSize;
// 		new_label_bucket.push_back((HLabel*)0 );
// 		HLabel* new_label_space = allocator_.allocate<HLabel>();
// 		HLabel* new_label_obj = new (new_label_space)HLabel;
// 		new_label_bucket[new_label] = new_label_obj;
// 		new_label_bucket[new_label]->label_id = new_label;
// 		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];
// 
// 		//��ʼ����
// 		map<IndexType,HVertex*> unMakePs;
// 
// 		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
// 		{
// 			if (iter->second->next_corr->label_parent[graphLevel])
// 			{
// 				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; //�ò������µ�label_parent��ַ.
// 
// 				IndexType vtx_id = iter->first;
// 
// 				if (nextVtx_label == curEdgeStart)
// 				{
// 					iter->second->label_parent.resize(srGraphSize + 1);
// 					iter->second->label_parent[srGraphSize] = new_label_bucket[edgePsCorNode];
// 					++iter;
// 
// 				}else if(nextVtx_label == curEdgeEnd)
// 				{
// 					iter->second->label_parent.resize(srGraphSize + 1);
// 
// 					iter->second->label_parent[srGraphSize] = new_label_bucket[new_label];
// 
// 					new_label_bucket[new_label]->vertex_bucket.insert(*iter);
// 
// 					hier_componets_[srFrame].label_of_vtx[vtx_id ] = new_label;
// 
// 					iter = splitedLabel->vertex_bucket.erase(iter);
// 
// 				}else//һЩ��ȷ��label�ĵ�
// 				{
// 					unMakePs.insert(*iter);
// 					iter = splitedLabel->vertex_bucket.erase(iter);
// 				}
// 			}
// 
// 		} //������Ҫ���ѿ��ÿ����
// 
// 		//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
// 		determinateUnmarkPoints(srFrame,unMakePs,new_label_bucket,edgePsCorNode,new_label,srGraphSize);
// 
// 		//������������мӱ߲���,
// 		//node<--->new_node
// 		GraphEdgeProperty newEP;
// 		newEP.start_ = edgePsCorNode;
// 		newEP.end_ = newGraphEdgeSize;
// 		newEP.index = newGraphEdgeSize;
// 		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
// 		newEP.edgePoints[edgeKey].insert(edgeCorrPrevVtx.begin(),edgeCorrPrevVtx.end());
// 		boost::add_edge(edgePsCorNode,newGraphEdgeSize,newEP,*shouldSplitGraph);
// 
// 
// 		//�϶����������ڵ��������ڵ�������߲���ֻ����� recordColapseEdges.size()����.
// 		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
// 		{
// 			GraphEdgeProperty glueEdge;
// 			glueEdge = (*iter);
// 			IndexType  sEdgeId = glueEdge.start_;
// 			IndexType  eEdgeId = glueEdge.end_;
// 
// 			map<IndexType,HVertex*>  edgePoints;
// 
// 			auto bIter = glueEdge.edgePoints.begin();
// 			edgePoints.insert(bIter->second.begin(),bIter->second.end() );
// 
// 			map<IndexType,HVertex*> startVtx = new_label_bucket[sEdgeId]->vertex_bucket;
// 			map<IndexType,HVertex*> endVtx = new_label_bucket[eEdgeId]->vertex_bucket;
// 			map<IndexType,HVertex*> nodeVtx = new_label_bucket[edgePsCorNode]->vertex_bucket;
// 			map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;
// 
// 			ScalarType minNode = 0.0;
// 			ScalarType minSize = 0.0;
// 
// 
// 			if (sEdgeId < edgePsCorNode) //nodeidΪ�յ�
// 			{
// 				minDistBeTwoParts(srFrame,startVtx,nodeVtx,minNode);
// 				minDistBeTwoParts(srFrame,startVtx,nSizeVtx,minSize);
// 
// 				if (minNode < minSize)
// 				{
// 					glueEdge.end_ = edgePsCorNode;
// 				}else
// 				{
// 					glueEdge.end_ = nSize;
// 				}
// 
// 			}else//nodeidΪ���
// 			{
// 				minDistBeTwoParts(srFrame,endVtx,nodeVtx,minNode);
// 				minDistBeTwoParts(srFrame,endVtx,nSizeVtx,minSize);
// 
// 				if (minNode < minSize)
// 				{
// 					glueEdge.start_ = edgePsCorNode;
// 					glueEdge.end_ = eEdgeId;
// 
// 				}else
// 				{
// 					glueEdge.start_ = eEdgeId;
// 					glueEdge.end_ = nSize;
// 				}
// 
// 			}
// 
// 			IndexType eKey = frame_index_to_key(glueEdge.start_,glueEdge.end_);
// 
// 			glueEdge.edgePoints[eKey] = edgePoints;
// 
// 			boost::add_edge(glueEdge.start_,glueEdge.end_,glueEdge,*shouldSplitGraph);
// 
// 	}

}

void DualwayPropagation::split_nest_graph_prev(IndexType startFrame,IndexType srFrame, IndexType tgFrame)
{
	if(srFrame == startFrame)
	{
		split_twoAjacent_graph_prev(srFrame,tgFrame);
		return;
	}

	while (srFrame > startFrame -1 && (startFrame > -1) )
	{
		split_twoAjacent_graph_prev(srFrame,tgFrame);

		srFrame --;

		tgFrame --;
	}

	return;
}

void DualwayPropagation::splitAllSquenceGraph(IndexType iterN)
{
	IndexType fSize = hier_componets_.size();

	map<IndexType,HFrame>::iterator  cIter= hier_componets_.begin();

	map<IndexType,HFrame>::iterator  cEnd = hier_componets_.end();

	//++cIter;
//  	IndexType startF = 6;
//  	while (startF -- > 0)
//  	{
//  		++cIter;
//  	}
	//cIter;
	--cEnd;

	IndexType iterNum = 0;

	IndexType startFrameId = cIter->first;



	for (; cIter != cEnd && iterNum < iterN; ++ cIter++, ++iterNum)
	{
		//���һ���ڵ㲻��������
		IndexType srFrame = cIter->first;
		IndexType tgFrame = srFrame + 1;

		split_twoAjacent_graph_next(srFrame,tgFrame );

		//show_corresponding(srFrame);

		//split_twoAjacent_graph_next_order(srFrame,tgFrame );

		split_nest_graph_prev(startFrameId,srFrame,tgFrame);

	}

	//show_corresponding(1);

	//���ÿ֡�ķָ����

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); ++ fIter)
	{
		IndexType gLevel = fIter->second.hier_label_bucket.size();
		IndexType labeSize = fIter->second.hier_label_bucket[gLevel - 1].size();
		Logger<<"  ��"<<fIter->first<<"֡���ָ��"<<labeSize<<"����.\n";
	}

}

void DualwayPropagation::buildMainPatchMatching()
{

	Logger<<"  Build Graph Matching.\n";

	//auto bIter = hier_componets_.begin();
	//auto eIter = hier_componets_.end();
	//eIter --;
	//IndexType itNum = 0;
	//IndexType itTot = 8;
	//for (; bIter != eIter && itNum < itTot; bIter ++, itNum ++)
	//{
	//	IndexType fId = bIter->first;

	//	Logger<<fId<<"֡��"<<fId + 1<<"֡�Ķ�Ӧ.\n";

	//	buildNextPatchCorr(fId,fId + 1);

	//	//buildPrevPatchCorr(fId,fId + 1);
	//}


	Logger<<" Using Cosegmentation Labels to build Patches Correspondences.\n";

	GraphMatch pairFrameSimilar(SampleSet::get_instance(), hier_componets_, 1, 2);
	pairFrameSimilar.buildPatchCorrespondenceByLabel();

	//pairFrameSimilar.calculateSimilar2Frame();

	pairFrameSimilar.setMatchset(0,0);

	Logger<<"  End Graph Matching.\n";
}

void DualwayPropagation::buildNextPatchCorr(IndexType srFrame,IndexType tgFrame)
{
	auto vv = hier_componets_.find(tgFrame);

	if (vv == hier_componets_.end())
	{
		Logger<<" No next frame.\n";
		return;
	}

	GraphMatch pairFrameSimilar(SampleSet::get_instance(), hier_componets_, srFrame, tgFrame);

	pairFrameSimilar.calculateSimilar2Frame();
}

void DualwayPropagation::buildPrevPatchCorr(IndexType srFrame,IndexType tgFrame)
{

}
IndexType DualwayPropagation::checkNextLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor)
{
	assert(edgePs.size() > 0);

	set<IndexType > labelS;
	labelS.clear();

	for (auto iter = edgePs.begin(); iter != edgePs.end(); iter ++)
	{
		IndexType corPId = iter->second->next_corr->vtx_id;
		//

		//IndexType TestnextVtxLabid = hier_components_[2].label_of_vtx[corPId];
		IndexType labParSize = iter->second->next_corr->label_parent.size();
		IndexType updateLevel = labParSize - 1;

		IndexType nextVtxLabid = iter->second->next_corr->label_parent[updateLevel]->label_id;

		labelS.insert(nextVtxLabid);

		HVertex* nextP = iter->second->next_corr;

		edgePsCor[corPId] = nextP;

	}

	if (labelS.size() > 1)
	{
		Logger<<"�߽��ָ�����¸��������Label!.\n";

	}

	return (*labelS.begin() );
}

bool DualwayPropagation::checkPrevLabelBucket(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor, IndexType& isSplit)
{
	assert(edgePs.size() > 0);

	set<IndexType > labelS;
	labelS.clear();

	for (auto iter = edgePs.begin(); iter != edgePs.end(); iter ++)
	{
		IndexType labParSize =  iter->second->prev_corr->label_parent.size();

		IndexType updateLevel = labParSize - 1;

		IndexType corPId = iter->second->prev_corr->vtx_id;

		IndexType nextVtxLabid = iter->second->prev_corr->label_parent[updateLevel]->label_id;

		labelS.insert(nextVtxLabid);

		HVertex* nextP = iter->second->prev_corr;

		edgePsCor[corPId] = nextP;

	}

	if (labelS.size() > 1)
	{
		Logger<<"�߽��ָ�����¸��������Label!.\n";

	}

	isSplit = (* labelS.begin() );

	return (labelS.size() == 1);
}

void DualwayPropagation::getNextCorVtx(map<IndexType,HVertex*>& edgePs, map<IndexType,HVertex*>& edgePsCor)
{
	// 	for (auto iter = edgePs.begin(); iter != edgePs.end(); iter ++ )
	// 	{
	// 		HVertex& nextP = iter->sec
	// 	}
}
void DualwayPropagation::read_data(char *label_name,char *corr_name)
{
	//read_label_file(label_name); //��ԭʼ�ĺ���

	//read_label_file_hier(label_name); //��ȡ�����Ϣ

	//read_label_file_coseg(label_name); //��ȡ���ָ��label�ļ�


    //201501010

  	read_label_file_coseg(label_name); //��ȡ���ָ��label�ļ�,����֮֡�䲻����label����

  	read_corr_file_coseg(corr_name);

}



void DualwayPropagation::wirteSplitGraphLables(std::string filename)
{
	/*	FILE* outfile = fopen(filename.c_str(),"w");
	for ( auto frame_iter = hier_componets_.begin();
	frame_iter != hier_componets_.end();
	frame_iter++ )
	{
	HFrame& frame = frame_iter->second;
	IndexType frame_idx = frame_iter->first;
	for ( auto label_iter = frame.label_of_vtx.begin(); 
	label_iter!=frame.label_of_vtx.end(); label_iter++  )
	{
	IndexType vtx_id = label_iter->first;
	IndexType label_idx = label_iter->second;
	fprintf(outfile, "%d %d %d\n", frame_idx, label_idx, vtx_id);

	}
	}

	fclose(outfile);*/



	FILE* outfile = fopen(filename.c_str(),"w");

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
	{
		IndexType gLevel = fIter->second.hier_label_bucket.size();//������߲��label_bucket
		IndexType fId = fIter->first;
		vector<HLabel*>& label_buctet = fIter->second.hier_label_bucket[gLevel - 1];

		for (auto lIter = label_buctet.begin(); lIter != label_buctet.end(); ++ lIter)
		{

			IndexType pId = (*lIter)->label_id;
			
			map<IndexType,HVertex*>& vtx_bucket = (*lIter)->vertex_bucket;

			for (auto vIt = vtx_bucket.begin(); vIt != vtx_bucket.end(); ++vIt)
			{
				IndexType vtx_id = vIt->first;

				fprintf(outfile, "%d %d %d\n", fId, pId, vtx_id);
			}
		}
	}

	fclose(outfile);
}


void DualwayPropagation::init_labeles_graph()
{
	for (auto citer = components_.begin(); citer!=components_.end(); citer++)
	{
		IndexType nodeSize = citer->second.label_bucket.size();

		LabelsGraph* new_labelGraph_space = allocator_.allocate<LabelsGraph>();

		LabelsGraph* new_labelGraph = new (new_labelGraph_space)LabelsGraph;


		new_labelGraph->added_vertex(nodeSize);

		for (IndexType i = 0; i < nodeSize-1; i ++)
		{
			for (IndexType j = i + 1; j < nodeSize; j++)
			{

				CLabel* fir = citer->second.label_bucket[i];

				CLabel* sec = citer->second.label_bucket[j];

				map<IndexType,CVertex*> fir_vtx = fir->vertex_bucket;

				map<IndexType,CVertex*> sec_vtx = sec->vertex_bucket;

				//����������֮�����̾���
				ScalarType minDis = 1000000;
				SampleSet& sample_set = SampleSet::get_instance();
				map<IndexType,CVertex*>::iterator biter1 , eiter1 ,biter2,eiter2;
				eiter1 = fir_vtx.end();
				eiter2 = sec_vtx.end();
				ScalarType dia ;
				for( biter1 = fir_vtx.begin() ;biter1 != eiter1 ;++biter1){
					for( biter2 = sec_vtx.begin() ;biter2 != eiter2 ;++biter2){
						IndexType index1 = biter1->first;
						Sample& s = sample_set[citer->first];
						dia = s.getBox().diag();

						PointType point1 =  s.vertices_matrix().col(index1);

						IndexType index2 = biter2->first;

						PointType point2 =  s.vertices_matrix().col(index2);

						ScalarType distance = (point1 - point2).norm();
						if(distance < minDis)minDis = distance;

					}

				}

				if (minDis <  dia * 0.0588 )
				{

					boost::add_edge(i,j,*new_labelGraph);
				}

			}
		}

		citer->second.labelGraph = new_labelGraph;

	}
}

void DualwayPropagation::minDistBeTwoParts(IndexType cFrame, map<IndexType,HVertex*>& fPart,map<IndexType,HVertex*>& sPart, ScalarType& misDis)
{

	Sample& smp = SampleSet::get_instance()[cFrame];

	ScalarType diag = smp.getBox().diag();

	vector<ScalarType> totDis;

	for (auto fIter = fPart.begin(); fIter != fPart.end(); ++fIter)
	{
		for (auto sIter = sPart.begin(); sIter != sPart.end(); ++sIter)
		{
			PointType fVtx = smp.vertices_matrix().col(fIter->first);
			PointType sVtx = smp.vertices_matrix().col(sIter->first);

			ScalarType dis = (fVtx - sVtx).norm();

			if (dis < 0.5 * diag)
			{
				totDis.push_back(dis);
			}
		}
	}

	if (totDis.size() == 0)
	{
		misDis =  0.5*diag;
	}else
	{
	    misDis = * min_element(totDis.begin(),totDis.end() );
	}


}

void DualwayPropagation::determinateUnmarkPoints(IndexType cFrame,map<IndexType,HVertex*>& unMarkPs,vector<HLabel*> oriLabelBucket,IndexType nodeId,IndexType newLabe,IndexType tgSize)
{
	if (unMarkPs.size() < 1)
	{
		Logger<<"û��Unmark ��.\n";
		return;
	}

	Sample& smp = SampleSet::get_instance()[cFrame];

	for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)
	{
		PointType pCoor = smp.vertices_matrix().col(iter->first);

// 		ScalarType d2node = p2PatchAvgDis(cFrame,pCoor,oriLabelBucket[nodeId]->vertex_bucket);//���ȡֵ����Сֵ
// 		ScalarType d2size = p2PatchAvgDis(cFrame,pCoor,oriLabelBucket[newLabe]->vertex_bucket );


// 		ScalarType d2node = p2PatchMinDis(cFrame,pCoor,oriLabelBucket[nodeId]->vertex_bucket);//�����ϸ���С
// 		ScalarType d2size = p2PatchMinDis(cFrame,pCoor,oriLabelBucket[newLabe]->vertex_bucket );


		//����߾���
		ScalarType d2node = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[nodeId]->vertex_bucket);//�����ϸ���С
		ScalarType d2size = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[newLabe]->vertex_bucket );

		if (d2node <= d2size)
		{

			oriLabelBucket[nodeId]->vertex_bucket.insert(*iter);

			hier_componets_[cFrame].label_of_vtx[iter->first] = nodeId;

			iter->second->label_parent.resize(tgSize + 1);
			iter->second->label_parent[tgSize] = oriLabelBucket[nodeId];

		}else
		{

			oriLabelBucket[newLabe]->vertex_bucket.insert(*iter);

			hier_componets_[cFrame].label_of_vtx[iter->first] = newLabe;

			//iter->second->label_parent.push_back(oriLabelBucket[newLabe] );
			iter->second->label_parent.resize(tgSize + 1);
			iter->second->label_parent[tgSize] = oriLabelBucket[newLabe] ;
		}

		iter = unMarkPs.erase(iter); 
	}

}
ScalarType DualwayPropagation::p2PatchAvgDis(IndexType cFrame, PointType& pCoor,map<IndexType,HVertex*>& parthPs)
{
	ScalarType avgDis = 1e5;

	ScalarType randDis = 0.0;

	IndexType pSize = parthPs.size();

	assert(pSize > 0);

	IndexType randPsSize = 100;
	
	Sample& smp = SampleSet::get_instance()[cFrame];


	auto eIter = parthPs.end();

	while(randPsSize -- > 0 )
	{
	    auto bIter = parthPs.begin();

		IndexType randPId = rand()%pSize;

		std::advance(bIter,randPId);

		if ( bIter !=  eIter)
		{
			IndexType tgPId = bIter->first;

			PointType tgCoor = smp.vertices_matrix().col(tgPId);

			randDis = (pCoor - tgCoor).norm();

			//avgDis += (pCoor - tgCoor).norm(); //ƽ��ֵ������

			if(randDis < avgDis)
			{
				avgDis = randDis;
			}
		}
	}


	return  1e2* avgDis;

}

ScalarType DualwayPropagation::p2PatchMinDis(IndexType cFrame, PointType& pCoor,map<IndexType,HVertex*>& parthPs)
{
	ScalarType avgDis = 1e5;

	ScalarType randDis = 0.0;

	IndexType pSize = parthPs.size();

	assert(pSize > 0);

	Sample& smp = SampleSet::get_instance()[cFrame];

	auto bIter = parthPs.begin();

	auto eIter = parthPs.end();

	for (; bIter != eIter; ++bIter)
	{
		IndexType tgPId = bIter->first;

		PointType tgCoor = smp.vertices_matrix().col(tgPId);

		randDis = (pCoor - tgCoor).norm();

		if(randDis < avgDis)
		{
			avgDis = randDis;
		}

	}

	return  1e2* avgDis;
}

ScalarType DualwayPropagation::p2PatchGeoDis(IndexType cFrame, HVertex& oriP,map<IndexType,HVertex*>& parthPs)
{

	ScalarType avgDis = 1e5;

	ScalarType randDis = 0.0;

	IndexType pSize = parthPs.size();

	assert(pSize > 0);

	IndexType srPId = oriP.vtx_id;
	IndexType nodeId = hier_componets_[cFrame].gId_of_vtx[srPId];

	ScalarType geoDis ;

	PCloudGraph* pg = hier_componets_[cFrame].pcGraph;

	IndexType nSize = boost::num_vertices(*pg);

	//����������
	pair<pcVertexIterator, pcVertexIterator> vi = boost::vertices(*pg);

	pcVertexIterator nodeIter = (vi.first + nodeId);

	pcVertexDescriptor nodeDesc = *nodeIter;

	vector<ScalarType> djDis(boost::num_vertices(*pg), 0.);
	vector<pcVertexDescriptor> parents(boost::num_vertices(*pg) );
	//djDis.resize();

	 auto p_map = boost::make_iterator_property_map(&parents[0], boost::get(boost::vertex_index, *pg));
	 auto w_map = boost::get(&PCEdgeProperty::dist,*pg);
	 auto d_map = boost::make_iterator_property_map(&djDis[0],boost::get(boost::vertex_index,*pg));

	boost::dijkstra_shortest_paths(*pg,nodeDesc,boost::weight_map(w_map).predecessor_map(p_map).distance_map(d_map) );



	auto bIter = parthPs.begin();

	auto eIter = parthPs.end();


	for (; bIter != eIter; ++bIter)
	{
		IndexType tgPId = bIter->first;
		IndexType tgNodeId = hier_componets_[cFrame].gId_of_vtx[tgPId];
		
		randDis = djDis[tgNodeId];

		if(randDis < avgDis)
		{
			avgDis = randDis;
		}

	}

	return  avgDis;
}


void DualwayPropagation::copyLabelBucket(vector<HLabel*>& leftLabels, const vector<HLabel*>& oriLabels)
{
	assert(oriLabels.size() > 0);

	for (auto iter = oriLabels.begin(); iter != oriLabels.end(); iter ++)
	{
		IndexType LId = (*iter)->label_id;

		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label = new (new_label_space)HLabel((*iter)->label_id,(*iter)->frame_parent,(*iter)->vertex_bucket,(*iter)->prev_corr,(*iter)->next_corr);

// 		for (auto cvIter = new_label->vertex_bucket.begin(); cvIter!= new_label->vertex_bucket.end(); cvIter++)
// 		{
// 			cvIter->second->label_parent = new_label;
// 		}

		leftLabels.push_back(new_label);
	}
}

// ȷ��ÿ�������labparentָ����Ŀ������ͬ
void DualwayPropagation::checkPsNewLabelParentPtr(vector<HLabel*> oriLabelBucket,IndexType labParSize)
{
	IndexType vId = 0;
	for (auto vIter = oriLabelBucket.begin(); vIter != oriLabelBucket.end(); vIter ++,vId ++)
	{
		 map<IndexType,HVertex*>& vtxBucket = (*vIter)->vertex_bucket;

		 HVertex* fVtx = (*vtxBucket.begin()).second;

		 IndexType vlfSize = fVtx->label_parent.size();

		 if ( vlfSize < labParSize)
		 {
			 for (auto lIter =vtxBucket.begin(); lIter != vtxBucket.end(); lIter ++ )
			 {
				 (*lIter).second->label_parent.push_back(oriLabelBucket[vId] );
			 }
		 }

	}
}

void DualwayPropagation::smoothAfterSplit()
{
	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter++)
	{
		smoothSingleFrame(fIter->first);
	}

}
void DualwayPropagation::smoothSingleFrame(IndexType frame_idx)//û���޸ĸ���label�������������㼯
{
	Sample& orig_smp = SampleSet::get_instance()[frame_idx];
	Sample* downsmp_ptr = new Sample;
	map<IndexType, IndexType> idx_mapper;
	IndexType i=0;
	for ( auto viter = hier_componets_[frame_idx].label_of_vtx.begin();
		viter != hier_componets_[frame_idx].label_of_vtx.end();
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

	const IndexType k = 30;
	IndexType neighbours[k];
	ScalarType dist[k];
	i = 0;

	//ֻ�������һ���label_bucket����
	map<IndexType, IndexType> old_label_map = hier_componets_[frame_idx].label_of_vtx;

	IndexType graphLevel = hier_componets_[frame_idx].hier_label_vtxBucket_index.size();
	IndexType hLevel_idx = graphLevel - 1;
	//IndexType beforLabelSize = hier_componets_[frame_idx].hier_label_bucket[hLevel_idx].size();

	IndexType beforLabelSize = 25;

	for ( i=0; i<idx_mapper.size(); i++)
	{
		vector<int> label_count( beforLabelSize, 0 );
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

		IndexType curVtxlabelIdx = hier_componets_[frame_idx].hier_label_vtxBucket_index[hLevel_idx][cur_vtx_label];

		if ( max_freq_label!=cur_vtx_label )
		{
			IndexType maxLabelIndex = hier_componets_[frame_idx].hier_label_vtxBucket_index[hLevel_idx][max_freq_label];
			//change label
			auto vv = hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][curVtxlabelIdx]->vertex_bucket.find(real_vtx_idx);
			
			assert(vv!=hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][curVtxlabelIdx]->vertex_bucket.end());

			(*vv).second->label_parent[hLevel_idx] = hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][maxLabelIndex ];
	 
			hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][maxLabelIndex]->vertex_bucket.insert(*vv);

			hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][curVtxlabelIdx]->vertex_bucket.erase(vv);

			hier_componets_[frame_idx].label_of_vtx[real_vtx_idx] = max_freq_label;

		}
	}

	delete downsmp_ptr;

	//for ( i=0; i<idx_mapper.size(); i++)
	//{
	//	vector<int> label_count( beforLabelSize, 0 );
	//	IndexType real_vtx_idx = idx_mapper[i];
	//	downsmp_ptr->neighbours(i, k, neighbours, dist );
	//	for ( IndexType neig_id = 0; neig_id<k; neig_id++ )
	//	{
	//		IndexType real_neig_id = idx_mapper[neighbours[neig_id]];
	//		IndexType neig_label = old_label_map[real_neig_id];
	//		label_count[neig_label]++;
	//	}
	//	IndexType max_freq_label = max_element(label_count.begin(), label_count.end()) - label_count.begin();
	//	IndexType cur_vtx_label = old_label_map[real_vtx_idx];

	//	IndexType curVtxlabelIdx = hier_componets_[frame_idx].hier_label_vtxBucket_index[graphLevel][cur_vtx_label];

	//	if ( max_freq_label!=cur_vtx_label )
	//	{
	//		//change label
	//		auto vv = hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][cur_vtx_label]->vertex_bucket.find(real_vtx_idx);

	//		assert(vv!=hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][cur_vtx_label]->vertex_bucket.end());

	//		(*vv).second->label_parent[graphLevel - 1] = hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][max_freq_label ];

	//		hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][max_freq_label]->vertex_bucket.insert(*vv);

	//		hier_componets_[frame_idx].hier_label_bucket[hLevel_idx][cur_vtx_label]->vertex_bucket.erase(vv);

	//		hier_componets_[frame_idx].label_of_vtx[real_vtx_idx] = max_freq_label;

	//	}
	//}




}

void DualwayPropagation::buildSquenceUniqLabel(std::string filename)
{
	map<IndexType,IndexType> patchTrajId;//����ͬ��patch����ͬ��id
	map<IndexType,bool>      isPatchTrav;//��¼���Ƿ���ʹ�--key = frameId  & labelId
	
	IndexType totLabelId = 0;
	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
	{
		IndexType gLevel = fIter->second.hier_label_bucket.size();//������߲��label_bucket
		IndexType fId = fIter->first;
		vector<HLabel*>& label_buctet = fIter->second.hier_label_bucket[gLevel - 1];
        
		IndexType lId = 0;
		for (auto lIter = label_buctet.begin(); lIter != label_buctet.end(); lIter ++, lId ++)
		{
			IndexType flKey = frame_label_to_key(fId,lId);
			if (isPatchTrav[flKey])
			{
				continue;
			}

			isPatchTrav[flKey] = true;

			patchTrajId[flKey] = totLabelId;

			HLabel* nextPatchPtr = label_buctet[lId]->next_corr;

			while(nextPatchPtr != NULL)
			{
				IndexType nFId = nextPatchPtr->frame_parent->frame_id;

				IndexType nLId = nextPatchPtr->label_id;

				IndexType nextflKey = frame_label_to_key(nFId,nLId);

				isPatchTrav[nextflKey] = true;

				patchTrajId[nextflKey] = totLabelId;

				nextPatchPtr = nextPatchPtr->next_corr;
			}

			totLabelId ++;
		}
	}

	//���ȫ�ֵ�Label
	writeSquenceLabel(patchTrajId,filename);

}

void DualwayPropagation::writeSquenceLabel(map<IndexType,IndexType>& patchTrajId,std::string filename)
{
		FILE* outfile = fopen(filename.c_str(),"w");

		for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
		{
			IndexType gLevel = fIter->second.hier_label_bucket.size();//������߲��label_bucket
			IndexType fId = fIter->first;
			vector<HLabel*>& label_buctet = fIter->second.hier_label_bucket[gLevel - 1];

			IndexType lId = 0;
			for (auto lIter = label_buctet.begin(); lIter != label_buctet.end(); lIter ++, lId ++)
			{
				IndexType flKey = frame_label_to_key(fId,lId);
				 auto totId = patchTrajId.find(flKey);
				 IndexType pId = 0;
				if (totId != patchTrajId.end() )
				{
					pId = patchTrajId[flKey];
				}else
				{
					pId =  15;
				}
				
				map<IndexType,HVertex*>& vtx_bucket = (*lIter)->vertex_bucket;

				for (auto vIt = vtx_bucket.begin(); vIt != vtx_bucket.end(); ++vIt)
				{
					IndexType vtx_id = vIt->first;
					fprintf(outfile, "%d %d %d\n", fId, pId, vtx_id);
				}
			}
		}

		fclose(outfile);
}

void DualwayPropagation::mergePatchesAfterCoSeg() //0831
{
	
    //����merge��ֹ���� ,��ȡCo-segmentation��label�ļ�֮������ú���

	buildPatchCorrespondenceByLabel();

	IndexType itNum = 5;
	IndexType i = 0;

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end() &&  i < itNum ; ++ fIter,++i)
	{
		IndexType frameId = fIter->first;
		if (hier_componets_.find(frameId + 1) != hier_componets_.end() )
		{
			IndexType combineNum = 1;//�ϲ�����

			while (combineNum -- > 0) 
			{
				set<IndexType> srBestSet,tgBestSet;
				srBestSet.clear();
				tgBestSet.clear();

				IndexType srLevel = hier_componets_[frameId].hier_graph.size();
				IndexType tgLevel = hier_componets_[frameId + 1].hier_graph.size();

				GraphMatch pairFrameSimilar(SampleSet::get_instance(), hier_componets_, frameId, frameId + 1);

				//pairFrameSimilar.findBestPatches(srLevel - 1, tgLevel - 1, srBestSet,tgBestSet);

				pairFrameSimilar.mergePatchesAfterCoSeg(srLevel - 1, tgLevel - 1,srBestSet,tgBestSet);

			}
		}
	}

}

void DualwayPropagation::buildPatchCorrespondenceByLabel()
{
	for (auto fiter = hier_componets_.begin(); fiter != hier_componets_.end(); ++fiter )
	{
		IndexType frame_id = fiter->first;
		IndexType gLevel = fiter->second.hier_label_bucket.size();
		--gLevel;

		vector<HLabel*>& labelBucket = fiter->second.hier_label_bucket[gLevel];

		if (hier_componets_.find(frame_id + 1) != hier_componets_.end() )
		{

			map<IndexType,IndexType> sr_labelIndexMap = fiter->second.hier_label_vtxBucket_index[gLevel];

			IndexType tLevel = hier_componets_[frame_id + 1].hier_label_vtxBucket_index.size();
			--tLevel;

			map<IndexType,IndexType> tg_labelIndexMap = hier_componets_[frame_id + 1].hier_label_vtxBucket_index[tLevel];

			vector<HLabel*>& tg_labelBucket = hier_componets_[frame_id + 1].hier_label_bucket[tLevel];

			IndexType i = 0;
			for ( auto mIter = sr_labelIndexMap.begin(); mIter != sr_labelIndexMap.end(); ++mIter,++i )
			{
				IndexType labelId = mIter->first;
				if (tg_labelIndexMap.find(labelId) != tg_labelIndexMap.end() )//������ͬ��label
				{
					IndexType corLabelId = tg_labelIndexMap[labelId];

					labelBucket[i]->next_corr = tg_labelBucket[corLabelId]; //next correspondence

					tg_labelBucket[corLabelId]->prev_corr = labelBucket[i]; //prev correspondence
				}
			}

		}
	}
}

void DualwayPropagation::mergeSingleTinyPatches()
{

	GraphMatch pairFrameSimilar(SampleSet::get_instance(), hier_componets_, 8, 9);

	for (auto fiter = hier_componets_.begin(); fiter != hier_componets_.end(); ++ fiter)
	{
		IndexType fId = fiter->first;
		if (hier_componets_.find(fId + 1) != hier_componets_.end())
		{
	             pairFrameSimilar.mergeTinyPatches(fId,0,20);
		}
	}



}


//�Ա�����֮����з���

void DualwayPropagation::split_twoAjacent_graph_next_order(IndexType srFrame, IndexType tgFrame)
{
	//��ǰ����

	//get the new graph of tgGrame--��Ҫ���

	Logger<<" .......\n";
	Logger<<"  Start next split.\n";
	IndexType tgGraSize = hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* oriGra = hier_componets_[tgFrame].hier_graph[tgGraSize - 1];
	LabelsGraph* new_graph = new LabelsGraph(*oriGra);

	//new_graph = oriGra;

	//vector<HLabel* > new_label_bucket =  hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1];
	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1] );

	//
	IndexType gLevel = 0;

	IndexType srGraSize = hier_componets_[srFrame].hier_graph.size();

	gLevel  = srGraSize - 1;//��ȡ���µĲ�
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

	IndexType labParentsize = tgGraSize + 1; //���ɵĲ���

	Logger<<srFrame<<"֡�ĵ�"<<gLevel<<"��߽�ָ�"<<tgFrame<<"��"<<tgGraSize - 1 <<"��"<<endl;

	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*srGraLat);

	generateOrderededges(srFrame,tgFrame);

	while (orderedEdgeQ.size() > 0)
	{
		EdgeSplitOrder oEdge = orderedEdgeQ.top();
		orderedEdgeQ.pop();

		EdgeDescriptor ed = oEdge.EdgeDec;

		GraphEdgeProperty& ep = (*srGraLat)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 1)
		{
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType newGraphEdgeSize = new_graph->m_edges.size();

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��
		//IndexType nodeId = edgeCorrNextVtx.size();


		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		Logger<<"  �ߵ����Ϊ"<<eS<<"�յ�Ϊ"<<eE<<endl;

// 		IndexType recordS = 0;       
// 		IndexType recordE = 0;
// 
// 		IndexType vtxBSzie = splitedLabel->vertex_bucket.size();
// 
//  		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
//  		{
//  			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //
//  
//  			if (prev_id == eS)
//  			{
//  				recordS ++;
//  
//  			}else if(prev_id == eE)
//  			{
//  				recordE ++;
//  			}
//  		}
 
 		//ScalarType ration = (ScalarType)(recordS + recordE)/vtxBSzie;
 
 		//�����ѳ����ĵ������һ�����ݺ���,��ñ߲����ѱ�
 
 		if ( oEdge.srCorNum < 10 || oEdge.tgCorNum < 10 )
 		{
 			Logger<<"�߽��̫����,����Ҫ����.\n";
 			continue;
 		}

		if ( oEdge.unMarkedRation < 0.2)
		{
			Logger<<"Unmark���ֵ̫��,��ʱ������.\n";
			continue;
		}

		//����nodeId �������ӵı�

		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);

		VertexIterator nodeIter = (vi.first + nodeId);

		VertexDescriptor nodeDesc = *nodeIter;

		//�ڵ��Ӧ�����г���
		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(nodeDesc,*new_graph);

		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;

		set<GraphEdgeProperty> collapseEdges;

		OutEdgeIterator oit,nextIt;

		oit = nodeEiter.first;

		for (nextIt = oit; oit != nodeEiter.second; oit = nextIt )
		{
			++nextIt;

			EdgeDescriptor nextEdgeD = *oit;

			GraphEdgeProperty& nextEP = (*new_graph)[nextEdgeD];

			collapseEdges.insert(nextEP);

			boost::remove_edge(*oit,*new_graph);//ɾ��������
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

		// 		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���

		map<IndexType,HVertex*> unMakePs;

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			IndexType vtx_id = iter->first;
			if (vtx_id == 562|| vtx_id == 22176||vtx_id == 913 || vtx_id == 15073|| vtx_id == 15073 || vtx_id == 913)
			{
				Logger<<"stop!.\n";

			}

			IndexType pId = iter->second->prev_corr->vtx_id;

			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //�ò������µ�label_parent��ַ.




			if (prev_id == eS)
			{

				iter->second->label_parent.resize(tgGraSize + 1);
				iter->second->label_parent[tgGraSize] =  new_label_bucket[nodeId];
				++iter;

			}else if(prev_id == eE)
			{

				iter->second->label_parent.resize(tgGraSize + 1);
				iter->second->label_parent[tgGraSize] = new_label_bucket[new_label] ;
				new_label_bucket[new_label]->vertex_bucket.insert(*iter);
				hier_componets_[tgFrame].label_of_vtx[vtx_id ] = new_label;
				iter = splitedLabel->vertex_bucket.erase(iter);


			}else//һЩ��ȷ��label�ĵ�
			{
				unMakePs.insert(*iter);
				iter = splitedLabel->vertex_bucket.erase(iter);
			}

		}

		if (new_label_bucket[new_label]->vertex_bucket.size() == 0 )
		{
			Logger<<"���ѳ����ĵ�̫��.\n";
			new_label_bucket.pop_back();
			continue;
		}

		//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
		determinateUnmarkPoints(tgFrame,unMakePs,new_label_bucket,nodeId,new_label,tgGraSize);


		//������������мӱ߲���,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = nodeId;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;

		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		//��û��index��ֵ

		newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());


		boost::add_edge(nodeId,nSize,newEP,*new_graph);

		//�϶����������ڵ��������ڵ�������߲���ֻ����� recordColapseEdges.size()����.
		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
		{

			GraphEdgeProperty glueEdge;

			glueEdge = (*iter);
			IndexType  sEdgeId = glueEdge.start_;
			IndexType  eEdgeId = glueEdge.end_;

			map<IndexType,HVertex*>  edgePoints;

			auto bIter = glueEdge.edgePoints.begin();
			edgePoints.insert(bIter->second.begin(),bIter->second.end() );

			map<IndexType,HVertex*> startVtx = new_label_bucket[sEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> endVtx = new_label_bucket[eEdgeId]->vertex_bucket;
			map<IndexType,HVertex*> nodeVtx = new_label_bucket[nodeId]->vertex_bucket;
			map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;

			ScalarType minNode = 0.0;
			ScalarType minSize = 0.0;


			if (sEdgeId < nodeId) //nodeidΪ�յ�
			{
				minDistBeTwoParts(tgFrame,startVtx,nodeVtx,minNode);
				minDistBeTwoParts(tgFrame,startVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.end_ = nodeId;
				}else
				{
					glueEdge.end_ = nSize;
				}

			}else//nodeidΪ���
			{
				minDistBeTwoParts(tgFrame,endVtx,nodeVtx,minNode);
				minDistBeTwoParts(tgFrame,endVtx,nSizeVtx,minSize);

				if (minNode < minSize)
				{
					glueEdge.start_ = nodeId;
					glueEdge.end_ = eEdgeId;

				}else
				{
					glueEdge.start_ = eEdgeId;
					glueEdge.end_ = nSize;
				}

			}

			IndexType eKey = frame_index_to_key(glueEdge.start_,glueEdge.end_);

			glueEdge.edgePoints[eKey] = edgePoints;

			boost::add_edge(glueEdge.start_,glueEdge.end_,glueEdge,*new_graph);

		}//����collapse�ı�
	} //end  while


	//0920 �ڱ���target frame�ϵ�ÿ���ڵ�

	//0920 �ڱ���target frame�ϵ�ÿ���ڵ�


	checkPsNewLabelParentPtr(new_label_bucket,labParentsize);//next dirction


	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//�������µ�graph

	hier_componets_[tgFrame].hier_label_vtxBucket_index.push_back(labelIndex);




	Logger<<"  End next split.\n";

	Logger<<" .......\n";



//////////------------------------------------------------------------------------------------------

// 	for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
// 	{
// 		EdgeDescriptor ed = *eit;
// 
// 		GraphEdgeProperty& ep = (*srGraLat)[ed];
// 
// 		map<IndexType,HVertex*> edgePoints;
// 
// 		auto ePsIt = ep.edgePoints.begin();
// 
// 		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  
// 
// 		if (edgePoints.size() < 3)
// 		{
// 			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
// 			continue;
// 		}
// 
// 		map<IndexType,HVertex*> edgeCorrNextVtx;
// 
// 		IndexType newGraphEdgeSize = new_graph->m_edges.size();
// 
// 		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��
// 		//IndexType nodeId = edgeCorrNextVtx.size();
// 
// 
// 		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
// 		HLabel* splitedLabel = new_label_bucket[nodeId];
// 
// 		IndexType eS = ep.start_;
// 		IndexType eE = ep.end_;
// 
// 		Logger<<"  �ߵ����Ϊ"<<eS<<"�յ�Ϊ"<<eE<<endl;
// 
// 		IndexType recordS = 0;       
// 		IndexType recordE = 0;
// 
// 		IndexType vtxBSzie = splitedLabel->vertex_bucket.size();
// 		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
// 		{
// 			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //
// 
// 			if (prev_id == eS)
// 			{
// 				recordS ++;
// 
// 			}else if(prev_id == eE)
// 			{
// 				recordE ++;
// 			}
// 		}
// 
// 		ScalarType ration = (ScalarType)(recordS + recordE)/vtxBSzie;
// 
// 		//�����ѳ����ĵ������һ�����ݺ���,��ñ߲����ѱ�
// 
// 		if ( recordE < 5 || recordS < 5 )
// 		{
// 			Logger<<"�߽��̫����,����Ҫ����.\n";
// 			continue;
// 		}
// 
// 		if (ration < 0.2)
// 		{
// 			Logger<<"Unmark���ֵ̫��,��ʱ������.\n";
// 			continue;
// 		}
// 
// 		//����nodeId �������ӵı�
// 
// 		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);
// 
// 		VertexIterator nodeIter = (vi.first + nodeId);
// 
// 		VertexDescriptor nodeDesc = *nodeIter;
// 
// 		//�ڵ��Ӧ�����г���
// 		pair<OutEdgeIterator,OutEdgeIterator> nodeEiter = boost::out_edges(nodeDesc,*new_graph);
// 
// 		map<IndexType,map<IndexType,HVertex*> > recordColapseEdges;
// 
// 		set<GraphEdgeProperty> collapseEdges;
// 
// 		OutEdgeIterator oit,nextIt;
// 
// 		oit = nodeEiter.first;
// 
// 		for (nextIt = oit; oit != nodeEiter.second; oit = nextIt )
// 		{
// 			++nextIt;
// 
// 			EdgeDescriptor nextEdgeD = *oit;
// 
// 			GraphEdgeProperty& nextEP = (*new_graph)[nextEdgeD];
// 
// 			collapseEdges.insert(nextEP);
// 
// 			boost::remove_edge(*oit,*new_graph);//ɾ��������
// 		}
// 
// 		//����һ���ڵ�
// 
// 		IndexType nSize = boost::num_vertices(*new_graph);
// 
// 		GraphVertexProperty vp(nSize,-1,-1);
// 
// 		boost::add_vertex(vp,*new_graph);
// 
// 
// 		//���·ָ����Ϣ,�����ӵ�Label���ΪnSize. �����ѵĵ�ΪnodeId
// 		IndexType new_label = nSize;
// 
// 		new_label_bucket.push_back((HLabel*)0 );
// 		HLabel* new_label_space = allocator_.allocate<HLabel>();
// 		HLabel* new_label_obj = new (new_label_space)HLabel;
// 		new_label_bucket[new_label] = new_label_obj;
// 		new_label_bucket[new_label]->label_id = new_label;
// 		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];
// 
// 		// 		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
// 
// 		map<IndexType,HVertex*> unMakePs;
// 
// 		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
// 		{
// 			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //�ò������µ�label_parent��ַ.
// 			IndexType vtx_id = iter->first;
// 
// 			if (prev_id == eS)
// 			{
// 
// 				iter->second->label_parent.resize(tgGraSize + 1);
// 				iter->second->label_parent[tgGraSize] =  new_label_bucket[nodeId];
// 				++iter;
// 
// 			}else if(prev_id == eE)
// 			{
// 
// 				iter->second->label_parent.resize(tgGraSize + 1);
// 				iter->second->label_parent[tgGraSize] = new_label_bucket[new_label] ;
// 				new_label_bucket[new_label]->vertex_bucket.insert(*iter);
// 				hier_componets_[tgFrame].label_of_vtx[vtx_id ] = new_label;
// 				iter = splitedLabel->vertex_bucket.erase(iter);
// 
// 
// 			}else//һЩ��ȷ��label�ĵ�
// 			{
// 				unMakePs.insert(*iter);
// 				iter = splitedLabel->vertex_bucket.erase(iter);
// 			}
// 
// 		}
// 
// 		//�����ȡ���������С�������жϲ�ȷ���������ĸ���.unmark Ҫô����nodeid Ҫô����new_label
// 		determinateUnmarkPoints(tgFrame,unMakePs,new_label_bucket,nodeId,new_label,tgGraSize);
// 
// 
// 		//������������мӱ߲���,
// 		//node<--->new_node
// 		GraphEdgeProperty newEP;
// 		newEP.start_ = nodeId;
// 		newEP.end_ = nSize;
// 		newEP.index = newGraphEdgeSize;
// 
// 		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
// 		//��û��index��ֵ
// 
// 		newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());
// 
// 
// 		boost::add_edge(nodeId,nSize,newEP,*new_graph);
// 
// 		//�϶����������ڵ��������ڵ�������߲���ֻ����� recordColapseEdges.size()����.
// 		for (auto iter = collapseEdges.begin(); iter != collapseEdges.end(); iter ++)
// 		{
// 
// 			GraphEdgeProperty glueEdge;
// 
// 			glueEdge = (*iter);
// 			IndexType  sEdgeId = glueEdge.start_;
// 			IndexType  eEdgeId = glueEdge.end_;
// 
// 			map<IndexType,HVertex*>  edgePoints;
// 
// 			auto bIter = glueEdge.edgePoints.begin();
// 			edgePoints.insert(bIter->second.begin(),bIter->second.end() );
// 
// 			map<IndexType,HVertex*> startVtx = new_label_bucket[sEdgeId]->vertex_bucket;
// 			map<IndexType,HVertex*> endVtx = new_label_bucket[eEdgeId]->vertex_bucket;
// 			map<IndexType,HVertex*> nodeVtx = new_label_bucket[nodeId]->vertex_bucket;
// 			map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;
// 
// 			ScalarType minNode = 0.0;
// 			ScalarType minSize = 0.0;
// 
// 
// 			if (sEdgeId < nodeId) //nodeidΪ�յ�
// 			{
// 				minDistBeTwoParts(tgFrame,startVtx,nodeVtx,minNode);
// 				minDistBeTwoParts(tgFrame,startVtx,nSizeVtx,minSize);
// 
// 				if (minNode < minSize)
// 				{
// 					glueEdge.end_ = nodeId;
// 				}else
// 				{
// 					glueEdge.end_ = nSize;
// 				}
// 
// 			}else//nodeidΪ���
// 			{
// 				minDistBeTwoParts(tgFrame,endVtx,nodeVtx,minNode);
// 				minDistBeTwoParts(tgFrame,endVtx,nSizeVtx,minSize);
// 
// 				if (minNode < minSize)
// 				{
// 					glueEdge.start_ = nodeId;
// 					glueEdge.end_ = eEdgeId;
// 
// 				}else
// 				{
// 					glueEdge.start_ = eEdgeId;
// 					glueEdge.end_ = nSize;
// 				}
// 
// 			}
// 
// 			IndexType eKey = frame_index_to_key(glueEdge.start_,glueEdge.end_);
// 
// 			glueEdge.edgePoints[eKey] = edgePoints;
// 
// 			boost::add_edge(glueEdge.start_,glueEdge.end_,glueEdge,*new_graph);
// 
// 		}//����collapse�ı�
// 
// 
// 	}//����ÿ����,ÿ���߶���ʹ��new_graph����һ���µĽڵ�
// 
// 
// 	checkPsNewLabelParentPtr(new_label_bucket,labParentsize);//next dirction
// 
// 	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);
// 
// 	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//�������µ�graph
// 
// 
// 	Logger<<"  End next split.\n";
// 
// 	Logger<<" .......\n";
}

void DualwayPropagation::generateOrderededges(IndexType srFrame, IndexType tgFrame)
{
	Logger<<"�Ա߽�������.\n";

	while (!orderedEdgeQ.empty() )
	{
		orderedEdgeQ.pop();
	}

	//��ǰ����

	//get the new graph of tgGrame--��Ҫ���


	Logger<<" .......\n";
	Logger<<"  Start next split.\n";
	IndexType tgGraSize = hier_componets_[tgFrame].hier_graph.size();
	LabelsGraph* oriGra = hier_componets_[tgFrame].hier_graph[tgGraSize - 1];
	LabelsGraph* new_graph = new LabelsGraph(*oriGra);

	//new_graph = oriGra;

	//vector<HLabel* > new_label_bucket =  hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1];
	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[tgFrame].hier_label_bucket[tgGraSize - 1] );


	//
	IndexType gLevel = 0;

	IndexType srGraSize = hier_componets_[srFrame].hier_graph.size();

	gLevel  = srGraSize - 1;//��ȡ���µĲ�
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

	IndexType labParentsize = tgGraSize + 1; //���ɵĲ���

	//Logger<<srFrame<<"֡�ĵ�"<<gLevel<<"��߽�ָ�"<<tgFrame<<"��"<<tgGraSize - 1 <<"��"<<endl;

	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*srGraLat);


	EdgeSplitOrder tempEdge;

	for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
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

		IndexType newGraphEdgeSize = new_graph->m_edges.size();

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//��ñ߽������һ֡��Ӧ�Ŀ�Ͷ�Ӧ��
		//IndexType nodeId = edgeCorrNextVtx.size();


		//��Ӧ��ȥ,��ǩ����ʼ��,�򱣳ֲ���
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		//Logger<<"  �ߵ����Ϊ"<<eS<<"�յ�Ϊ"<<eE<<endl;

		IndexType recordS = 0;       
		IndexType recordE = 0;

		IndexType vtxBSzie = splitedLabel->vertex_bucket.size();
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
		{
			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //

			if (prev_id == eS)
			{
				recordS ++;

			}else if(prev_id == eE)
			{
				recordE ++;
			}
		}

		ScalarType ration = (ScalarType)(recordS + recordE)/vtxBSzie;


		tempEdge.EdgeDec = ed;
		tempEdge.unMarkedRation = ration;
		tempEdge.srCorNum = recordS;
		tempEdge.tgCorNum = recordE;


		orderedEdgeQ.push(tempEdge);
	}

	Logger<<"����������.\n";
	//��Ҫɾ��������ڴ�!!!
}


void DualwayPropagation::generateOrderPrevEdges(IndexType srFrame, IndexType tgFrame)
{
	while (!orderedEdgeQ.empty() )
	{
		orderedEdgeQ.pop();
	}

	//��ȡ��Ҫ���µ�graph
	Logger<<" .......\n";
	Logger<<"  Start prev split.\n";

	IndexType srGraphSize = hier_componets_[srFrame].hier_graph.size();
	IndexType tgGraphSize = hier_componets_[tgFrame].hier_graph.size();

	assert(srGraphSize > 0 && tgGraphSize > 0);

	LabelsGraph* oriSpGra = hier_componets_[srFrame].hier_graph[srGraphSize - 1];
	LabelsGraph* shouldSplitGraph = new LabelsGraph(*oriSpGra);

	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[srFrame].hier_label_bucket[srGraphSize - 1] );

	//��ȡָ���ָ��graph
	LabelsGraph* guideSplitGraph;
	IndexType graphLevel = 0;

	assert(tgGraphSize > 1);

	if (tgGraphSize > 2)
	{
		graphLevel = tgGraphSize - 1;
	}else
	{
		graphLevel = tgGraphSize - 2;
	}
	guideSplitGraph = hier_componets_[tgFrame].hier_graph[graphLevel];

	//�Ա߽������ȼ�����

	// guideSplitGraph��ÿ��������һ�ηָ�
	pair<EdgeIterator,EdgeIterator> ei = boost::edges(*guideSplitGraph);
	EdgeSplitOrder tempEdge;

	for (EdgeIterator eit = ei.first; eit != ei.second; ++eit)
	{
		EdgeDescriptor ed = *eit;

		GraphEdgeProperty& ep = (*guideSplitGraph)[ed];

		map<IndexType,HVertex*> edgePoints;

		auto ePsIt = ep.edgePoints.begin();

		edgePoints.insert(ePsIt->second.begin(),ePsIt->second.end() );  

		if (edgePoints.size() < 1)
		{
			Logger<<"���ϵĶ�����̫��,�޷�����.\n";
			continue;
		}

		Logger<<tgFrame<<"֡�ı����"<<ep.start_<<"���յ�"<<ep.end_<<"��ʼ����"<<endl;

		map<IndexType,HVertex*> edgeCorrPrevVtx;

		IndexType nSize = boost::num_vertices(*shouldSplitGraph);

		IndexType newGraphEdgeSize = shouldSplitGraph->m_edges.size(); //Ϊ�˸������ӵı�������

		bool isSplit = true;

		IndexType edgePsCorNode = 0;

		isSplit = checkPrevLabelBucket(edgePoints,edgeCorrPrevVtx,edgePsCorNode);

		if ( edgePsCorNode < 0 || edgePsCorNode > nSize - 1)
		{
			Logger<<"�߽��ҵ��Ŀ�Խ��!.\n";
			break;
		}

		if (!isSplit)
		{
			Logger<<"�߽��ָ����������,���ָ�.\n";
			continue;
		}

		HLabel* splitedLabel = new_label_bucket[edgePsCorNode]; //������Ҫ���ѵĿ�

		IndexType curEdgeStart = ep.start_;
		IndexType curEdgeEnd   = ep.end_;
		IndexType strCorPsSzie = 0;
		IndexType endCorPsSize = 0;

		//����ǰ��һ���򵥵�Ԥ�ж�
		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); iter ++)
		{
			if (iter->second->next_corr->label_parent[graphLevel] != NULL)
			{

				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; 

				if (nextVtx_label == curEdgeStart)
				{
					strCorPsSzie ++;

				}else if(nextVtx_label == curEdgeEnd)
				{
					endCorPsSize ++;
				}
			}

		}

		IndexType  vtxBSize = splitedLabel->vertex_bucket.size();

		ScalarType ration = (ScalarType)(strCorPsSzie + endCorPsSize)/vtxBSize;

		tempEdge.EdgeDec = ed;
		tempEdge.unMarkedRation = ration;
		tempEdge.srCorNum = strCorPsSzie;
		tempEdge.tgCorNum = endCorPsSize;


		orderedEdgeQ.push(tempEdge);

	}

}


void DualwayPropagation::constructPCloudGraph()
{
	Logger<<"  Begin initialize point cloud graphs.\n";

	for (auto citer = hier_componets_.begin(); citer!=hier_componets_.end(); citer++)
	{
		
		IndexType nodeSize = citer->second.label_of_vtx.size();//.size();

		PCloudGraph* new_pcGraph_space = allocator_.allocate<PCloudGraph>();

		PCloudGraph* new_pcGraph = new (new_pcGraph_space)PCloudGraph;

		addGraphVertex(*new_pcGraph,citer->first);

		addGraphEdge(*new_pcGraph,citer->first);

		citer->second.pcGraph = new_pcGraph;
	    
	}

   Logger<<"  End initialize point cloud graphs.\n";
}

void DualwayPropagation::addGraphVertex(PCloudGraph& pcGraph, IndexType frameId)
{
	IndexType nSize = hier_componets_[frameId].label_of_vtx.size();

	auto vIter = hier_componets_[frameId].label_of_vtx.begin();
	auto vEnd = hier_componets_[frameId].label_of_vtx.end();

	IndexType gIndex = 0;
	map<IndexType,IndexType> gNodeIdx;

	for (; vIter != vEnd; ++ vIter,++ gIndex)
	{
		IndexType labelId = (*vIter).second;
		IndexType vtxId = (*vIter).first;

		IndexType idxm = hier_componets_[frameId].hier_label_vtxBucket_index[0][labelId];
		HVertex* vtx = hier_componets_[frameId].hier_label_bucket[0][idxm]->vertex_bucket[vtxId];

		gNodeIdx[vtxId] = gIndex;

		PCVertexProperty vp;
		vp.index = gIndex;
		vp.vtxSite = vtx;
		boost::add_vertex(vp,pcGraph);	
	}

	hier_componets_[frameId].gId_of_vtx = gNodeIdx;

}
void DualwayPropagation::addGraphEdge(PCloudGraph& pcGraph, IndexType frameId)
{
	IndexType nSize = hier_componets_[frameId].label_of_vtx.size();

	auto vIter = hier_componets_[frameId].label_of_vtx.begin();
	auto vEnd = hier_componets_[frameId].label_of_vtx.end();

	unordered_map<IndexType,bool> recordEdges;

	buildKdTree(frameId);

	IndexType gIndex = 0;
	const IndexType k = 6;
	IndexType neighbours[k];
	ScalarType dist[k];

	IndexType edNum = 0;
	for (; vIter != vEnd; ++ vIter,++ gIndex)
	{
	    IndexType vtxId = (*vIter).first;
		downSample->neighbours(gIndex,k,neighbours,dist);

		for (IndexType i = 1; i < k; ++i)
		{
		  bool temp = recordEdges[frame_index_to_key(gIndex,neighbours[i]) ];

		  if (!temp)
		  {
			  PCEdgeProperty ep;
			  ep.index = edNum;
			  ep.start_ = gIndex;
			  ep.end_ = neighbours[i];
			  ep.dist = dist[i];

			  boost::add_edge(gIndex,neighbours[i],ep,pcGraph);

			  recordEdges[frame_index_to_key(neighbours[i],gIndex)] = true;

			  ++ edNum;
		  }
		}
	}

}

void DualwayPropagation::show_corresponding(int f)
{


	for ( IndexType l = 0; l<hier_componets_[f].hier_label_bucket[0].size(); l++ )
	{
		HLabel& label = *hier_componets_[f].hier_label_bucket[0][l];
		IndexType mSize = label.vertex_bucket.size();
		IndexType i = 5;
		for ( auto viter = label.vertex_bucket.begin();
			viter!=label.vertex_bucket.end() && i < mSize;
			i += 5, advance(viter,5))
		{
			HVertex& vtx = *(viter->second);
			if ( vtx.next_corr )
			{
				Tracer::get_instance().add_record(f, vtx.vtx_id, f+1, vtx.next_corr->vtx_id);
			}

// 			if (vtx.prev_corr)
// 			{
// 				Tracer::get_instance().add_record(f, vtx.vtx_id, f-1, vtx.prev_corr->vtx_id);
// 			}

		}
	}

	//old version
// 	for ( IndexType l = 0; l<components_[f].label_bucket.size(); l++ )
// 	{
// 		CLabel& label = *components_[f].label_bucket[l];
// 		for ( auto viter = label.vertex_bucket.begin();
// 			viter!=label.vertex_bucket.end();
// 			viter++)
// 		{
// 			CVertex& vtx = *(viter->second);
// 			if ( vtx.next_corr )
// 			{
// 				Tracer::get_instance().add_record(f, vtx.vtx_id, f+1, vtx.next_corr->vtx_id);
// 			}
// 			if (vtx.prev_corr)
// 			{
// 				Tracer::get_instance().add_record(f, vtx.vtx_id, f-1, vtx.prev_corr->vtx_id);
// 			}
// 		}
// 	}
}

