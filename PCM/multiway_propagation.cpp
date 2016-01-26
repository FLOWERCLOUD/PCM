#include"multiway_propagation.h"


#define  INF_LOCAL 1000000
#define  COUT_DEBUG 0

#define frame_index_to_key(f,i) ((f<<16)|i)
#define frame_label_to_key(f,l) ((f<<16)|l)
#define get_index_from_key(k) (k&0xffff)
#define get_frame_from_key(k) (k>>16)
#define get_current_dir _getcwd


PoolAllocator DualwayPropagation::allocator_;

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


	//删掉label_buckert中为空的元素
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
			
		} ;//它就是边界点
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


void DualwayPropagation:: init_labeles_graph_hier(ScalarType distThr)
{

	Logger<<"  Begin initialize graphs.\n";

	if (distThr >= 1.0 || distThr < 0.0)
	{
		distThr = 0.05;
	}

	for (auto citer = hier_componets_.begin(); citer!=hier_componets_.end(); citer++)
	{

		IndexType lbsize = citer->second.hier_label_bucket.size();
		assert( lbsize > 0);
		IndexType nodeSize = citer->second.hier_label_bucket[lbsize -1].size();//.size();

		LabelsGraph* new_labelGraph_space = allocator_.allocate<LabelsGraph>();

		LabelsGraph* new_labelGraph = new (new_labelGraph_space)LabelsGraph;

		GraphVertexProperty gvp;

		map<IndexType,IndexType> labelLevel = citer->second.hier_label_vtxBucket_index[lbsize - 1];

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

				//计算两个块之间的最短距离
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

						//用dijstra最短距离来计算两点之间的距离

						if(distance < minDis)minDis = distance;
						IndexType label1 = hier_componets_[citer->first].label_of_vtx[ index1];
						IndexType label2 = hier_componets_[citer->first].label_of_vtx[ index2];
						if( distance < dia*distThr)
						{
							pointdistance tpd( index1 ,  label1 , index2 ,label2 ,distance );
							PriQue.push( tpd);
						}

					}

				}

				if (minDis <  dia * distThr )//0.0588
				{
					getEdgeVertexs2( citer->first ,PriQue ,gep.edgePoints);

					//边界点集个数太少也不需要添加,可能是噪声点.

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
						// i == j 时不处理
					}

				}

			}

		}

		//citer->second.hier_graph.clear();

		citer->second.hier_graph.push_back( new_labelGraph);

	}//帧遍历结束

	Logger<<"  End initialize graphs.\n";
}

void DualwayPropagation::split_twoAjacent_graph_next(IndexType srFrame, IndexType tgFrame)
{
	//向前分裂

	//get the new graph of tgGrame--需要深拷贝


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

	gLevel  = srGraSize - 1;//获取最新的层
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

    IndexType labParentsize = tgGraSize + 1; //生成的层数

	Logger<<srFrame<<"帧的第"<<gLevel<<"层边界分割"<<tgFrame<<"的"<<tgGraSize - 1 <<"层"<<endl;

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
			Logger<<"边上的顶点数太少,无法分裂.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType newGraphEdgeSize = new_graph->m_edges.size();

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//获得边界点在下一帧对应的块和对应点
		//IndexType nodeId = edgeCorrNextVtx.size();


		//对应回去,标签是起始点,则保持不变
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		Logger<<"  边的起点为"<<eS<<"终点为"<<eE<<endl;

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

		//若分裂出来的点个数有一个数据很少,则该边不做裂变

		if ( recordE < 5 || recordS < 5 )
		{
			Logger<<"边界点太靠近,不需要分裂.\n";
			continue;
		}

		if (ration < 0.2)
		{
			Logger<<"Unmark点比值太大,暂时不分裂.\n";
			continue;
		}

		//遍历nodeId 上相连接的边

		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);

		VertexIterator nodeIter = (vi.first + nodeId);

		VertexDescriptor nodeDesc = *nodeIter;

		//节点对应的所有出边
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
			
			boost::remove_edge(*oit,*new_graph);//删除这条边
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
 
// 		//对应回去,标签是起始点,则保持不变

		 map<IndexType,HVertex*> unMakePs;

		// Sample& smp = SampleSet::get_instance()[tgFrame];
		 //ScalarType egThreshold = smp.getBox().diag();

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //得不到最新的label_parent地址.
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


 			}else//一些待确定label的点
 			{
                   unMakePs.insert(*iter);
                   iter = splitedLabel->vertex_bucket.erase(iter);
 			}

		}

		//用随机取点产生的最小距离来判断不确定点属于哪个类.unmark 要么属于nodeid 要么属于new_label
		determinateUnmarkPoints(tgFrame,unMakePs,new_label_bucket,nodeId,new_label,tgGraSize);


		//对这两个点进行加边操作,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = nodeId;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;

		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		//还没对index赋值

		newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());


		boost::add_edge(nodeId,nSize,newEP,*new_graph);

		//断定查找两个节点与其它节点进行连边操作只会出现 recordColapseEdges.size()次数.
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


			if (sEdgeId < nodeId) //nodeid为终点
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
							
			}else//nodeid为起点
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

		}//遍历collapse的边


	}//遍历每条边,每条边都会使得new_graph增加一个新的节点


	

	checkPsNewLabelParentPtr(new_label_bucket,labParentsize);//next dirction


	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//保存最新的graph

	hier_componets_[tgFrame].hier_label_vtxBucket_index.push_back(labelIndex);

	Logger<<"  End next split.\n";
	Logger<<" .......\n";
}

void DualwayPropagation::split_twoAjacent_graph_prev(IndexType srFrame, IndexType tgFrame)
{

	//获取需要更新的graph
	Logger<<" .......\n";
	Logger<<"  Start prev split.\n";

	IndexType srGraphSize = hier_componets_[srFrame].hier_graph.size();
	IndexType tgGraphSize = hier_componets_[tgFrame].hier_graph.size();

	assert(srGraphSize > 0 && tgGraphSize > 0);

	LabelsGraph* oriSpGra = hier_componets_[srFrame].hier_graph[srGraphSize - 1];
	LabelsGraph* shouldSplitGraph = new LabelsGraph(*oriSpGra);

	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[srFrame].hier_label_bucket[srGraphSize - 1] );

	//获取指导分割的graph
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

	//对边进行优先级排序




	// guideSplitGraph的每条边引导一次分割
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
			Logger<<"边上的顶点数太少,无法分裂.\n";
			continue;
		}

		Logger<<tgFrame<<"帧的边起点"<<ep.start_<<"与终点"<<ep.end_<<"开始分裂"<<endl;

		map<IndexType,HVertex*> edgeCorrPrevVtx;

		IndexType nSize = boost::num_vertices(*shouldSplitGraph);

		IndexType newGraphEdgeSize = shouldSplitGraph->m_edges.size(); //为了给新增加的边添加序号

		bool isSplit = true;

		IndexType edgePsCorNode = 0;

	    isSplit = checkPrevLabelBucket(edgePoints,edgeCorrPrevVtx,edgePsCorNode);

		if ( edgePsCorNode < 0 || edgePsCorNode > nSize - 1)
		{
			Logger<<"边界找到的块越界!.\n";
			break;
		}

		if (!isSplit)
		{
			Logger<<"边界点指向了两个块,不分割.\n";
			continue;
		}

		HLabel* splitedLabel = new_label_bucket[edgePsCorNode]; //可能需要分裂的块

		IndexType curEdgeStart = ep.start_;
		IndexType curEdgeEnd   = ep.end_;
		IndexType strCorPsSzie = 0;
		IndexType endCorPsSize = 0;

		//分类前做一个简单的预判断
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
			Logger<<"边界点靠近边界(对应点个数),不需要分裂.\n";
			continue;
		}

		if ( ration < 0.2)
		{
			Logger<<"四两拨千斤?算了吧!.\n";
			continue;
		}

		//用eit 边开始分裂
		pair<VertexIterator, VertexIterator> vi = boost::vertices(*shouldSplitGraph);

		VertexIterator curNodeIter = (vi.first + edgePsCorNode);

		VertexDescriptor curNodeDesc = *curNodeIter;

		//节点对应的所有出边
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

			boost::remove_edge(*oit,*shouldSplitGraph);//删除这条边
		}

		//增加一个节点

		GraphVertexProperty vp(nSize,-1,-1);
		boost::add_vertex(vp,*shouldSplitGraph);

		//更新分割块信息,新增加的Label标号为nSize. 被分裂的点为edgePsCorNode
		IndexType new_label = nSize;
		new_label_bucket.push_back((HLabel*)0 );
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label_obj = new (new_label_space)HLabel;
		new_label_bucket[new_label] = new_label_obj;
		new_label_bucket[new_label]->label_id = new_label;
		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];

		//开始分裂
		map<IndexType,HVertex*> unMarkPs;

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			if (iter->second->next_corr->label_parent[graphLevel])
			{
				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; //得不到最新的label_parent地址.

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

				}else//一些待确定label的点
				{
					unMarkPs.insert(*iter);
					iter = splitedLabel->vertex_bucket.erase(iter);
				}
			}

		} //遍历需要分裂块的每个点

		if (!unMarkPs.empty())
		{
			if ( (!new_label_bucket[new_label]->vertex_bucket.empty() ) && (!new_label_bucket[edgePsCorNode]->vertex_bucket.empty() ) )
			{
				//用随机取点产生的最小距离来判断不确定点属于哪个类.unmark 要么属于nodeid 要么属于new_label
				determinateUnmarkPoints(srFrame,unMarkPs,new_label_bucket,edgePsCorNode,new_label,srGraphSize);

			}else if (new_label_bucket[new_label]->vertex_bucket.empty())
			{
				Logger<<"分裂出来的点太少.\n";

				new_label_bucket.pop_back();

				for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)//放回原来的块中
				{
					new_label_bucket[edgePsCorNode]->vertex_bucket.insert(*iter);

					hier_componets_[tgFrame].label_of_vtx[iter->first] = edgePsCorNode;

					iter->second->label_parent.resize(srGraphSize + 1);
					iter->second->label_parent[srGraphSize] = new_label_bucket[edgePsCorNode];

					iter = unMarkPs.erase(iter); 

				}

				continue;

			}else if(new_label_bucket[edgePsCorNode]->vertex_bucket.empty() )
			{
				Logger<<"ori parch empty!.\n";

				break;
			}

		}

		//用随机取点产生的最小距离来判断不确定点属于哪个类.unmark 要么属于nodeid 要么属于new_label
		//determinateUnmarkPoints(srFrame,unMakePs,new_label_bucket,edgePsCorNode,new_label,srGraphSize);


		//对这两个点进行加边操作,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = edgePsCorNode;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;
		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		newEP.edgePoints[edgeKey].insert(edgeCorrPrevVtx.begin(),edgeCorrPrevVtx.end());
		boost::add_edge(edgePsCorNode,nSize,newEP,*shouldSplitGraph);


		//断定查找两个节点与其它节点进行连边操作只会出现 recordColapseEdges.size()次数.
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


			if (sEdgeId < edgePsCorNode) //nodeid为终点
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

			}else//nodeid为起点
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

		}//遍历collapse的边

	} //遍历引导分割图的每条边

	checkPsNewLabelParentPtr(new_label_bucket,srGraphSize + 1);

	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[srFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[srFrame].hier_graph.push_back(shouldSplitGraph);//保存最新的graph

	hier_componets_[srFrame].hier_label_vtxBucket_index.push_back(labelIndex);

    Logger<<"  End prev split.\n";
	Logger<<" .......\n";
}

void DualwayPropagation::split_twoAjacent_graph_prev_order(IndexType srFrame, IndexType tgFrame)
{
// 	//获取需要更新的graph
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
// 	//获取指导分割的graph
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
// 	// guideSplitGraph的每条边引导一次分割
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
// 			Logger<<"边上的顶点数太少,无法分裂.\n";
// 			continue;
// 		}
// 
// 		map<IndexType,HVertex*> edgeCorrNextVtx;
// 
// 		IndexType newGraphEdgeSize = shouldSplitGraph->m_edges.size();
// 
// 		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//获得边界点在下一帧对应的块和对应点
// 		//IndexType nodeId = edgeCorrNextVtx.size();
// 
// 
// 		//对应回去,标签是起始点,则保持不变
// 		HLabel* splitedLabel = new_label_bucket[nodeId];
// 
// 		IndexType eS = ep.start_;
// 		IndexType eE = ep.end_;
// 
// 		Logger<<"  边的起点为"<<eS<<"终点为"<<eE<<endl;
// 
// 		if ( oEdge.srCorNum < 10 || oEdge.tgCorNum < 10 )
// 		{
// 			Logger<<"边界点太靠近,不需要分裂.\n";
// 			continue;
// 		}
// 
// 		if ( oEdge.unMarkedRation < 0.2)
// 		{
// 			Logger<<"Unmark点比值太大,暂时不分裂.\n";
// 			continue;
// 		}
// 
// 		//用eit 边开始分裂
// 		pair<VertexIterator, VertexIterator> vi = boost::vertices(*shouldSplitGraph);
// 
// 		VertexIterator curNodeIter = (vi.first + edgePsCorNode);
// 
// 		VertexDescriptor curNodeDesc = *curNodeIter;
// 
// 		//节点对应的所有出边
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
// 			boost::remove_edge(*oit,*shouldSplitGraph);//删除这条边
// 		}
// 
// 		//增加一个节点
// 
// 		GraphVertexProperty vp(newGraphEdgeSize,-1,-1);
// 		boost::add_vertex(vp,*shouldSplitGraph);
// 
// 		//更新分割块信息,新增加的Label标号为nSize. 被分裂的点为edgePsCorNode
// 		IndexType new_label = newGraphEdgeSize;
// 		new_label_bucket.push_back((HLabel*)0 );
// 		HLabel* new_label_space = allocator_.allocate<HLabel>();
// 		HLabel* new_label_obj = new (new_label_space)HLabel;
// 		new_label_bucket[new_label] = new_label_obj;
// 		new_label_bucket[new_label]->label_id = new_label;
// 		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];
// 
// 		//开始分裂
// 		map<IndexType,HVertex*> unMakePs;
// 
// 		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
// 		{
// 			if (iter->second->next_corr->label_parent[graphLevel])
// 			{
// 				IndexType nextVtx_label = iter->second->next_corr->label_parent[graphLevel]->label_id; //得不到最新的label_parent地址.
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
// 				}else//一些待确定label的点
// 				{
// 					unMakePs.insert(*iter);
// 					iter = splitedLabel->vertex_bucket.erase(iter);
// 				}
// 			}
// 
// 		} //遍历需要分裂块的每个点
// 
// 		//用随机取点产生的最小距离来判断不确定点属于哪个类.unmark 要么属于nodeid 要么属于new_label
// 		determinateUnmarkPoints(srFrame,unMakePs,new_label_bucket,edgePsCorNode,new_label,srGraphSize);
// 
// 		//对这两个点进行加边操作,
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
// 		//断定查找两个节点与其它节点进行连边操作只会出现 recordColapseEdges.size()次数.
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
// 			if (sEdgeId < edgePsCorNode) //nodeid为终点
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
// 			}else//nodeid为起点
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


// 	IndexType startF = 4;
// 	while (startF -- > 0)
// 	{
// 		++cIter;
// 	}

	--cEnd;

	IndexType iterNum = 0;

	IndexType startFrameId = cIter->first;



	for (; cIter != cEnd && iterNum < iterN; ++ cIter++, ++iterNum)
	{
		//最后一个节点不可以运算
		IndexType srFrame = cIter->first;
		IndexType tgFrame = srFrame + 1;

		//split_twoAjacent_graph_next(srFrame,tgFrame );

		//show_corresponding(srFrame);

		split_twoAjacent_graph_next_order(srFrame,tgFrame );

		split_nest_graph_prev(startFrameId,srFrame,tgFrame);

	}

	//show_corresponding(1);

	//输出每帧的分割块数

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); ++ fIter)
	{
		IndexType gLevel = fIter->second.hier_label_bucket.size();
		IndexType labeSize = fIter->second.hier_label_bucket[gLevel - 1].size();
		Logger<<"  第"<<fIter->first<<"帧共分割出"<<labeSize<<"个块.\n";
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

	//	Logger<<fId<<"帧到"<<fId + 1<<"帧的对应.\n";

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
		Logger<<"边界点指向了下个块的两个Label!.\n";

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
		Logger<<"边界点指向了下个块的两个Label!.\n";

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
	//read_label_file(label_name); //最原始的函数

	//read_label_file_hier(label_name); //读取层次信息

	//read_label_file_coseg(label_name); //读取共分割的label文件


    //201501010

  	read_label_file_coseg(label_name); //读取共分割的label文件,处理帧之间不连续label情形

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
		IndexType gLevel = fIter->second.hier_label_bucket.size();//访问最高层的label_bucket
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

				//计算两个块之间的最短距离
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
		Logger<<"没有Unmark 点.\n";
		return;
	}

	Sample& smp = SampleSet::get_instance()[cFrame];

	for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)
	{
		PointType pCoor = smp.vertices_matrix().col(iter->first);

// 		ScalarType d2node = p2PatchAvgDis(cFrame,pCoor,oriLabelBucket[nodeId]->vertex_bucket);//随机取值的最小值
// 		ScalarType d2size = p2PatchAvgDis(cFrame,pCoor,oriLabelBucket[newLabe]->vertex_bucket );


// 		ScalarType d2node = p2PatchMinDis(cFrame,pCoor,oriLabelBucket[nodeId]->vertex_bucket);//距离严格最小
// 		ScalarType d2size = p2PatchMinDis(cFrame,pCoor,oriLabelBucket[newLabe]->vertex_bucket );


// 		//测地线距离
// 		ScalarType d2node = 0.;
// 		ScalarType d2size = 0.;
// 
// 		if (oriLabelBucket[nodeId]->vertex_bucket.empty() )
// 		{
// 			d2node = 1e5;
// 		}else
// 		{
// 			d2node = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[nodeId]->vertex_bucket);//距离严格最小
// 		}
// 
// 		if (oriLabelBucket[newLabe]->vertex_bucket.empty() )
// 		{
// 			d2size = 1e5;
// 		}else
// 		{
// 			d2size = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[newLabe]->vertex_bucket );
// 		}

		ScalarType d2node = p2PatchGeoDis(cFrame,*(*iter).second,oriLabelBucket[nodeId]->vertex_bucket);//距离严格最小
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

			//avgDis += (pCoor - tgCoor).norm(); //平均值不合理

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

	//顶点描述符
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

// 确保每个顶点的labparent指针数目层数相同
void DualwayPropagation::checkPsNewLabelParentPtr(vector<HLabel*> oriLabelBucket,IndexType labParSize)
{
	IndexType vId = 0;
	for (auto vIter = oriLabelBucket.begin(); vIter != oriLabelBucket.end(); vIter ++,vId ++)
	{
		 map<IndexType,HVertex*>& vtxBucket = (*vIter)->vertex_bucket;

		 if (vtxBucket.empty())
		 {
			 continue;
		 }

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
void DualwayPropagation::smoothSingleFrame(IndexType frame_idx)//没有修改根据label号码来索引顶点集
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

	//只处理最后一层的label_bucket数据
	map<IndexType, IndexType> old_label_map = hier_componets_[frame_idx].label_of_vtx;

	IndexType graphLevel = hier_componets_[frame_idx].hier_label_vtxBucket_index.size();
	IndexType hLevel_idx = graphLevel - 1;
	//IndexType beforLabelSize = hier_componets_[frame_idx].hier_label_bucket[hLevel_idx].size();

	IndexType beforLabelSize = 100;

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
	map<IndexType,IndexType> patchTrajId;//给相同的patch赋相同的id
	map<IndexType,bool>      isPatchTrav;//记录块是否访问过--key = frameId  & labelId
	
	IndexType totLabelId = 0;
	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
	{
		IndexType gLevel = fIter->second.hier_label_bucket.size();//访问最高层的label_bucket
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

	//输出全局的Label
	writeSquenceLabel(patchTrajId,filename);

}

void DualwayPropagation::writeSquenceLabel(map<IndexType,IndexType>& patchTrajId,std::string filename)
{
		FILE* outfile = fopen(filename.c_str(),"w");

		for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
		{
			IndexType gLevel = fIter->second.hier_label_bucket.size();//访问最高层的label_bucket
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
	
    //给定merge终止条件 ,读取Co-segmentation的label文件之后操作该函数

	buildPatchCorrespondenceByLabel();

	IndexType itNum = 9;
	IndexType i = 0;

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end() &&  i < itNum ; ++ fIter,++i)
	{
		IndexType frameId = fIter->first;
		if (hier_componets_.find(frameId + 1) != hier_componets_.end() )
		{
			IndexType combineNum = 1;//合并次数

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


			for ( auto mIter = sr_labelIndexMap.begin(); mIter != sr_labelIndexMap.end(); ++mIter )
			{
				IndexType labelId = mIter->first;
				IndexType realSrLabIdx = mIter->second;

				if (tg_labelIndexMap.find(labelId) != tg_labelIndexMap.end() )//存在相同的label
				{
					IndexType corLabelId = tg_labelIndexMap[labelId];

					labelBucket[realSrLabIdx]->next_corr = tg_labelBucket[corLabelId]; //next correspondence

					tg_labelBucket[corLabelId]->prev_corr = labelBucket[realSrLabIdx]; //prev correspondence
				}
			}

		}
	}
}
 //---------------------------------------

void DualwayPropagation::mergeSingleTinyPatches(IndexType vSize)
{

	GraphMatch pairFrameSimilar(SampleSet::get_instance(), hier_componets_, 8, 9);

	for (auto fiter = hier_componets_.begin(); fiter != hier_componets_.end(); ++ fiter)
	{
		IndexType fId = fiter->first;
		if (hier_componets_.find(fId + 1) != hier_componets_.end())
		{
	        pairFrameSimilar.mergeTinyPatches(fId,vSize);
		}
	}

	Logger<<"After merge tiny patches.\n";

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); ++ fIter)
	{
		IndexType gLevel = fIter->second.hier_label_bucket.size();
		IndexType labeSize = fIter->second.hier_label_bucket[gLevel - 1].size();
		Logger<<"  第"<<fIter->first<<"帧共分割出"<<labeSize<<"个块.\n";
	}

}


//对边排序之后进行分裂

void DualwayPropagation::split_twoAjacent_graph_next_order(IndexType srFrame, IndexType tgFrame)
{
	//向前分裂

	//get the new graph of tgGrame--需要深拷贝

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

	gLevel  = srGraSize - 1;//获取最新的层
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

	IndexType labParentsize = tgGraSize + 1; //生成的层数

	Logger<<srFrame<<"帧的第"<<gLevel<<"层边界分割"<<tgFrame<<"的"<<tgGraSize - 1 <<"层"<<endl;

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
			Logger<<"边上的顶点数太少,无法分裂.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType newGraphEdgeSize = new_graph->m_edges.size();

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//获得边界点在下一帧对应的块和对应点
		//IndexType nodeId = edgeCorrNextVtx.size();


		//对应回去,标签是起始点,则保持不变
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		Logger<<"  边的起点为"<<eS<<"终点为"<<eE<<endl;

		
 		//若分裂出来的点个数有一个数据很少,则该边不做裂变
 
 		if ( oEdge.srCorNum < 10 || oEdge.tgCorNum < 10 )
 		{
 			Logger<<"边界点太靠近,不需要分裂.\n";
 			continue;
 		}

		if ( oEdge.unMarkedRation < 0.2)
		{
			Logger<<"Unmark点比值太大,暂时不分裂.\n";
			continue;
		}

		//遍历nodeId 上相连接的边

		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);

		VertexIterator nodeIter = (vi.first + nodeId);

		VertexDescriptor nodeDesc = *nodeIter;

		//节点对应的所有出边
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

			boost::remove_edge(*oit,*new_graph);//删除这条边
		}

// 		//增加一个节点// 增加顶点操作放在前面.
// 
 		IndexType nSize = boost::num_vertices(*new_graph);
// 
// 		GraphVertexProperty vp(nSize,-1,-1);
// 
// 		boost::add_vertex(vp,*new_graph);

		//更新分割块信息,新增加的Label标号为nSize. 被分裂的点为nodeId
		IndexType new_label = nSize;

		new_label_bucket.push_back((HLabel*)0 );
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label_obj = new (new_label_space)HLabel;
		new_label_bucket[new_label] = new_label_obj;
		new_label_bucket[new_label]->label_id = new_label;
		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];

		// 		//对应回去,标签是起始点,则保持不变

		map<IndexType,HVertex*> unMarkPs;

		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
		{
			IndexType vtx_id = iter->first;

			IndexType pId = iter->second->prev_corr->vtx_id;

			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //得不到最新的label_parent地址.

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


			}else//一些待确定label的点
			{
				unMarkPs.insert(*iter);
				iter = splitedLabel->vertex_bucket.erase(iter);
			}

		}

		if (!unMarkPs.empty())
		{
			if ( (!new_label_bucket[new_label]->vertex_bucket.empty() ) && (!new_label_bucket[nodeId]->vertex_bucket.empty() ) )
			{
				//用随机取点产生的最小距离来判断不确定点属于哪个类.unmark 要么属于nodeid 要么属于new_label
				determinateUnmarkPoints(tgFrame,unMarkPs,new_label_bucket,nodeId,new_label,tgGraSize);

			}else if (new_label_bucket[new_label]->vertex_bucket.empty())
			{
				Logger<<"!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!分裂出来的点太少.\n";

				new_label_bucket.pop_back();

				for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)//放回原来的块中--并没增加节点
				{
					new_label_bucket[nodeId]->vertex_bucket.insert(*iter);

					hier_componets_[tgFrame].label_of_vtx[iter->first] = nodeId;

					iter->second->label_parent.resize(tgGraSize + 1);
					iter->second->label_parent[tgGraSize] = new_label_bucket[nodeId];

					iter = unMarkPs.erase(iter); 

				}

				continue;

			}else if(new_label_bucket[nodeId]->vertex_bucket.empty() )
			{
				Logger<<"ori parch empty!.\n";
				//unmarked和new_label都放回到原来的块中.

				for (auto iter = unMarkPs.begin(); iter != unMarkPs.end();)//放回原来的块中--并没增加节点
				{
					new_label_bucket[nodeId]->vertex_bucket.insert(*iter);

					hier_componets_[tgFrame].label_of_vtx[iter->first] = nodeId;

					iter->second->label_parent.resize(tgGraSize + 1);
					iter->second->label_parent[tgGraSize] = new_label_bucket[nodeId];

					iter = unMarkPs.erase(iter); 

				}

				//把新生成的块也放回原块中

				auto bMapVtx = new_label_bucket[new_label]->vertex_bucket.begin();
				auto eMapVtx = new_label_bucket[new_label]->vertex_bucket.end();

				for (; bMapVtx != eMapVtx; ++ bMapVtx)
				{
					new_label_bucket[nodeId]->vertex_bucket.insert(*bMapVtx);

					hier_componets_[tgFrame].label_of_vtx[bMapVtx->first] = nodeId;

					bMapVtx->second->label_parent.resize(tgGraSize + 1);

					bMapVtx->second->label_parent[tgGraSize] = new_label_bucket[nodeId];

				}

				new_label_bucket.pop_back();

				break;
			}

		}


		//增加一个节点 --把增加顶点操作放在确定增加一个类上面.

		//IndexType nSize = boost::num_vertices(*new_graph);

		GraphVertexProperty vp(nSize,-1,-1);

		boost::add_vertex(vp,*new_graph);


		//对这两个点进行加边操作,
		//node<--->new_node
		GraphEdgeProperty newEP;
		newEP.start_ = nodeId;
		newEP.end_ = nSize;
		newEP.index = newGraphEdgeSize;

		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		//还没对index赋值

		newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());

		boost::add_edge(nodeId,nSize,newEP,*new_graph);


		////12-12-2015
		////两点之间边距离满足某一阈值则进行加边操作，不做边保持不变的假设。
		////遍历所有的节点除了nodeid和nSize;
		//map<IndexType,HVertex*> nodeVtx = new_label_bucket[nodeId]->vertex_bucket;
		//map<IndexType,HVertex*> nSizeVtx = new_label_bucket[nSize]->vertex_bucket;

		//pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);
		//for (VertexIterator fIt = vi.first; fIt != vi.second; ++fIt)
		//{
		//	VertexDescriptor vD = (*fIt);
		//	GraphVertexProperty& vP = (*new_graph)[vD];

		//	IndexType nId = vP.index;
		//	if (nodeId == nId || nSize == nId)
		//	{
		//		continue;
		//	}

		//	map<IndexType,HVertex*> curVtx = new_label_bucket[nId]->vertex_bucket;

		//	ScalarType minNode = 0.0;
		//	ScalarType minSize = 0.0;

		//	minDistBeTwoParts(tgFrame,curVtx,nodeVtx,minNode);
		//	minDistBeTwoParts(tgFrame,curVtx,nSizeVtx,minSize);

		//	if(minNode < 0.05 )
		//	{
		//		GraphEdgeProperty newEP;

		//		if (nId < nodeId)
		//		{
		//			newEP.start_ = nId;
		//			newEP.end_ = nodeId;
		//		} 
		//		else
		//		{
		//			newEP.start_ = nodeId;
		//			newEP.end_ = nId;
		//		}
		//		newGraphEdgeSize = new_graph->m_edges.size();
		//		newEP.index = newGraphEdgeSize;

		//		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		//		//还没对degePoints赋值,只有赋上合理的值才能在下次分裂时起作用。

		//		//newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());

		//		boost::add_edge(nId,nodeId,newEP,*new_graph);
		//	}

		//	if (minSize < 0.05)
		//	{
		//		GraphEdgeProperty newEP;

		//		if (nId < nSize)
		//		{
		//			newEP.start_ = nId;
		//			newEP.end_ = nSize;
		//		} 
		//		else
		//		{
		//			newEP.start_ = nSize;
		//			newEP.end_ = nId;
		//		}
		//		newGraphEdgeSize = new_graph->m_edges.size();
		//		newEP.index = newGraphEdgeSize;

		//		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
		//		//还没对degePoints赋值,只有赋上合理的值才能在下次分裂时起作用。

		//		//newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());

		//		boost::add_edge(nId,nodeId,newEP,*new_graph);
		//	}

		//}//end for vertex

  //     //结束加入循环判断


		//断定查找两个节点与其它节点进行连边操作只会出现 recordColapseEdges.size()次数.
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


			if (sEdgeId < nodeId) //nodeid为终点
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

			}else//nodeid为起点
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

		}//遍历collapse的边


	} //end  while

	checkPsNewLabelParentPtr(new_label_bucket,labParentsize);//next dirction

	map<IndexType,IndexType> labelIndex;
	IndexType kk=0;
	for (auto iter = new_label_bucket.begin(); iter != new_label_bucket.end(); ++ iter,++kk)
	{
		IndexType label = (*iter)->label_id;
		labelIndex[label] = kk;
	}

	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);

	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//保存最新的graph

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
// 			Logger<<"边上的顶点数太少,无法分裂.\n";
// 			continue;
// 		}
// 
// 		map<IndexType,HVertex*> edgeCorrNextVtx;
// 
// 		IndexType newGraphEdgeSize = new_graph->m_edges.size();
// 
// 		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//获得边界点在下一帧对应的块和对应点
// 		//IndexType nodeId = edgeCorrNextVtx.size();
// 
// 
// 		//对应回去,标签是起始点,则保持不变
// 		HLabel* splitedLabel = new_label_bucket[nodeId];
// 
// 		IndexType eS = ep.start_;
// 		IndexType eE = ep.end_;
// 
// 		Logger<<"  边的起点为"<<eS<<"终点为"<<eE<<endl;
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
// 		//若分裂出来的点个数有一个数据很少,则该边不做裂变
// 
// 		if ( recordE < 5 || recordS < 5 )
// 		{
// 			Logger<<"边界点太靠近,不需要分裂.\n";
// 			continue;
// 		}
// 
// 		if (ration < 0.2)
// 		{
// 			Logger<<"Unmark点比值太大,暂时不分裂.\n";
// 			continue;
// 		}
// 
// 		//遍历nodeId 上相连接的边
// 
// 		pair<VertexIterator, VertexIterator> vi = boost::vertices(*new_graph);
// 
// 		VertexIterator nodeIter = (vi.first + nodeId);
// 
// 		VertexDescriptor nodeDesc = *nodeIter;
// 
// 		//节点对应的所有出边
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
// 			boost::remove_edge(*oit,*new_graph);//删除这条边
// 		}
// 
// 		//增加一个节点
// 
// 		IndexType nSize = boost::num_vertices(*new_graph);
// 
// 		GraphVertexProperty vp(nSize,-1,-1);
// 
// 		boost::add_vertex(vp,*new_graph);
// 
// 
// 		//更新分割块信息,新增加的Label标号为nSize. 被分裂的点为nodeId
// 		IndexType new_label = nSize;
// 
// 		new_label_bucket.push_back((HLabel*)0 );
// 		HLabel* new_label_space = allocator_.allocate<HLabel>();
// 		HLabel* new_label_obj = new (new_label_space)HLabel;
// 		new_label_bucket[new_label] = new_label_obj;
// 		new_label_bucket[new_label]->label_id = new_label;
// 		new_label_bucket[new_label]->frame_parent = &hier_componets_[tgFrame];
// 
// 		// 		//对应回去,标签是起始点,则保持不变
// 
// 		map<IndexType,HVertex*> unMakePs;
// 
// 		for (auto iter = splitedLabel->vertex_bucket.begin(); iter != splitedLabel->vertex_bucket.end(); )
// 		{
// 			IndexType prev_id = iter->second->prev_corr->label_parent[gLevel]->label_id; //得不到最新的label_parent地址.
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
// 			}else//一些待确定label的点
// 			{
// 				unMakePs.insert(*iter);
// 				iter = splitedLabel->vertex_bucket.erase(iter);
// 			}
// 
// 		}
// 
// 		//用随机取点产生的最小距离来判断不确定点属于哪个类.unmark 要么属于nodeid 要么属于new_label
// 		determinateUnmarkPoints(tgFrame,unMakePs,new_label_bucket,nodeId,new_label,tgGraSize);
// 
// 
// 		//对这两个点进行加边操作,
// 		//node<--->new_node
// 		GraphEdgeProperty newEP;
// 		newEP.start_ = nodeId;
// 		newEP.end_ = nSize;
// 		newEP.index = newGraphEdgeSize;
// 
// 		IndexType edgeKey = frame_index_to_key(newEP.start_,newEP.end_);
// 		//还没对index赋值
// 
// 		newEP.edgePoints[edgeKey].insert(edgeCorrNextVtx.begin(),edgeCorrNextVtx.end());
// 
// 
// 		boost::add_edge(nodeId,nSize,newEP,*new_graph);
// 
// 		//断定查找两个节点与其它节点进行连边操作只会出现 recordColapseEdges.size()次数.
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
// 			if (sEdgeId < nodeId) //nodeid为终点
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
// 			}else//nodeid为起点
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
// 		}//遍历collapse的边
// 
// 
// 	}//遍历每条边,每条边都会使得new_graph增加一个新的节点
// 
// 
// 	checkPsNewLabelParentPtr(new_label_bucket,labParentsize);//next dirction
// 
// 	hier_componets_[tgFrame].hier_label_bucket.push_back(new_label_bucket);
// 
// 	hier_componets_[tgFrame].hier_graph.push_back(new_graph);//保存最新的graph
// 
// 
// 	Logger<<"  End next split.\n";
// 
// 	Logger<<" .......\n";
}

void DualwayPropagation::generateOrderededges(IndexType srFrame, IndexType tgFrame)
{
	Logger<<"对边进行排序.\n";

	while (!orderedEdgeQ.empty() )
	{
		orderedEdgeQ.pop();
	}

	//向前分裂

	//get the new graph of tgGrame--需要深拷贝


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

	gLevel  = srGraSize - 1;//获取最新的层
	LabelsGraph* srGraLat = hier_componets_[srFrame].hier_graph[gLevel];

	IndexType labParentsize = tgGraSize + 1; //生成的层数

	//Logger<<srFrame<<"帧的第"<<gLevel<<"层边界分割"<<tgFrame<<"的"<<tgGraSize - 1 <<"层"<<endl;

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
			Logger<<"边上的顶点数太少,无法分裂.\n";
			continue;
		}

		map<IndexType,HVertex*> edgeCorrNextVtx;

		IndexType newGraphEdgeSize = new_graph->m_edges.size();

		//并没有更新图结构

		IndexType nodeId = checkNextLabelBucket(edgePoints,edgeCorrNextVtx);//获得边界点在下一帧对应的块和对应点
		//IndexType nodeId = edgeCorrNextVtx.size();


		//对应回去,标签是起始点,则保持不变
		HLabel* splitedLabel = new_label_bucket[nodeId];

		IndexType eS = ep.start_;
		IndexType eE = ep.end_;

		//Logger<<"  边的起点为"<<eS<<"终点为"<<eE<<endl;

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

	Logger<<"结束边排序.\n";
	//需要删除申请的内存!!!
}


void DualwayPropagation::generateOrderPrevEdges(IndexType srFrame, IndexType tgFrame)
{
	while (!orderedEdgeQ.empty() )
	{
		orderedEdgeQ.pop();
	}

	//获取需要更新的graph
	Logger<<" .......\n";
	Logger<<"  Start prev split.\n";

	IndexType srGraphSize = hier_componets_[srFrame].hier_graph.size();
	IndexType tgGraphSize = hier_componets_[tgFrame].hier_graph.size();

	assert(srGraphSize > 0 && tgGraphSize > 0);

	LabelsGraph* oriSpGra = hier_componets_[srFrame].hier_graph[srGraphSize - 1];
	LabelsGraph* shouldSplitGraph = new LabelsGraph(*oriSpGra);

	vector<HLabel* > new_label_bucket;
	copyLabelBucket(new_label_bucket,hier_componets_[srFrame].hier_label_bucket[srGraphSize - 1] );

	//获取指导分割的graph
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

	//对边进行优先级排序

	// guideSplitGraph的每条边引导一次分割
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
			Logger<<"边上的顶点数太少,无法分裂.\n";
			continue;
		}

		Logger<<tgFrame<<"帧的边起点"<<ep.start_<<"与终点"<<ep.end_<<"开始分裂"<<endl;

		map<IndexType,HVertex*> edgeCorrPrevVtx;

		IndexType nSize = boost::num_vertices(*shouldSplitGraph);

		IndexType newGraphEdgeSize = shouldSplitGraph->m_edges.size(); //为了给新增加的边添加序号

		bool isSplit = true;

		IndexType edgePsCorNode = 0;

		isSplit = checkPrevLabelBucket(edgePoints,edgeCorrPrevVtx,edgePsCorNode);

		if ( edgePsCorNode < 0 || edgePsCorNode > nSize - 1)
		{
			Logger<<"边界找到的块越界!.\n";
			break;
		}

		if (!isSplit)
		{
			Logger<<"边界点指向了两个块,不分割.\n";
			continue;
		}

		HLabel* splitedLabel = new_label_bucket[edgePsCorNode]; //可能需要分裂的块

		IndexType curEdgeStart = ep.start_;
		IndexType curEdgeEnd   = ep.end_;
		IndexType strCorPsSzie = 0;
		IndexType endCorPsSize = 0;

		//分类前做一个简单的预判断
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
			//增加两点法向夹角判断,若不大于180°,则添加一条边;否则不添加边操作.

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

			if (vtx.prev_corr)
			{
				Tracer::get_instance().add_record(f, vtx.vtx_id, f-1, vtx.prev_corr->vtx_id);
			}

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

void DualwayPropagation::mergePatchTraj()
{

   //calculateSimilar2Componet();

	vector<PatchTraj> ptnodes;


	buildPatchCorrespondenceByLabel();

	generTrajNodes(ptnodes);//

	vector<IndexType> Labels;

	Labels.resize(ptnodes.size(),0);

	graphCuts(ptnodes,Labels);

	mergeSquences(Labels);
}


void DualwayPropagation::generTrajNodes(vector<PatchTraj>& pNodes)
{

	map<IndexType,bool>  isPatchTrav;//记录块是否访问过--key = frameId  & labelId

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
	{
	    IndexType gLevel = fIter->second.hier_label_bucket.size();//访问最高层 图结构

		--gLevel;

		IndexType fId = fIter->first;
		vector<HLabel*>& label_buctet = fIter->second.hier_label_bucket[gLevel];

		map<IndexType,IndexType> labelIdx = fIter->second.hier_label_vtxBucket_index[gLevel];

		IndexType lId = 0;

		for (auto lIter = label_buctet.begin(); lIter != label_buctet.end(); lIter ++)
		{
			lId = (*lIter)->label_id;

			IndexType flKey = frame_label_to_key(fId,lId);

			if (isPatchTrav[flKey])
			{
				continue;
			}

			isPatchTrav[flKey] = true;

			IndexType beforFrame = fId;

			PatchTraj tempNode;

			tempNode.label_id = lId;

			tempNode.startFrame = fId;

			tempNode.endFrame = fId;

			IndexType labId = labelIdx[lId];

			HLabel* nextPatchPtr = label_buctet[labId]->next_corr;

			while (nextPatchPtr != NULL)
			{
			   IndexType nFId = nextPatchPtr->frame_parent->frame_id;

			   IndexType nextflKey = frame_label_to_key(nFId,lId);

			   isPatchTrav[nextflKey] = true;

			   tempNode.endFrame = nFId;

			   Matrix34 toTrans;
			   Matrix34 backTrans;

			   toTrans.setZero();
			   backTrans.setZero();

			   calculateTrans(lId, beforFrame,nFId,toTrans,backTrans);

			   tempNode.fNode.push_back(toTrans);

			   tempNode.bNode.push_back(backTrans);

			   beforFrame = nFId;

			   nextPatchPtr = nextPatchPtr->next_corr;
			}

			pNodes.push_back(tempNode);
		}

	}
}

void DualwayPropagation::graphCuts(vector<PatchTraj>& pNodes, vector<IndexType>& labels)
{
   IndexType nodeSize = pNodes.size();

   GCoptimizationGeneralGraph* segGraphC = new GCoptimizationGeneralGraph(nodeSize,nodeSize,&hier_componets_);

   setSegNeihbor(pNodes,*segGraphC);//设置边界,同时设置数据项  //setSegDataItem(*segGraphC);

   //calculateSimilar2Componet();

   setSegSmoothItem(*segGraphC);//用形状统计图--默认值为0

   segGraphC->expansion(2);

   segGraphC->swap(1);

   getSegLabels(*segGraphC,labels);


}

void DualwayPropagation::mergeSquences(vector<IndexType>& labes)
{

}

void DualwayPropagation::setSegNeihbor(vector<PatchTraj>& pNodes, GCoptimizationGeneralGraph& segGraphC)
{
	IndexType i = 0;
	IndexType j = 0;
	for (auto vIter = pNodes.begin(); vIter != (pNodes.end() - 1); ++ vIter, ++ i)
	{
		auto vvIter = (vIter + 1);

		segGraphC.setLabel(i,i);

		ScalarType dValue = 1000.0;
		segGraphC.setDataCost(i,i,0.);

		j = i + 1;
        for (; vvIter != pNodes.end(); ++ vvIter,++j )
        {
			if (isAdjInSeq(pNodes[i],pNodes[j]) )
			{
				segGraphC.setNeighbors(i,j);

				dValue = motionSimilarityBetw2Nodes(i,j,pNodes);

				segGraphC.setDataCost(i,j,dValue);

			}else
			{
				segGraphC.setDataCost(i,j,1e5);
			}
        }
	}

	segGraphC.setLabel(i,i); //the end() node

	segGraphC.printfNeig();

}

void DualwayPropagation::setSegDataItem(GCoptimizationGeneralGraph& segGraphC)
{
	segGraphC.setDataCost(0,1,0.5);
	
}

void DualwayPropagation::setSegSmoothItem(GCoptimizationGeneralGraph& segGraphC)
{
//运用SDF值来确定两块之间的相似性
}

bool DualwayPropagation::isAdjInSeq(PatchTraj& nodeA, PatchTraj& nodeB)
{
    //只要在共同帧上存在相邻关系,则认为这两条轨迹相邻.

	IndexType lab_A = nodeA.label_id;
	IndexType lab_B = nodeB.label_id;

	IndexType com_str = max(nodeA.startFrame, nodeB.startFrame);
	IndexType com_end = min(nodeA.endFrame, nodeB.endFrame);

	if (com_str > com_end)
	{
		return  false;
	}

	for (IndexType i = com_str; i < com_end; ++ i)
	{
		IndexType gLevel = hier_componets_[i].hier_graph.size();
		--gLevel;

		LabelsGraph& curGraph = *hier_componets_[i].hier_graph[gLevel];

		IndexType aIdx = hier_componets_[i].hier_label_vtxBucket_index[gLevel][lab_A];
		IndexType bIdx = hier_componets_[i].hier_label_vtxBucket_index[gLevel][lab_B];

		pair<VertexIterator, VertexIterator> vi = boost::vertices(curGraph);

		VertexIterator nodeIterA = (vi.first + aIdx);
		VertexIterator nodeIterB = (vi.first + bIdx);

		VertexDescriptor nodeDescA = *nodeIterA;
		VertexDescriptor nodeDescB = *nodeIterB;

		pair<EdgeDescriptor,bool> isAdj = edge(nodeDescA,nodeDescB,curGraph);

		if (isAdj.second) 
		{
			return true;
		}

	}
	

	return false;
}

void DualwayPropagation::getSegLabels(GCoptimizationGeneralGraph& segGraphC, vector<IndexType>& labels)
{
	IndexType vSize = labels.size();
	for (IndexType i = 0; i < vSize; ++ i)
	{
		labels[i] = segGraphC.whatLabel(i);

		printf("%d Node's label is %d.\n",i,labels[i]);

	}
}

ScalarType DualwayPropagation::motionSimilarityBetw2Nodes(IndexType i, IndexType j, vector<PatchTraj>& oriTraj)
{
	//Frobenius
	IndexType strF = max(oriTraj[i].startFrame, oriTraj[j].startFrame);

	IndexType endF = min(oriTraj[i].endFrame, oriTraj[j].endFrame);

	if (strF >= endF)
	{
		return 1e5;
	}

	ScalarType fDis = -1e5;

	IndexType iStr = oriTraj[i].startFrame;
	IndexType jStr = oriTraj[j].startFrame;

	for (IndexType fid = strF; fid < endF; ++ fid)
	{
		Matrix34 iTrans;
		Matrix34 jTrans;

		iTrans = oriTraj[i].fNode[fid - iStr];
		jTrans = oriTraj[j].fNode[fid - jStr];

		ScalarType tDis = (iTrans - jTrans).norm();

		if (tDis > fDis)
		{
			fDis = tDis;
		}
	}

	return fDis;
}

void DualwayPropagation::calculateTrans(IndexType lab,IndexType sFrame, IndexType tFrame, Matrix34& toTrans, Matrix34& backTrans)
{
	Matrix3X srCoor;
	Matrix3X tgCoor;


    getCoorByVtxBucket(lab,sFrame,srCoor);
	getCoorByVtxBucket(lab,tFrame,tgCoor);

	IndexType oriSize = srCoor.cols();
	IndexType trgSize = tgCoor.cols();


	DeformableRegistration nonrigid;

	MatrixXXi vtxMap;
	Matrix3X oriCopyCoor = srCoor;
	Matrix3X tgCopyCoor = tgCoor;

	vtxMap.resize(1,oriSize );
	SICP::Parameters pa(false,2.0,10,1.2,1e5,20,20,1,1e-5); 

	SICP::point_to_point(oriCopyCoor,tgCoor,vtxMap,pa);

	nonrigid.alignTargetCoorChangeSize(tgCoor,vtxMap);

	nonrigid.calculateAffineTrans(srCoor,tgCoor,toTrans);

	//需要经过配准才可以计算仿射变换
	vtxMap.resize(1,trgSize);

	tgCoor.resize(3,trgSize);

	tgCoor = tgCopyCoor;

	SICP::point_to_point(tgCopyCoor,srCoor,vtxMap,pa);

	nonrigid.alignTargetCoorChangeSize(srCoor,vtxMap);

	nonrigid.calculateAffineTrans(tgCoor,srCoor,backTrans);


}

void DualwayPropagation::getCoorByVtxBucket(IndexType lab, IndexType frameId, Matrix3X& vtxCoor)
{
	IndexType gLev = hier_componets_[frameId].hier_graph.size();
	--gLev;

	map<IndexType,IndexType> labMap = hier_componets_[frameId].hier_label_vtxBucket_index[gLev];

	HLabel* label_buctet = hier_componets_[frameId].hier_label_bucket[gLev][labMap[lab]];

	Sample& smp = SampleSet::get_instance()[frameId];

	auto vIter = label_buctet->vertex_bucket.begin();

	auto vEnd = label_buctet->vertex_bucket.end();

	IndexType vSize = label_buctet->vertex_bucket.size();

	vtxCoor.resize(3,vSize);

	IndexType i = 0;
	for (; vIter != vEnd; ++ vIter, ++ i)
	{
		IndexType vId = vIter->first;

		vtxCoor.col(i) = smp.vertices_matrix().col(vId);
	}

}

void DualwayPropagation::initGraphAfterCoseg()
{
	//HFrame 里面有vtx_bucket 和 vtxBucketMap, 但没有graph; 
	//coseg之后图不变,但图对应的节点已经变了,需要调整vector的顺序


	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
	{
		IndexType bLevel = fIter->second.hier_label_bucket.size();//访问最高层 图结构
		--bLevel;

		IndexType gLevel = fIter->second.hier_graph.size();  //此时,graph比bucket少一层
	    --gLevel;

		LabelsGraph* nGraph =  new LabelsGraph(  *(fIter->second.hier_graph[gLevel] ) );//获取当前最新的图结构 

		fIter->second.hier_graph.push_back(nGraph); //赋旧图结构

		vector<HLabel*> temp = fIter->second.hier_label_bucket[bLevel]; //最上层的bucket

		vector<HLabel*> secBucket = fIter->second.hier_label_bucket[bLevel - 1]; //倒数第二层的bucket

		map<IndexType, IndexType> secbucketIdx = fIter->second.hier_label_vtxBucket_index[bLevel - 1];

		map<IndexType, IndexType> lasterbucketIdx = fIter->second.hier_label_vtxBucket_index[bLevel];

		vector<HLabel*> midBucket;
		midBucket.resize(secBucket.size() );

		map<IndexType, IndexType> midBucketIdx;
		midBucketIdx.clear();

		IndexType i = 0;
		for (auto hiter = secBucket.begin(); hiter != secBucket.end(); ++ hiter,++ i)
		{

			HVertex* tV= (*(*hiter)->vertex_bucket.begin()).second;

			IndexType vtxId = tV->vtx_id;

			IndexType lasterLab = fIter->second.label_of_vtx[vtxId];

			IndexType labIdx = lasterbucketIdx[lasterLab];

			//IndexType interIdxM = secbucketIdx[iterLab];
// 			midBucketIdx[lasterLab] = interIdxM; 
// 			midBucket[interIdxM] = temp[labIdx];

			midBucketIdx[lasterLab] = i;

			midBucket.push_back(temp[labIdx] );
		}

		fIter->second.hier_label_bucket[bLevel] = midBucket;

		fIter->second.hier_label_vtxBucket_index[bLevel] = midBucketIdx;

	}
}

void DualwayPropagation::calculateSimilar2Componet()
{
	patchSimilarMea.clear();

	Engine* ep;

	if (! (ep = engOpen(NULL)) )
	{
		Logger<< "Can't not start Matlab engine.\n";
		return;
	}

	//Get the executable file's path
	char cur_path[FILENAME_MAX];
	if (!get_current_dir(cur_path, sizeof(cur_path)))
	{
		Logger<<"read error current path.\n";
		return;
	}

	cur_path[sizeof(cur_path) - 1] = '\0';
	strcat(cur_path,"\\sdf_dist");
	char cd_cur_path[FILENAME_MAX + 3] = "cd ";
	strcat(cd_cur_path, cur_path);
	engEvalString(ep, cd_cur_path );

	//函数调用 "sdf_dist_matrix( points, normals, patches, patch_to_do, outliers, idx, flat_strictness )"

	for (auto fIter = hier_componets_.begin(); fIter != hier_componets_.end(); fIter ++)
	{
		IndexType bLevel = fIter->second.hier_label_bucket.size();//访问最高层 图结构

		--bLevel;

		Sample& tempS = SampleSet::get_instance()[fIter->first];

		map<IndexType,IndexType>& vId2Graph = hier_componets_[fIter->first].gId_of_vtx;

		map<IndexType,IndexType>& labVtx = hier_componets_[fIter->first].label_of_vtx;

		IndexType sVSize = vId2Graph.size();
		IndexType oriVSize = tempS.num_vertices();

		MatrixXXDR sCoorMat;

		MatrixXXDR sNormMat;

// 		sCoorMat.resize(3,oriVSize);
// 		sNormMat.resize(3,oriVSize);

// 		for (auto mBeg = vId2Graph.begin(); mBeg != vId2Graph.end(); ++ mBeg)
// 		{
// 			IndexType vIdOri = mBeg->first;
// 			IndexType vIdGra = mBeg->second;
// 			sCoorMat.col(vIdGra) << tempS.vertices_matrix()(0,vIdOri),tempS.vertices_matrix()(1,vIdOri),tempS.vertices_matrix()(2,vIdOri);
// 			sNormMat.col(vIdGra) << tempS.nor_matrix()(0,vIdOri),tempS.nor_matrix()(1,vIdOri),tempS.nor_matrix()(2,vIdOri);
// 		}


		//传入原始点云数据
		sCoorMat.resize(3,oriVSize);

		sNormMat.resize(3,oriVSize);

		sCoorMat = tempS.vertices_matrix().transpose().cast<double>();

		sNormMat = tempS.nor_matrix().transpose().cast<double>();

		mxArray* mxCoor = mxCreateDoubleMatrix(oriVSize, 3, mxREAL);

		mxArray* mxNorm = mxCreateDoubleMatrix(oriVSize, 3, mxREAL);

		memcpy( (char*)mxGetPr(mxCoor), (char*)sCoorMat.data(), 3*oriVSize*sizeof(double) );

		engPutVariable(ep,"points",mxCoor);

		memcpy((char*)mxGetPr(mxNorm), (char*)sNormMat.data(), 3*oriVSize*sizeof(double) );

		engPutVariable(ep,"normals",mxNorm);


		vector<HLabel*> label_buctet = hier_componets_[fIter->first].hier_label_bucket[bLevel];

		IndexType vBSize = label_buctet.size();

		auto vIter = label_buctet.begin();

		auto vEnd = label_buctet.end();

		mwSize * te = new mwSize[vBSize];

		mwSize* ndim = new mwSize[1];

		ndim[0] = vBSize;

		mxArray *mxPatches = mxCreateCellArray(1, ndim);//(维数,各个维度的长度)

		Eigen::VectorXd patchTodoVS;

		patchTodoVS.resize(vBSize,1);

		Eigen::VectorXd vIdAll;

		vIdAll.resize(sVSize,1);

		IndexType iV = 0;

		for (; vIter != vEnd; ++ vIter, ++ iV)
		{
			auto vtxBktBeg =  (*vIter)->vertex_bucket.begin();

			auto vtxBktEnd =  (*vIter)->vertex_bucket.end();

			IndexType vSize = (*vIter)->vertex_bucket.size();

			te[iV] = vSize;
		    
			patchTodoVS[iV] = (iV + 1);

			Eigen::MatrixXd tempId;
			tempId.resize(vSize,1);

			mxArray* tempMat = mxCreateDoubleMatrix(vSize,1,mxREAL);

			IndexType iId = 0;
			for (; vtxBktBeg != vtxBktEnd; ++ vtxBktBeg, ++iId)
			{
				IndexType vtxId = vtxBktBeg->first;
				tempId(iId,0) = (vtxId + 1);//原始点云索引

				//IndexType gVtxId = vId2Graph[vtxId];
				//vIdAll[gVtxId] = (iV + 1);//采样点的分类结果
				//tempId(iId,0) = (gVtxId + 1);//graph 上点索引
			}
			
			memcpy( (char*)mxGetPr(tempMat), (char*)tempId.data(), vSize * sizeof(double) );
			
			mxSetCell(mxPatches,iV,tempMat);

		}

		engPutVariable(ep,"patches",mxPatches);

	    mxArray* mxPatchTodo  = mxCreateDoubleMatrix(1, vBSize,mxREAL);

		memcpy((char*)mxGetPr(mxPatchTodo), (char*)patchTodoVS.data(), vBSize * sizeof(double) );

		engPutVariable(ep,"patch_to_do",mxPatchTodo);

		mxArray* outlier = mxCreateDoubleMatrix(1,1,mxREAL);

		engPutVariable(ep,"outliers",outlier);

//      采样点分类结果
// 		mxArray* pointId = mxCreateDoubleMatrix(1, sVSize,mxREAL);
// 
// 		memcpy((char*)mxGetPr(pointId), (char*)vIdAll.data(), sVSize * sizeof(double) );
// 
// 		engPutVariable(ep,"idx",pointId);


		vector<IndexType> smpVtxId;
		vector<IndexType> smpvtxLab;
		vector<double> oriVtxId;
		oriVtxId.resize(oriVSize,0);

		for (auto vIter = labVtx.begin(); vIter != labVtx.end(); ++vIter)
		{
			smpVtxId.push_back(vIter->first);
			smpvtxLab.push_back(vIter->second);
		}

		propagateLabel2Orignal(tempS,smpVtxId,smpvtxLab,oriVtxId);

		mxArray* pointId = mxCreateDoubleMatrix(1, oriVSize,mxREAL);

		memcpy((double*)mxGetPr(pointId), (double*)oriVtxId.data(), oriVSize * sizeof(double) );

		engPutVariable(ep,"idx",pointId);

		mxArray* flat_st = mxCreateDoubleScalar(0.5);

		engPutVariable(ep,"flat_strictness",flat_st);

		engEvalString(ep,"dist_matrix = sdf_dist_matrix( points, normals, patches, patch_to_do, outliers, idx, flat_strictness );" );

		mxArray* mxSimilar = NULL;

		mxSimilar = engGetVariable(ep,"dist_matrix");

		double* patchSimi = (double*)mxGetData(mxSimilar);

		Eigen::MatrixXd sigFrameSimi;

		sigFrameSimi.resize(vBSize,vBSize);

		memcpy((char*)(&sigFrameSimi), (char*)(patchSimi), vBSize * vBSize * sizeof(double) );

		patchSimilarMea.push_back(sigFrameSimi);

		mxDestroyArray(mxCoor);
		mxDestroyArray(mxNorm);
		//mxDestroyArray(mxPatches);
		mxDestroyArray(mxPatchTodo);
		mxDestroyArray(mxSimilar);//?
	}
		engClose(ep);

}

void DualwayPropagation::propagateLabel2Orignal(Sample& oriPC,vector<IndexType>& sampleVtxId,vector<IndexType>& label_smp,vector<double>& label_ori)
{
	Logger<<"Start paopa label to ori.\n";

	IndexType nCluster = 35;

	map<IndexType,IndexType> smpLabel;
	map<IndexType,IndexType>::iterator IsValidIter;
	for (int i = 0; i < label_smp.size(); i++)
	{
		smpLabel.insert(make_pair(sampleVtxId[i],label_smp[i]));
	}

	const IndexType k = 250;
	IndexType neighbours[k];
	ScalarType dist[k];

	IndexType vtx_num = oriPC.num_vertices();
	IndexType result_label;


	for(IndexType vtx_id = 0; vtx_id < vtx_num; vtx_id ++)
	{
		vector<IndexType> recordLabelTime(nCluster,0);
		result_label = -1;
		oriPC.neighbours(vtx_id, k, neighbours, dist);
		for(IndexType neig_id = 0; neig_id < k; neig_id ++)
		{
			IsValidIter = smpLabel.find(neighbours[neig_id]);
			if(IsValidIter != smpLabel.end())
			{
				recordLabelTime[IsValidIter->second] += 1;
			}
		}
		for (int i = 0; i < nCluster; i++)
		{
			if(result_label < recordLabelTime[i])
			{
				//label_ori[vtx_id] = i;

				label_ori[vtx_id] = (i + 1.);//为了把数据传输给matlab,label标签加1;
				result_label = recordLabelTime[i];
			}
		}

	}

	smpLabel.clear();

	Logger<<"End propagation label to ori .\n";
}