#include "GraphMatching.h"

void GraphMatch::calculateSimilar2Frame()
{
	if( components_.find(srFrame+1) == components_.end())return;

	HFrame& cur_frame = components_[srFrame];

	IndexType srLevel = components_[srFrame].hier_graph.size();
	IndexType tgLevel = components_[srFrame].hier_graph.size();

	srLevel--;
	tgLevel--;

	LabelsGraph* labelGraphOfCFrame = cur_frame.hier_graph[srLevel];

	//���������
	LabelsGraph::vertex_iterator CFvertex_biter ,CFvertex_eiter;

	LabelsGraph::adjacency_iterator CFadj_bitr , CFadj_eitr;

	tie( CFvertex_biter ,CFvertex_eiter) = boost::vertices(* labelGraphOfCFrame);

	for(; CFvertex_biter != CFvertex_eiter ; ++CFvertex_biter)
	{
		//�Ӷ���������ж�

		if (!equeue_.empty())
		{
			while (!equeue_.empty())
			{
				equeue_.pop();
			}
		}

		IndexType vtx_id = *CFvertex_biter;


		IndexType vtxSize = components_[srFrame].hier_label_bucket[srLevel][vtx_id]->vertex_bucket.size();

		if ( vtxSize < 5) //�����̫����ȥ�Ҷ�Ӧ��
		{
			Logger<<"  �鱾ʡ̫С,��ʱ���Ҷ�Ӧ!.\n";
			continue;
		}

		boost::tie(CFadj_bitr , CFadj_eitr) = boost::adjacent_vertices( *CFvertex_biter ,*labelGraphOfCFrame ); // ����ڽӵ�ĵ�����

		std::set<set<IndexType> > unique_pair_set;

		//�����ڵ� ƥ��
		std::set<IndexType> CFcmp_node ,NFcmp_node;  //�ӿڶ��� �����ڴ������ڼ����������Ͻڵ�

		std::set<set<IndexType> > comp_set;

		CFcmp_node.clear();

		CFcmp_node.insert(vtx_id);

		matchNextFrameOnePt(CFcmp_node,srLevel,tgLevel);

		//�ҵ�����ƥ��ֱ�Ӷ�HLabel�Ķ�Ӧ��ֵ
		if (equeue_.empty())
		{
			continue;
		}

		NFcmp_node = equeue_.top().tgPatches;

		ScalarType val = equeue_.top().Er;

		if (val > 10)
		{
			Logger<<" û���ҵ����ʵĿ�.\n";
			continue;
		}

		if (NFcmp_node.size() > 1)
		{
			Logger<<"Merge this patches to match label"<<vtx_id<<endl;

		}

		//build patch correspondece		

		IndexType corLabelId = (*NFcmp_node.begin());

		Logger<<srFrame<<"֡�Ŀ�"<<vtx_id<<"��֡"<<tgFrame<<"�Ŀ�"<<corLabelId<<"��Ӧ.\n";

		HLabel* nextPtPtr = components_[tgFrame].hier_label_bucket[tgLevel][corLabelId]; 

		components_[srFrame].hier_label_bucket[srLevel][vtx_id]->next_corr = nextPtPtr;

		// 			comp_set.insert(CFcmp_node);
		// 
		// 			unique_pair_set.insert(CFcmp_node);
		// 
		// 			CFrameSetSet_.insert(comp_set.begin() , comp_set.end() ); 
		// 			sadj_map_.insert( std::make_pair(* CFvertex_biter ,CFcmp_node));
		// 			stotal_map_.insert( std::make_pair( * CFvertex_biter , comp_set));
	}
}

void GraphMatch::matchNextFrameOnePt(set<IndexType>& CFcmp_node, IndexType srGraLevel,IndexType tgGraLevel)
{
	LabelsGraph::vertex_iterator NFvertex_biter ;
	LabelsGraph::vertex_iterator NFvertex_eiter;
	LabelsGraph::adjacency_iterator  NFadj_bitr ,NFadj_eitr;

	IndexType gLevel = components_[tgFrame].hier_graph.size();

	LabelsGraph* labelGraphOfNFrame = components_[tgFrame].hier_graph[tgGraLevel];

	set<IndexType> NFcmp_node;

	tie( NFvertex_biter ,NFvertex_eiter) = boost::vertices(* labelGraphOfNFrame);

	for(; NFvertex_biter != NFvertex_eiter; ++NFvertex_biter)
	{
		IndexType vtx_num = components_[tgFrame].hier_label_bucket[tgGraLevel][*NFvertex_biter]->vertex_bucket.size();
		
		if (vtx_num < 5)
		{
			continue;
		}

		boost::tie( NFadj_bitr , NFadj_eitr) = boost::adjacent_vertices( *NFvertex_biter ,*labelGraphOfNFrame ); // ����ڽӵ�ĵ�����

		std::set<set<IndexType> > comp_set;	
		NFcmp_node.clear();
		NFcmp_node.insert(*NFvertex_biter);
		comp_set.insert(NFcmp_node);
		NFrameSetSet_.insert(comp_set.begin() , comp_set.end() );
		tadj_map_.insert( std::make_pair(* NFvertex_biter ,NFcmp_node));
		ttotal_map_.insert( std::make_pair( * NFvertex_biter , comp_set));

		//ScalarType  error = distance2PatchesLevel(CFcmp_node, NFcmp_node, srGraLevel, tgGraLevel);

		ScalarType error = distBetween2Frame4Items(CFcmp_node, NFcmp_node, srGraLevel, tgGraLevel);

		if (error < INF)
		{
	    	equeue_.push( PatchMatch( CFcmp_node ,NFcmp_node, error ) );
		}

	}
}

ScalarType GraphMatch::distBetween2Frame4Items(const set<IndexType>& sPatches, const set<IndexType>& tPatches, IndexType srLevel, IndexType tgLevel)
{
	assert(sPatches.size() > 0 && tPatches.size() > 0);

	map<IndexType,HVertex*> srVertex;

	set<IndexType>::iterator sIter = sPatches.begin();


	for (; sIter != sPatches.end(); sIter ++)
	{

		HLabel* temp = components_[srFrame].hier_label_bucket[srLevel][*sIter];

		srVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
	}

	map<IndexType,HVertex*> tgVertex;

	set<IndexType>::iterator tIter = tPatches.begin();

	for (; tIter != tPatches.end(); tIter++)
	{
		HLabel* temp = components_[tgFrame].hier_label_bucket[tgLevel][*tIter];
		tgVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
	}

	map<IndexType,HVertex*> toComVtx; // ��¼��Щ�����ǵĶ�Ӧ�����ڶ�Ӧ����

	map<IndexType,HVertex*> backComVtx; // �໥��¼

	bool isTo = true;
	computer_common(srVertex,tgVertex,toComVtx,isTo);

	isTo = false;
	computer_common(srVertex, tgVertex,backComVtx,isTo);


	ScalarType srV_size = srVertex.size();
	ScalarType tgV_size = tgVertex.size();

	IndexType toCom_size = toComVtx.size();
	IndexType backCom_size = backComVtx.size();

	ScalarType corValue =  0.5 * (toCom_size/srV_size + backCom_size/tgV_size );

	ScalarType scaleShapeValue =  abs( srV_size - tgV_size) / max(srV_size, tgV_size );


	Logger<<" Source combine patches.\n";
	for (auto  iter = sPatches.begin(); iter != sPatches.end(); ++ iter)
	{
		Logger<<(*iter)<<"  ";
	}
	Logger<<".\n";



	Logger<<" Target combine patches.\n";
	for (auto  iter = tPatches.begin(); iter != tPatches.end(); ++ iter)
	{
		Logger<<(*iter)<<"  ";
	}
	Logger<<".\n";

	//ScalarType totDis = corValue + scaleShapeValue + deformationValue;

	Logger<<" Two condition type distances.\n";
	Logger<<" corValue = "<<corValue<<endl;
	Logger<<" scaleShapeValue = "<<scaleShapeValue<<endl;
	//Logger<<" deformationvalue = "<<deformationValue<<endl;

	if ( scaleShapeValue > 0.2 || corValue < 0.6 || toComVtx.size() < 5 || backComVtx.size() < 5)
	{
		Logger<<" The distance is INF.\n";

		return INF;

	}else
	{
		//ScalarType deformationValue  = deformationDis(toComVtx,backComVtx);//���ö�Ӧ��
		ScalarType deformationValue  = deformationDis(srVertex, tgVertex); //����ȫ���ĵ�

		Logger<<" deformationvalue = "<<deformationValue<<endl;

		return deformationValue + scaleShapeValue  +  (1 - corValue);

		//return  (1. - corValue) + deformationValue;
	}
}

ScalarType GraphMatch::deformationDis(map<IndexType,HVertex*>& toComVtx, map<IndexType,HVertex*>& backComVtx)
{
     //ScalarType dis = toDeformationErr(toComVtx) + backDeformationErr(backComVtx);//ֻ���ö�Ӧ��

	 ScalarType dis = toDeformationErrAllVtx(toComVtx,backComVtx) + backDeformationErrAllVtx(toComVtx,backComVtx);//local ICP

	 return dis;
}

ScalarType GraphMatch::toDeformationErr(map<IndexType,HVertex*>& toComVtx)
{
	Sample& oriFrame = sample_set_[srFrame];
	Sample& tarFrame = sample_set_[tgFrame];

	IndexType toPsSize = toComVtx.size();


	Matrix3X s_coord, t_coord;
	s_coord.setZero(3, toPsSize);
	t_coord.setZero(3, toPsSize);

	IndexType i = 0;
	for (auto toIter = toComVtx.begin(); toIter != toComVtx.end(); ++toIter,++i)
	{

		IndexType sPsId = toIter->first;
		IndexType corPsId = toIter->second->next_corr->vtx_id;

		s_coord.col(i) = oriFrame.vertices_matrix().col(sPsId);
		t_coord.col(i) = tarFrame.vertices_matrix().col(corPsId);

	}

	Matrix33 rot_mat;
	MatrixXX tran_vec;

	point2point(s_coord, t_coord, rot_mat, tran_vec);

	//����ÿ���㵽��Ӧ�����
	MatrixXX traSCoor = rot_mat * s_coord;

	traSCoor.colwise() += tran_vec.col(0);

	MatrixXX errMat = t_coord - traSCoor;

	ScalarType totDis = errMat.colwise().norm().sum();

	return totDis;
	//return totDis/i;
}

ScalarType GraphMatch::backDeformationErr(map<IndexType,HVertex*>& backeComVtx)
{
	Sample& oriFrame = sample_set_[tgFrame];

	Sample& tarFrame = sample_set_[srFrame];


	IndexType backPsSize = backeComVtx.size();

	Matrix3X s_coord, t_coord;
	s_coord.setZero(3, backPsSize);
	t_coord.setZero(3, backPsSize);

	IndexType i = 0;
	for (auto backIter = backeComVtx.begin(); backIter != backeComVtx.end(); ++backIter,++i)
	{

		IndexType sPsId = backIter->first;
		IndexType corPsId = backIter->second->prev_corr->vtx_id;

		s_coord.col(i) = oriFrame.vertices_matrix().col(sPsId);
		t_coord.col(i) = tarFrame.vertices_matrix().col(corPsId);

	}

	Matrix33 rot_mat;
	MatrixXX tran_vec;

	point2point(s_coord, t_coord, rot_mat, tran_vec);

	//����ÿ���㵽��Ӧ�����
	MatrixXX traSCoor = rot_mat * s_coord;

	traSCoor.colwise() += tran_vec.col(0);

	MatrixXX errMat = t_coord - traSCoor;

	ScalarType totDis = errMat.colwise().norm().sum();

	return totDis;//ƫ�������С�Ŀ�
	//return totDis/i;//ÿ�����λ��ͬ
}

void  GraphMatch::point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec)
{
	Eigen::Vector3f X_mean, Y_mean;

	for(int i=0; i<3; ++i) //���������Ƶľ�ֵ
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
		transformation.linear().noalias() = svd.matrixV()*svd.matrixU().transpose();//������ת����
	}

	transVec = Y_mean - transformation.linear()*X_mean;
	rotMat = transformation.linear() ;

	srCloud.colwise() += X_mean;
	tgCloud.colwise() += Y_mean;
}

ScalarType GraphMatch::distance2PatchesLevel(const set<IndexType>& sPatches,const set<IndexType>& tPatches,IndexType srLevel, IndexType tgLevel)
											 
{
	assert(sPatches.size() > 0 && tPatches.size() > 0);

	map<IndexType,HVertex*> srVertex;

	set<IndexType>::iterator sIter = sPatches.begin();

	for (; sIter != sPatches.end(); sIter ++)
	{
		HLabel* temp = components_[srFrame].hier_label_bucket[srLevel][*sIter];

		srVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
	}

	map<IndexType,HVertex*> tgVertex;

	set<IndexType>::iterator tIter = tPatches.begin();

	for (; tIter != tPatches.end(); tIter++)
	{
		HLabel* temp = components_[tgFrame].hier_label_bucket[tgLevel][*tIter];
		tgVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
	}

	map<IndexType,HVertex*> toComVtx; // ��¼��Щ�����ǵĶ�Ӧ�����ڶ�Ӧ����

	map<IndexType,HVertex*> backComVtx; // �໥��¼

	bool isTo = true;
	computer_common(srVertex,tgVertex,toComVtx,isTo);

	isTo = false;
	computer_common(srVertex, tgVertex,backComVtx,isTo);


	//ScalarType corValue =  0.5 * (toComVtx.size()/srVertex.size() + backComVtx.size()/tgVertex.size() );

	//ScalarType scaleShape =  abs((int) (srVertex.size() - tgVertex.size() )) / std::max(srVertex.size(), tgVertex.size() );

	ScalarType srV_size = srVertex.size();
	ScalarType tgV_size = tgVertex.size();

	IndexType toCom_size = toComVtx.size();
	IndexType backCom_size = backComVtx.size();

	ScalarType corValue =  0.5 * (toCom_size/srV_size + backCom_size/tgV_size );

	ScalarType scaleShape =  abs( srV_size - tgV_size) / max(srV_size, tgV_size );

	if ( scaleShape > 0.5 || corValue < 0.5)
	{
		return INF;
	}else
	{
		return  1. - corValue;
	}

}

void GraphMatch::computer_common(map<IndexType,HVertex*>& srPatchesVtx, map<IndexType,HVertex*>& tgPatchesVtx, 
								 map<IndexType,HVertex*>& comVtx,bool isTo)
{
	if (isTo) // -->
	{
		map<IndexType,HVertex*>::iterator iter = srPatchesVtx.begin();

		for (; iter != srPatchesVtx.end(); iter ++)
		{
			if( tgPatchesVtx.find( iter->second->next_corr->vtx_id) != tgPatchesVtx.end())
			{
				comVtx.insert(*iter);
			}

		}

	}else// <--
	{
		map<IndexType,HVertex*>::iterator iter = tgPatchesVtx.begin();

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

void GraphMatch::setMatchset(IndexType srLevel, IndexType tgLevel)
{
	if( components_.find(srFrame) == components_.end())return;

	HFrame& cur_frame = components_[srFrame];
	HFrame& next_frame = components_[tgFrame];

	LabelsGraph* curGraph = cur_frame.hier_graph[srLevel];
	LabelsGraph* nextGraph = next_frame.hier_graph[tgLevel];

	//���������
	LabelsGraph::vertex_iterator CFvertex_biter ,CFvertex_eiter, NFvertex_biter,NFvertex_eiter;
	LabelsGraph::adjacency_iterator CFadj_bitr , CFadj_eitr ,NFadj_bitr , NFadj_eitr;

	tie(CFvertex_biter, CFvertex_eiter) = boost::vertices(*curGraph);

	for(; CFvertex_biter != CFvertex_eiter ; ++CFvertex_biter)
	{
		boost::tie(CFadj_bitr , CFadj_eitr) = boost::adjacent_vertices( *CFvertex_biter, *curGraph); // ����ڽӵ�ĵ�����
		std::set<set<IndexType> > unique_pair_set;

		//�����ڵ� ƥ��
		std::set<IndexType> curNodeSet, nextNodeSet;  //�ӿڶ��� �����ڴ������ڼ����������Ͻڵ�
		std::set<set<IndexType> > comp_set;

		curNodeSet.clear();
		curNodeSet.insert(*CFvertex_biter);

		if( CFadj_bitr == CFadj_eitr )/** ˵���ýڵ�û���ڽӽڵ㣬�����ɿ�**/
		{

		 	std::set<set<IndexType> > comp_set;
		 
		 	curNodeSet.clear();
		 	curNodeSet.insert(*CFvertex_biter);
		 
		 	// ������һ֡��ѭ��
		 	matchNextFrameSet(curNodeSet,*nextGraph,srLevel,tgLevel);
		 
		 	comp_set.insert(curNodeSet);
		 	unique_pair_set.insert(curNodeSet);
		 	CFrameSetSet_.insert(comp_set.begin() , comp_set.end() ); 
		 	sadj_map_.insert( std::make_pair(* CFvertex_biter ,curNodeSet));
		 	stotal_map_.insert( std::make_pair( * CFvertex_biter , comp_set));
		 
		}else/*�˽ڵ�����ڽӽڵ�*/
		{ 	 
		 	set<set<IndexType> > comp_set; //��ϵ����
		 	set<IndexType> cmps_node; //
		 	comp_set.clear();
		 	cmps_node.clear();
		 	cmps_node.insert(*CFvertex_biter);   //�����ڵ�
		 	comp_set.insert(cmps_node);
		 	cmps_node.clear();
		 
		 	set<set<IndexType> > ::iterator CFcomsetbitr ,CFcomseteitr;
		 	for( ; CFadj_bitr != CFadj_eitr ; ++CFadj_bitr )
			{
		 		cmps_node.insert( *CFadj_bitr);  // ���ڴ洢�ڵ�
		 		CFcomsetbitr = comp_set.begin();
		 		CFcomseteitr = comp_set.end();
		 		if( CFcomsetbitr == CFcomseteitr  )/*��ʼΪ��*/
				{
		 			//������ڵ�
		 			std::set<IndexType> cmps_node;
		 			cmps_node.insert(*CFvertex_biter);   
		 			comp_set.insert(cmps_node);
		 		}
		 		CFcomsetbitr = comp_set.begin();
		 		CFcomseteitr = comp_set.end();
		 
		 		set< set<IndexType> > tmpsetset;
		 		tmpsetset.clear();
		 		tmpsetset = comp_set;
		 		for( ; CFcomsetbitr != CFcomseteitr ; ++ CFcomsetbitr)
				{
		 
		 			set<IndexType> new_set = *CFcomsetbitr; 
		 			new_set.insert( *CFadj_bitr );
		 			tmpsetset.insert( new_set); 
		 		}
		 		comp_set = tmpsetset;
		 	}

		 	CFrameSetSet_.insert(comp_set.begin() , comp_set.end() );
		 	sadj_map_.insert( std::make_pair(* CFvertex_biter ,cmps_node));
		 	stotal_map_.insert( std::make_pair( * CFvertex_biter , comp_set));
		 	CFcomsetbitr = comp_set.begin();
		 	CFcomseteitr = comp_set.end();

		 	if(COUT_DEBUG)std::cout<<"��:"<< *CFvertex_biter <<"���������:" <<std::endl;
		 	for( ; CFcomsetbitr != CFcomseteitr ; ++ CFcomsetbitr)
			{
		 		//test
		 		set<IndexType> testset = * CFcomsetbitr;	 
		 		// ������һ֡��ѭ��
		 		curNodeSet = *CFcomsetbitr;

				if (curNodeSet.size() < 3) //��ֹ�������
				{
		 		    matchNextFrameSet(curNodeSet,*nextGraph,srLevel,tgLevel);
				}

		 		unique_pair_set.insert(*CFcomsetbitr);
		 	}		 
		 
		}
	}
}

void  GraphMatch::matchNextFrameSet(set<IndexType>& curNodeSet, LabelsGraph& tgGraph,IndexType srGraLevel,IndexType tgGraLevel)
{
	if(COUT_DEBUG)std::cout<<"������һ֡�ڵ㣬�������뵱ǰ�ڵ�"<<std::endl;

	LabelsGraph::vertex_iterator NFvertex_biter ;
	LabelsGraph::vertex_iterator NFvertex_eiter;
	LabelsGraph::adjacency_iterator  NFadj_bitr ,NFadj_eitr;

	set<IndexType> nextNodeSet;

	tie( NFvertex_biter ,NFvertex_eiter) = boost::vertices(tgGraph);

	for(; NFvertex_biter != NFvertex_eiter ; ++NFvertex_biter)
	{

	    boost::tie( NFadj_bitr , NFadj_eitr) = boost::adjacent_vertices( *NFvertex_biter,tgGraph); // ����ڽӵ�ĵ�����

		if( NFadj_bitr == NFadj_eitr )/* ˵���ýڵ�û���ڽӽڵ㣬�����ɿ�*/
		{
			std::set<set<IndexType> > comp_set;	
			nextNodeSet.clear();
			nextNodeSet.insert(*NFvertex_biter);
			comp_set.insert(nextNodeSet);
			NFrameSetSet_.insert(comp_set.begin() , comp_set.end() );
			tadj_map_.insert( std::make_pair(* NFvertex_biter ,nextNodeSet));
			ttotal_map_.insert( std::make_pair( * NFvertex_biter , comp_set));

			/*   �˴����� ��Ͻڼ� error*/
			//ScalarType  error = distance2PatchesLevel(curNodeSet,nextNodeSet,srGraLevel,tgGraLevel);

			ScalarType  error = distBetween2Frame4Items(curNodeSet,nextNodeSet,srGraLevel,tgGraLevel);

			equeue_.push( PatchMatch( curNodeSet, nextNodeSet, error ) );

		}else/*�˽ڵ�����ڽӽڵ�*/
		{
			std::set<set<IndexType> > comp_set;//ʹ�õ����˽ڵ��Ӧ�ĵ����ݣ���Ѱ�� ��һ֡�е�ƥ��ڵ�
			std::set<IndexType> cmps_node;
			comp_set.clear();
			cmps_node.clear();
			cmps_node.insert(*NFvertex_biter);   //
			comp_set.insert( cmps_node);
			cmps_node.clear();

			//std::stack<set<IndexType> > cmps_node_stack;  //����ջ�����н���ʽ���//�˴������ظ��½ڵ�
			

			std::set<set<IndexType> > ::iterator NFcomsetbitr ,NFcomseteitr;
			for(; NFadj_bitr != NFadj_eitr ; NFadj_bitr++ )
			{
		     	cmps_node.insert(*NFadj_bitr);
				NFcomsetbitr = comp_set.begin();
				NFcomseteitr = comp_set.end();
				if( NFcomsetbitr == NFcomseteitr)/*��ʼΪ��*/ 
				{
					//������ڵ�
					std::set<IndexType> cmps_node;
					cmps_node.insert(*NFvertex_biter);   
					comp_set.insert(cmps_node);
				}

				NFcomsetbitr = comp_set.begin();
				NFcomseteitr = comp_set.end();
				set< set<IndexType> > tmpsetset;
				tmpsetset.clear();
				tmpsetset = comp_set;

				for(; NFcomsetbitr != NFcomseteitr ; ++ NFcomsetbitr)
				{
					set<IndexType> new_set = *NFcomsetbitr; 
					new_set.insert( *NFadj_bitr );
					tmpsetset.insert( new_set); 

					if(COUT_DEBUG)std::cout<<"after tmpsetset test:"<< std::endl;
				}

				comp_set = tmpsetset;  //�ؼ�����

				if(COUT_DEBUG)std::cout<<"�ڵ����:"<< std::endl;
			}

			NFrameSetSet_.insert(comp_set.begin() , comp_set.end() );
			tadj_map_.insert( std::make_pair(* NFvertex_biter , cmps_node));
			ttotal_map_.insert( std::make_pair( * NFvertex_biter , comp_set));
			NFcomsetbitr = comp_set.begin();
			NFcomseteitr = comp_set.end();

			if(COUT_DEBUG)std::cout<<"��һ֡�ڵ�:"<< *NFvertex_biter <<"���������:" <<std::endl;

			for(; NFcomsetbitr != NFcomseteitr; ++ NFcomsetbitr)
			{
				nextNodeSet.clear();

				nextNodeSet = *NFcomsetbitr;

				if (nextNodeSet.size() < 3) //��ֹ�������
				{
				  /* ScalarType  error = distance2PatchesLevel(curNodeSet,nextNodeSet,srGraLevel,tgGraLevel);*/
				   ScalarType  error = distBetween2Frame4Items(curNodeSet,nextNodeSet,srGraLevel,tgGraLevel);
				   if (error < INF)
				   {
				       equeue_.push( PatchMatch(curNodeSet, nextNodeSet, error ) );
				   }

				}

			}//�������

		}//���ӵ�ǿ�

	}// ����ÿ���ڵ�
}

void GraphMatch::mergePatches(IndexType srLevel,IndexType tgLevel)
{

	//��ȡ��õļ������

	IndexType Item = 30;
	
	while (Item -- > 0)
	{
         PatchMatch bestPatch = equeue_.top();
		 equeue_.pop();

		 set<IndexType> srPatches = bestPatch.srPatches;
		 set<IndexType> tgPateches = bestPatch.tgPatches;

		 Logger<<Item<<"�����ŵ����.\n";

		 if (srPatches.size() > 0)
		 {
			 Logger<<" Merge  srFrame  patches.\n";
			 for (auto  iter = srPatches.begin(); iter != srPatches.end(); ++ iter)
			 {
				 Logger<<(*iter)<<"  ";
			 }
			 Logger<<".\n";

		    // mergeSrPatches(srLevel,srPatches);
		 }

		 if (tgPateches.size() > 0)
		 {
			 Logger<<" Merge  tgFrame  patches.\n";
			 for (auto  iter = tgPateches.begin(); iter != tgPateches.end(); ++ iter)
			 {
				 Logger<<(*iter)<<"  ";
			 }
			 Logger<<".\n";

		     //mergeTgPatches(tgLevel,tgPateches);
		 }
	}
}

void GraphMatch::mergeSrPatches(IndexType srLevel,set<IndexType>& srBestCombine)
{

	mergePatchesOri(srLevel,srFrame,srBestCombine);

}

void GraphMatch::mergeTgPatches(IndexType tgLevel,set<IndexType>& tgBestCombine)
{
	
	mergePatchesOri(tgLevel,tgFrame,tgBestCombine);
}

//�����ĺϲ��㼯
IndexType GraphMatch::mergePatchesOri(IndexType depth,IndexType frameId, set<IndexType>& patches)
{

	if ( patches.size() < 2)
	{
		return (*patches.begin());
	}

	IndexType firstId  = (*patches.begin() );
	IndexType secondId = ( *(++patches.begin() ) );


	//assert �����ڵ�����

	HFrame& cframe = components_[frameId];

	//IndexType depth_ = 0;//ÿ�ζ���ȡ��ͬ�Ĳ�
// 	depth_ = cframe.hier_graph.size();
// 	--depth_;

	vector<LabelsGraph*>& chierGraph = cframe.hier_graph; //���е�ͼ��

    auto test = boost::edge(firstId,secondId,*chierGraph[depth]);

	if (!test.second)
	{
		Logger<<firstId<<" ��"<<secondId<<"������.\n";
		
		return firstId;
	}


	map<IndexType ,IndexType>& labelofvtx = cframe.label_of_vtx;   //��ʼ labelofvtx

	LabelsGraph* p =new LabelsGraph( *( chierGraph[depth]) );  //������һ���µ�         //��ʼ����ͼ
	LabelsGraph& oriGra = *p;//��������oriGra

	vector< vector<HLabel*> >& chierlabucket =  cframe.hier_label_bucket;   //���е�label_bucket

	assert( chierlabucket.size() > depth);

	vector<HLabel*>* p_t;

	p_t = DeepCopyHLableVec( &chierlabucket[depth] );    //��������

	//copyLabelBucket(p_t,chierlabucket[depth]);

	vector<HLabel*>& old_label_bucket = *p_t;//��������old label bucket

	//Sample& csmp = SampleSet::get_instance()[frameId];

	LabelsGraph* new_graph  ;  //��ͼ ָ��
	IndexType	oriNodeNum  = boost::num_vertices( oriGra);     
	IndexType  targetLabel;
	map<IndexType,IndexType> labelIdmap;
	new_graph = mergeIter(depth,frameId,oriGra , old_label_bucket ,patches ,targetLabel,labelIdmap);
	components_[frameId].hier_graph.push_back(new_graph);   //������һ�㱣���½��
	components_[frameId].hier_label_vtxBucket_index.push_back(labelIdmap);
	chierlabucket.push_back(old_label_bucket);
	Logger<<" label bucket ��С��"<< chierlabucket.size()<<std::endl;

	return targetLabel;

}

IndexType GraphMatch::mergeAllNodes(IndexType frameId, set<IndexType>& patches)
{
	if ( patches.size() < 2)
	{
		return (*patches.begin());
	}

	HFrame& cframe = components_[frameId];

	IndexType depth = 0;//ÿ�ζ���ȡ��ͬ�Ĳ�
	 depth = cframe.hier_graph.size();
	 -- depth;

	vector<LabelsGraph*>& chierGraph = cframe.hier_graph;

	map<IndexType ,IndexType>& labelofvtx = cframe.label_of_vtx;   //��ʼ labelofvtx

	LabelsGraph* p =new LabelsGraph( *( chierGraph[depth]) );  //������һ���µ�         //��ʼ����ͼ
	LabelsGraph& oriGra = *p;//��������oriGra

	vector< vector<HLabel*> >& chierlabucket =  cframe.hier_label_bucket;   //

	assert( chierlabucket.size() > depth);

	vector<HLabel*>* p_t;

	p_t = DeepCopyHLableVec( &chierlabucket[depth] );    //��������

	//copyLabelBucket(p_t,chierlabucket[depth]);

	vector<HLabel*>& old_label_bucket = *p_t;//��������old label bucket

	Sample& csmp = SampleSet::get_instance()[frameId];

	LabelsGraph* new_graph  ;  //��ͼ ָ��
	IndexType	oriNodeNum  = boost::num_vertices( oriGra);     
	IndexType  targetLabel;
	map<IndexType,IndexType> LabelIdIndex;
	new_graph = mergeIter(depth,frameId,oriGra , old_label_bucket ,patches ,targetLabel,LabelIdIndex);
	components_[frameId].hier_graph.push_back(new_graph);   //������һ�㱣���½��
	components_[frameId].hier_label_vtxBucket_index.push_back(LabelIdIndex);
	chierlabucket.push_back(old_label_bucket);

	Logger<<" label bucket ��С��"<< chierlabucket.size()<<std::endl;

	return targetLabel;
}

LabelsGraph* GraphMatch::mergeIter(IndexType gLevel,  IndexType frameId, LabelsGraph& _orIG ,
								   vector<HLabel*>& old_label_bucket, set<IndexType>& patches ,IndexType& new_label,map<IndexType,IndexType>& LabelIdIndex )
{
	if( patches.size() < 2)
	{
		set<IndexType>::iterator itr = patches.begin();
		new_label = *itr ;
		return &_orIG;

	}

	set<IndexType>* new_patches ;
	new_patches = &patches;
	LabelsGraph& oriGra = _orIG;   //��ʼ����ͼ

	//�ȱ���ÿ���ڵ� ��������е������
	IndexType	oriNodeNum ;
	IndexType	newNodeNum;
	newNodeNum = oriNodeNum  = boost::num_vertices( oriGra);   //��ʼ��ֵ
	assert( newNodeNum > patches.size());

	LabelsGraph*	newGraph = &_orIG;  // ��ʼ��ֵ
	VertexIterator vtxbitr ,vtxeitr;
	boost::tie(vtxbitr , vtxeitr) = boost::vertices( oriGra);

	VertexDescriptor cvd;

	for(; vtxbitr != vtxeitr; ++ vtxbitr)
	{
		HLabel& cLable = *(old_label_bucket[ *vtxbitr]); //ȥ����û����
	    

		new_patches = new set<int>();
		new_patches->clear();
		cvd = *vtxbitr;
		IndexType clabelSize = cLable.vertex_bucket.size();

		//Logger<<"label �ţ� "<< cLable.label_id <<"   ��ĸ��� :"<< clabelSize<<std::endl;

		if( patches.find( cvd) != patches.end() )
		{
			newGraph = new LabelsGraph();
			OutEdgeIterator outegbitr , outegeitr ;
			boost::tie( outegbitr , outegeitr) = boost::out_edges( *vtxbitr , oriGra);
			map<IndexType , map<IndexType , HVertex*> > recordColapseEdges;

			vector< VertexDescriptor> recordEdgeVD;

			set<GraphEdgeProperty> collapseEdgesProperty;
			OutEdgeIterator nextoutegbitr;  //���ڱ�����һ����������ַ������ɾ��֮�������ʧЧ

			set<IndexType> lowAdjNodeId;

			// ��������� �ڽڵ�
			IndexType bestNeibIndexInVector = 0 ;

			for(  nextoutegbitr = outegbitr ; outegbitr != outegeitr ; outegbitr = nextoutegbitr)
			{
				++nextoutegbitr;
				VertexDescriptor sd = boost::target( *outegbitr ,oriGra);

				lowAdjNodeId.insert(sd);

				if (patches.find(sd) != patches.end())
				{
					bestNeibIndexInVector = sd;
				}

				recordEdgeVD.push_back( sd);
				EdgeDescriptor& tobeRemovedEdge = *outegbitr;
				GraphEdgeProperty& nextEP = ( oriGra)[tobeRemovedEdge];
				collapseEdgesProperty.insert( nextEP );
				boost::remove_edge( tobeRemovedEdge , oriGra);   //ɾ��������
			}


			//���½���;;��ڵ㵽bestNode
			IndexType eId = 0;


			for (auto vIter = lowAdjNodeId.begin(); vIter != lowAdjNodeId.end(); ++vIter)
			{
				if ((*vIter) == bestNeibIndexInVector ) continue;

				GraphEdgeProperty edgeP;
				edgeP.index = eId;

				if ( (*vIter) < bestNeibIndexInVector)
				{
					edgeP.start_ = (*vIter);
					edgeP.end_ = bestNeibIndexInVector;

				}else
				{
					edgeP.start_ = bestNeibIndexInVector;
					edgeP.end_ = (*vIter);
				}

				boost::add_edge( edgeP.start_, edgeP.end_, edgeP, oriGra );

			}

// 			for( IndexType i = 0 ; i <recordNeibNode.size() ; ++ i)
// 			{
// 				if( i == bestNeibIndexInVector){ continue;}
// 				//��Ӧ�ü���property
// 				IndexType srLabelId = recordNeibNode[bestNeibIndexInVector]->label_id;
// 
// 				IndexType tgLabelId = recordNeibNode[i]->label_id;
// 
// 				IndexType srNodeId =  components_[frameId].hier_label_vtxBucket_index[gLevel][srLabelId];
// 
// 				IndexType tgNodeId =  components_[frameId].hier_label_vtxBucket_index[gLevel][tgLabelId];
// 
// 				GraphEdgeProperty edgeP;
// 				edgeP.index = eId;
// 				edgeP.start_ = srNodeId;
// 				edgeP.end_ = tgNodeId;
// 
// 				boost::add_edge( srNodeId, tgNodeId, edgeP, oriGra );
// 
// 				//boost::add_edge( recordNeibNode[bestNeibIndexInVector]->label_id , recordNeibNode[i]->label_id ,oriGra );
// 
// 			}

			HLabel& bestNeiLabel = *(old_label_bucket[ bestNeibIndexInVector] );

			mergeTwoLabelAndPushBack(cLable , bestNeiLabel);

			// old_label_bucket�Ѿ����µ��� ����������.

			processAndCreateNewGraphAndPushBack(oriGra , cvd ,old_label_bucket , components_[frameId].label_of_vtx ,*newGraph); //ϴ��,��֤���������
			

			//�������ɵ�labelId �ͽڵ����map���¸�ֵ
			//map<IndexType,IndexType> updateIndexMap;// = components_[frameId].hier_label_vtxBucket_index[gLevel];
			IndexType gLevelIndex = 0;
			for (auto iter = old_label_bucket.begin(); iter != old_label_bucket.end(); ++ iter,++ gLevelIndex)
			{
				IndexType labelId = (*iter)->label_id;
				LabelIdIndex[labelId] = gLevelIndex;
			}


			if( newGraph !=&oriGra)
			{
				delete &oriGra;
			} 

			newNodeNum = oriNodeNum -1;
			break;
		}


	}//���������

	set<IndexType>::iterator setbitr ,seteitr;
	setbitr = patches.begin();
	seteitr = patches.end();
	for( ; setbitr != seteitr ; ++setbitr)
	{
		if(*setbitr == cvd) continue;

		if(*setbitr< cvd)
		{
			(*new_patches).insert( *setbitr);
		}
		else
		{
			(*new_patches).insert( (*setbitr) -1);
		}
	}

	return newGraph;
	//return mergeIter(gLevel,frameId, *newGraph ,old_label_bucket, *new_patches ,new_label,LabelIdIndex);

}

IndexType GraphMatch::mergeAllSetOneTime(LabelsGraph& _orIG ,vector<HLabel*>& old_label_bucket , set<IndexType>& patches ,IndexType& new_label )
{

	return 0;
}

vector<HLabel*>* GraphMatch::DeepCopyHLableVec( const vector<HLabel*>* _p_oriHLabelvec )
{
	vector<HLabel*>* _p_newHLabelvec;
	_p_newHLabelvec = new vector<HLabel*>(*_p_oriHLabelvec);
	vector<HLabel*>::iterator bitr ,eitr;
	bitr = _p_newHLabelvec->begin();
	eitr = _p_newHLabelvec->end();
	for(; bitr!= eitr ; ++bitr)
	{
		HLabel* p_orihlabel = * bitr;
		HLabel* p_newhlabel;
		p_newhlabel = DeepCopyHLabel(p_orihlabel);
		*bitr = p_newhlabel;
	}

	return _p_newHLabelvec;
}

HLabel* GraphMatch::DeepCopyHLabel( const HLabel* _p_oriHlabel  )
{
	HLabel* _p_newHlabel = new HLabel( *_p_oriHlabel);
	//��� newHlabel �е� label_id �� frame_parent  ����Ҫ�䣬 prev_corr ��next_corr Ҳ����Ҫ�� �� ��vertex_bucket �е�Hvertex�豻����
	//map<IndexType , HVertex*> new_vertex_bucket;
	map<IndexType ,HVertex*>::iterator  hvtxbitr ,hvtxeitr;
	hvtxbitr = (*_p_newHlabel).vertex_bucket.begin();
	hvtxeitr = (*_p_newHlabel).vertex_bucket.end();
	for( ; hvtxbitr != hvtxeitr ;++hvtxbitr)
	{
		HVertex* p_orivtx = hvtxbitr->second;
		HVertex* p_newvtx;
		p_newvtx = DeepCopyVertex(p_orivtx  ,_p_newHlabel);
		hvtxbitr->second = p_newvtx;
	}
	return _p_newHlabel;

}
HVertex* GraphMatch::DeepCopyVertex( const HVertex* _p_orivtx, HLabel* new_label_parent )
{
	HVertex* _p_newvtx = new HVertex( *_p_orivtx);
	//����µ� HVertex �е� label_parent �����¸�ֵ��

	IndexType gSize = _p_newvtx->label_parent.size();

	_p_newvtx->label_parent[gSize - 1] = new_label_parent;

	return _p_newvtx;
}

void GraphMatch::mergeTwoLabelAndPushBack(HLabel& _cLabel , HLabel& _bestNeiLabel_  )
{
// #define mergeTwoLabel_tset 0
//  	// �� _cLabel_ �� ����bestNeiLabel �ϲ������Ѻϲ��Ľ������ ����clabel��
//  	map<IndexType , HVertex*>& lvtxbck = _cLabel.vertex_bucket;
//  	map<IndexType , HVertex*>& rvtxbck = _bestNeiLabel_.vertex_bucket;
//  
//  	//_cLabel_.label_id  =  10000;  //��ʾ����Ҫ������ϴ��
//  	IndexType	mergedLabel = _bestNeiLabel_.label_id;   //�ϲ���ʹ�õ�label��
//  	map<IndexType , HVertex*>::iterator l_vtxbkbitr ,l_vtxbkeitr;
//  	map<IndexType , HVertex*>::iterator r_vtxbkbitr ,r_vtxbkeitr;
//  	l_vtxbkbitr = lvtxbck.begin();
//  	l_vtxbkeitr = lvtxbck.end();
//  	for( l_vtxbkbitr = lvtxbck.begin() ; l_vtxbkbitr != l_vtxbkeitr ;++ l_vtxbkbitr)
//  	{
//  		HVertex&  vtx = *(l_vtxbkbitr->second) ;
//  		IndexType gSize = vtx.label_parent.size();
//  
//  		vtx.label_parent[gSize - 1] = &_bestNeiLabel_;
//  	}
//  
//  	rvtxbck.insert( lvtxbck.begin() , lvtxbck.end()); //�������ݿ�ʼת��
//  
//  	//���ԭ���ĵ㼯
//  	_cLabel.vertex_bucket.clear();//0901
//  
//  
//  	//      free_labelId_ = _bestNeiLabel_.label_id;
//  	// �� frame �е� label_of_vtx  ���д���
//  	map<IndexType ,IndexType>& labelofvtx= _bestNeiLabel_.frame_parent->label_of_vtx;
//  
//  	r_vtxbkbitr = rvtxbck.begin();
//  	r_vtxbkeitr = rvtxbck.end();
//  	for( ; r_vtxbkbitr!= r_vtxbkeitr ;++r_vtxbkbitr)
//  	{
//  		labelofvtx[ r_vtxbkbitr->first] = mergedLabel;
//  	}



 	///ѡ���ö���������Ϊ���յ�label����0901---yuanqing----
 	// �� _cLabel_ �� ����bestNeiLabel �ϲ������Ѻϲ��Ľ������ ����clabel��
 	map<IndexType , HVertex*>& lvtxbck = _cLabel.vertex_bucket;
 	map<IndexType , HVertex*>& rvtxbck = _bestNeiLabel_.vertex_bucket;
 
 	IndexType lsize = lvtxbck.size();
 	IndexType rsize = rvtxbck.size();
 
 	IndexType	mergedLabel = 0;
 	if (lsize > rsize)
 	{
 		mergedLabel = _cLabel.label_id;
 		_bestNeiLabel_.label_id = mergedLabel;
 
 	}else
 	{
 		mergedLabel = _bestNeiLabel_.label_id;
 	}
 
 	//IndexType	mergedLabel = _bestNeiLabel_.label_id;   //�ϲ���ʹ�õ�label��
 	map<IndexType , HVertex*>::iterator l_vtxbkbitr ,l_vtxbkeitr;
 	map<IndexType , HVertex*>::iterator r_vtxbkbitr ,r_vtxbkeitr;
 	l_vtxbkbitr = lvtxbck.begin();
 	l_vtxbkeitr = lvtxbck.end();
 	for( l_vtxbkbitr = lvtxbck.begin() ; l_vtxbkbitr != l_vtxbkeitr ;++ l_vtxbkbitr)
 	{
 		HVertex&  vtx = *(l_vtxbkbitr->second) ;
 		IndexType gSize = vtx.label_parent.size();
 
 		vtx.label_parent[gSize - 1] = &_bestNeiLabel_;
 	}
 
 	rvtxbck.insert( lvtxbck.begin() , lvtxbck.end()); //�������ݿ�ʼת��
 
 	//���ԭ���ĵ㼯
 	_cLabel.vertex_bucket.clear();//0901
 
 
 	//      free_labelId_ = _bestNeiLabel_.label_id;
 	// �� frame �е� label_of_vtx  ���д���
 	map<IndexType ,IndexType>& labelofvtx= _bestNeiLabel_.frame_parent->label_of_vtx;
 
 	r_vtxbkbitr = rvtxbck.begin();
 	r_vtxbkeitr = rvtxbck.end();
 	for( ; r_vtxbkbitr!= r_vtxbkeitr ;++r_vtxbkbitr)
 	{
 		labelofvtx[ r_vtxbkbitr->first] = mergedLabel;
 	}

}

void GraphMatch::processAndCreateNewGraphAndPushBack(LabelsGraph& _preGraph , VertexDescriptor& _targetVd ,vector<HLabel*>& _labelvec_ ,
										 map<IndexType ,IndexType>& labelofvtx_, LabelsGraph& newGraph_)
{

#define processAndCreateNewGraph_test 0

	//assert(_targetVd > =0 );
	assert( &_preGraph != &newGraph_);
	VertexIterator vtxbitr ,vtxeitr;
	vector<HLabel*> newLabelvec;
	boost::tie( vtxbitr ,vtxeitr) = boost::vertices( _preGraph);
	for( ; vtxbitr != vtxeitr ; ++vtxbitr)
	{
		if( *vtxbitr < _targetVd)
		{
			GraphVertexProperty newvp = _preGraph[*vtxbitr];
			add_vertex( newvp ,newGraph_);
			newLabelvec.push_back( _labelvec_[ *vtxbitr]);
		}
		else if( *vtxbitr == _targetVd)
		{


		}else if( *vtxbitr >_targetVd)
		{
			GraphVertexProperty newvp = _preGraph[*vtxbitr];
			newvp.index =  *vtxbitr-1;
			add_vertex( newvp ,newGraph_);

			//��Ҫ��label ��Ӧ��ʵ�崦��
			//_labelvec_[ *vtxbitr]->label_id = *vtxbitr-1;

			map<IndexType , HVertex*>::iterator vtxbckbitr = _labelvec_[ *vtxbitr]->vertex_bucket.begin();
			map<IndexType , HVertex*>::iterator vtxbckeitr = _labelvec_[ *vtxbitr]->vertex_bucket.end();
			for( ;vtxbckbitr != vtxbckeitr ;++vtxbckbitr)
			{
				//Logger<<"vtxbucket->first"<<vtxbckbitr->first<<"*vtxbitr"<<(*vtxbitr)-1<<std::endl;
				labelofvtx_[ vtxbckbitr->first] = (*vtxbitr) -1;
			}
			//Logger<<"labelvec_ size"<<_labelvec_.size()<<std::endl;
			newLabelvec.push_back( _labelvec_[ *vtxbitr]);

		}else
		{


		}

	}//���ű�����

	_labelvec_.clear();
	_labelvec_ = newLabelvec;

	IndexType egId = 0;
	GraphEdgeProperty newEdge;

	EdgeIterator	egbitr ,egeitr;
	boost::tie( egbitr ,egeitr) = boost::edges( _preGraph);
	for( ; egbitr != egeitr ; ++egbitr)
	{
		EdgeDescriptor		ed = *egbitr;
		VertexDescriptor	sd = boost::source( ed ,_preGraph);
		VertexDescriptor	td = boost::target( ed ,_preGraph);
		if( _targetVd == sd || _targetVd == td )
		{
			continue;
		}
		if( (sd<_targetVd) &&( td< _targetVd))
		{
			newEdge.index = egId;
			if (sd < td)
			{
			   newEdge.start_ = sd;
			   newEdge.end_ = td;
			}else
			{
				newEdge.start_ = td;
				newEdge.end_ = sd;
			}

			boost::add_edge(sd, td, newEdge, newGraph_);

			++ egId;

		}else if( (sd>_targetVd) &&( td< _targetVd))
		{
			newEdge.index = egId;

			if (sd - 1 < td)
			{
				newEdge.start_ = sd - 1;
				newEdge.end_ = td;
			}else
			{
				newEdge.start_ = td;
				newEdge.end_ = sd - 1;
			}

			boost::add_edge( sd-1 ,td , newEdge, newGraph_);

			++egId;

		}else if( (sd<_targetVd) &&( td> _targetVd))
		{
			newEdge.index = egId;

			if (sd < td - 1)
			{
				newEdge.start_ = sd;
				newEdge.end_ = td - 1;
			}else
			{
				newEdge.start_ = td -1;
				newEdge.end_ = sd;
			}

			boost::add_edge( sd , td-1 ,newEdge,newGraph_);

			++egId;

		}else if((sd>_targetVd) &&( td> _targetVd))
		{
			newEdge.index = egId;

			if (sd -1 < td - 1)
			{
				newEdge.start_ = sd - 1;
				newEdge.end_ = td - 1;
			}else
			{
				newEdge.start_ = td -1;
				newEdge.end_ = sd - 1;
			}
			boost::add_edge( sd -1 ,td -1 ,newEdge,newGraph_);

			++egId;
		}

	}


// 	if(processAndCreateNewGraph_test)
// 	{
// 		coutGraph(newGraph_);
// 	}
}

void GraphMatch::copyLabelBucket(vector<HLabel*>& leftLabels, const vector<HLabel*>& oriLabels)
{
	assert(oriLabels.size() > 0);

	leftLabels.resize(0);

	for (auto iter = oriLabels.begin(); iter != oriLabels.end(); iter ++)
	{
		IndexType LId = (*iter)->label_id;
		HLabel* new_label_space = allocator_.allocate<HLabel>();
		HLabel* new_label = new (new_label_space)HLabel((*iter)->label_id,(*iter)->frame_parent,(*iter)->vertex_bucket,(*iter)->prev_corr,(*iter)->next_corr);
		leftLabels.push_back(new_label);
	}
}

void GraphMatch::buildPatchCorrespondenceByLabel()
{
	for (auto fiter = components_.begin(); fiter != components_.end(); ++fiter )
	{
		 IndexType frame_id = fiter->first;
		 vector<HLabel*>& labelBucket = fiter->second.hier_label_bucket[0];
		 
		if (components_.find(frame_id + 1) != components_.end() )
		{

			map<IndexType,IndexType> sr_labelIndexMap = fiter->second.hier_label_vtxBucket_index[0];
			map<IndexType,IndexType> tg_labelIndexMap = components_[frame_id + 1].hier_label_vtxBucket_index[0];

			vector<HLabel*>& tg_labelBucket = components_[frame_id + 1].hier_label_bucket[0];

			IndexType i = 0;
			for ( auto mIter = sr_labelIndexMap.begin(); mIter != sr_labelIndexMap.end(); ++mIter,++i )
			{
				IndexType labelId = mIter->first;
				if (tg_labelIndexMap.find(labelId) != tg_labelIndexMap.end() )//������ͬ��label
				{
					IndexType corLabelId = tg_labelIndexMap[labelId];

					labelBucket[i]->next_corr = tg_labelBucket[corLabelId]; //next correspondence

					tg_labelBucket[corLabelId]->prev_corr = labelBucket[i];
				}
			}

		}
	}

	//printMatchingInformation();

}

void GraphMatch::printMatchingInformation()
{
	char* filename = "Correspondence Infor";

	FILE* outfile = fopen(filename, "w");

	map<IndexType,bool> isTrav;

	for (auto fiter = components_.begin(); fiter != components_.end(); ++fiter)
	{
		IndexType frame_id = fiter->first;
		vector<HLabel*>& labelBucket = fiter->second.hier_label_bucket[0];

		if (components_.find(frame_id + 1) != components_.end() )
		{
			IndexType lId = 0;
			for (auto lIter = labelBucket.begin(); lIter != labelBucket.end(); ++lIter, ++ lId )
			{
				IndexType labeId = (*lIter)->label_id;

				IndexType vtxNum = (*lIter)->vertex_bucket.size();

				HLabel* nextLabelBucket = (*lIter)->next_corr;

				IndexType flKey = frame_label_to_key(frame_id,labeId);
				
				if(!isTrav[flKey])
				{
					isTrav[flKey] = true;

					fprintf(outfile,"Frame Id %d have %d  points ,The label Id is %d\n", frame_id, vtxNum, labeId);

					//Logger<<" Frame"<<frame_id<<" have "<<vtxNum<<" points, and the Label Id is "<<labeId<<endl;
					if (nextLabelBucket != NULL)
					{
						while (nextLabelBucket != NULL)
						{
							IndexType nFId = nextLabelBucket->frame_parent->frame_id;

							IndexType nLId = nextLabelBucket->label_id;

							IndexType  nVtxNum = nextLabelBucket->vertex_bucket.size();

							IndexType nextFlKey = frame_label_to_key(nFId,nLId);
							isTrav[nextFlKey] = true;

							fprintf(outfile,"Correctness Next patch information is: \n");
							fprintf(outfile,"Frame Id %d have %d  points ,The label Id is %d\n", nFId, nVtxNum, nLId);
							//Logger<<" Correctness Next patch information is: "<<endl;
							//Logger<<" Frame id is "<<nFId<<" label Id is "<<nLId<<" and have "<<nVtxNum<<"points"<<endl;

							nextLabelBucket = nextLabelBucket->next_corr;

						}
					}else
					{
						fprintf(outfile," This patch has not next correspondence!!!!!!.\n");
						//Logger<<" This patch has not next correspondence!!!!!!.\n";
					}

				}else
				{
					fprintf(outfile," This patch was trav befor!.\n");
					//Logger<<" This patch was trav befor!.\n";
				}

			}
		}
	}

	fclose(outfile);
}

void GraphMatch::findBestPatches(IndexType srLevel, IndexType tgLevel,set<IndexType>& srBestSet, set<IndexType>& tgBestSet)
{
	set<IndexType> srTempBest;
	set<IndexType> tgTempBest;
	srTempBest.clear();
	tgTempBest.clear();

	if (!equeue_.empty())
	{
		equeue_.pop();
	}

	//����ԭͼ�ĵ㼯--����Ŀ��ͼ�ĵ㼯
	HFrame& cur_frame = components_[srFrame];
	HFrame& next_frame = components_[tgFrame];
	LabelsGraph* curGraph = cur_frame.hier_graph[srLevel];
	LabelsGraph* nextGraph = next_frame.hier_graph[tgLevel];

	VertexIterator vBeginIter,vEndIter;
	AdjacencyIterator vAdjBeginIter,vAdjEndIter;
	
	tie(vBeginIter,vEndIter) = boost::vertices(*curGraph);
	
	vector<HLabel*>& vtxBucket = components_[srFrame].hier_label_bucket[srLevel];

	for (; vBeginIter != vEndIter; ++ vBeginIter)
	{
		srTempBest.clear();
		tgTempBest.clear();

		//�����ڽӱ߿ղ���,�����㶼Ҫ�����ĵ����ж�
			srTempBest.insert(*vBeginIter);
             
			//�ҵ�����ڵ�����һ֡�Ķ�Ӧ��
			HLabel* nextPatch = vtxBucket[*vBeginIter]->next_corr;

			if (nextPatch != NULL)//��Ӧ��ǿ�
			{
				IndexType labelId = nextPatch->label_id;

				vector<HLabel*>& vtxNextBucket = components_[tgFrame].hier_label_bucket[tgLevel];

				IndexType nodeId = components_[tgFrame].hier_label_vtxBucket_index[tgLevel][labelId];

				pair<VertexIterator, VertexIterator> vi = boost::vertices(*nextGraph);
				VertexIterator nodeIter = (vi.first + nodeId);
				AdjacencyIterator vNextBegIter,vNextEndIter;
				boost::tie(vNextBegIter,vNextEndIter) = boost::adjacent_vertices(*nodeIter,*nextGraph);

				if (vNextBegIter == vNextEndIter)//��Ӧ��Ҳ�ǹ�����
				{
					tgTempBest.insert(*nodeIter);

					ScalarType error = distBetween2Setpatches(srTempBest,tgTempBest,srLevel,tgLevel);

					if (error < INF)
					{
						equeue_.push( PatchMatch(srTempBest,tgTempBest, error ) );
					}

					tgTempBest.erase(*nodeIter);

				}else
				{
					tgTempBest.insert(*nodeIter);

					IndexType vtxOriSize = vtxNextBucket[nodeId]->vertex_bucket.size();

					for (; vNextBegIter != vNextEndIter; ++ vNextBegIter)//���ʶ�Ӧ����ڱ�
					{
						IndexType neigVtxSize = vtxNextBucket[*vNextBegIter]->vertex_bucket.size();

// 						if (neigVtxSize < vtxOriSize)
// 						{
							tgTempBest.insert(*vNextBegIter);

							ScalarType error = distBetween2Setpatches(srTempBest,tgTempBest,srLevel,tgLevel);

							if (error < INF)
							{
								equeue_.push( PatchMatch(srTempBest,tgTempBest, error ) );
							}

							tgTempBest.erase(*vNextBegIter);

// 						}else
// 						{
// 							Logger<<" ��Ŀ��֡��,��ֹС���̲����!,�������!.\n";
// 						}
					}
				}

			}else
			{
				Logger<<"�ù�����û��next, ��ʱ������.\n";
			}//��Ӧ��Ϊ��

		//�����ڽӱ߿ղ���,�����㶼Ҫ�����ĵ����ж�


        //���ڽӱ�
        tgTempBest.clear();

		boost::tie(vAdjBeginIter,vAdjEndIter) = boost::adjacent_vertices(*vBeginIter,*curGraph);

		if (vAdjBeginIter != vAdjEndIter)//�ǹ����Ľڵ�
		{
           //srTempBest.insert(*vBeginIter);	   
		   IndexType vtxNumSrOrignal = vtxBucket[*vBeginIter]->vertex_bucket.size();

		   for (;vAdjBeginIter != vAdjEndIter; ++vAdjBeginIter)
		   {
			   IndexType vtxNumSrOriAdj = vtxBucket[*vAdjBeginIter]->vertex_bucket.size();

			   if (vtxNumSrOriAdj < vtxNumSrOrignal)
			   {
				   srTempBest.insert(*vAdjBeginIter);//��ʱ,srBestSet�����������ڵ�ԭʼ
				   //�ҵ����ԵĶ�Ӧ A,B,tgTempSet�ܹ��������{A} {B} {A,B}
				   //case1  ���vBeginIter�Ķ�Ӧ��
				   //�ҵ�����ڵ�����һ֡�Ķ�Ӧ��
				   HLabel* nextOriPatch = vtxBucket[*vBeginIter]->next_corr;
				   HLabel* nextAdjPatch = vtxBucket[*vAdjBeginIter]->next_corr; 

				   IndexType nextOriId = -1;
				   IndexType nextAdjId = -1;

				   if (nextOriPatch)   //{A}
				   {
					   IndexType labelId = nextOriPatch->label_id;

					   vector<HLabel*>& vtxNextBucket = components_[tgFrame].hier_label_bucket[tgLevel];

					   IndexType nodeId = components_[tgFrame].hier_label_vtxBucket_index[tgLevel][labelId];

					   tgTempBest.insert(nodeId);

					   nextOriId = nodeId;

					   ScalarType error = distBetween2Setpatches(srTempBest,tgTempBest,srLevel,tgLevel);

					   if (error < INF)
					   {
						   equeue_.push( PatchMatch(srTempBest,tgTempBest, error ) );
					   }

					   tgTempBest.erase(nodeId);
				   }

				   if (nextAdjPatch)  //{B}
				   {
					   IndexType labelId = nextAdjPatch->label_id;

					   vector<HLabel*>& vtxNextBucket = components_[tgFrame].hier_label_bucket[tgLevel];

					   IndexType nodeId = components_[tgFrame].hier_label_vtxBucket_index[tgLevel][labelId];

					   tgTempBest.insert(nodeId);

					   nextAdjId = nodeId;

					   ScalarType error = distBetween2Setpatches(srTempBest,tgTempBest,srLevel,tgLevel);

					   if (error < INF)
					   {
						   equeue_.push( PatchMatch(srTempBest,tgTempBest, error ) );
					   }

					   tgTempBest.erase(nodeId);
				   }

				   if (nextOriPatch && nextAdjPatch)  //{A,B}, A��B �п��ܲ�����.
				   {
					   tgTempBest.insert(nextOriId);
					   tgTempBest.insert(nextAdjId);

					   ScalarType error = distBetween2Setpatches(srTempBest,tgTempBest,srLevel,tgLevel);

					   if (error < INF)
					   {
						   equeue_.push( PatchMatch(srTempBest,tgTempBest, error ) );
					   }

					   tgTempBest.erase(nextOriId);
					   tgTempBest.erase(nextAdjId);

				   }

				   srTempBest.erase(*vAdjBeginIter);

			   }else
			   {
				   Logger<<" ��ԭʼ֡��, ��ֹС���̲����!,�������!.\n";
			   }
		   }


		}//�ýڵ������ӱ�

	}//����ԭͼ�����е�

	//
	//printNiceMatchingPatches();

	while (equeue_.size() > 0)
	{
		PatchMatch bestPatch = equeue_.top();

		if (bestPatch.srPatches.size() > 1 || bestPatch.tgPatches.size() > 1)
		{
			srBestSet = bestPatch.srPatches;
			tgBestSet = bestPatch.tgPatches;
			break;
		}

		equeue_.pop();
	}

// 	PatchMatch bestPatch = equeue_.top();
// 	srBestSet = bestPatch.srPatches;
// 	tgBestSet = bestPatch.tgPatches;

}

void GraphMatch::mergePatchesAfterCoSeg(IndexType srLevel, IndexType tgLevel,set<IndexType>& srBestSet, set<IndexType>& tgBestSet)
{
	IndexType fpId = 1;
	IndexType spid = 6;

	printf("Start to merge patch %d and %d.\n",fpId,spid);

	//�㼯Ϊһ��Ԫ���򲻺ϲ�,�����㲻����Ҳ���ϲ�..
	auto  isExistF_sr = components_[srFrame].hier_label_vtxBucket_index[srLevel].find(fpId);
	auto  isExistS_sr = components_[srFrame].hier_label_vtxBucket_index[srLevel].find(spid);

	map<IndexType,IndexType>::iterator idxEnd_sr = components_[srFrame].hier_label_vtxBucket_index[srLevel].end();

	if (isExistF_sr == idxEnd_sr || isExistS_sr == idxEnd_sr)
	{
		Logger<<"���������еĿ�.\n";
	}else
	{
		IndexType fnid_sr = components_[srFrame].hier_label_vtxBucket_index[srLevel][fpId];
		IndexType snid_sr = components_[srFrame].hier_label_vtxBucket_index[srLevel][spid];

		srBestSet.insert(fnid_sr);
		srBestSet.insert(snid_sr);

		if (srBestSet.size() > 1)
		{
			Logger<<" Start Merge source Graph  "<<srFrame<<endl;

			mergeSrPatches(srLevel,srBestSet);

			Logger<<" End Merge source Graph.\n";
		}
	}




	//merge patches 
 	auto  isExistF = components_[tgFrame].hier_label_vtxBucket_index[tgLevel].find(fpId);
 	auto  isExistS = components_[tgFrame].hier_label_vtxBucket_index[tgLevel].find(spid);
 
     map<IndexType,IndexType>::iterator idxEnd = components_[tgFrame].hier_label_vtxBucket_index[tgLevel].end();
 
 	if (isExistF == idxEnd || isExistS == idxEnd)
 	{
 		Logger<<"���������еĿ�.\n";
		return;
 	}

	IndexType fnid = components_[tgFrame].hier_label_vtxBucket_index[tgLevel][fpId];
	IndexType snid = components_[tgFrame].hier_label_vtxBucket_index[tgLevel][spid];

	tgBestSet.insert(fnid);
	tgBestSet.insert(snid);

	if (tgBestSet.size() > 1)
	{
		Logger<<" Start Merge target Graph"<<tgFrame<<endl;

	    mergeTgPatches(tgLevel,tgBestSet);

		Logger<<" End Merge target Graph.\n";
	}


	printf("End to merge patch %d and %d.\n",fpId,spid);
}

ScalarType GraphMatch::distBetween2Setpatches(const set<IndexType>& sPatches, const set<IndexType>& tPatches,IndexType srLevel, IndexType tgLevel)
{
	assert(sPatches.size() > 0 && tPatches.size() > 0);

	map<IndexType,HVertex*> srVertex;

	set<IndexType>::iterator sIter = sPatches.begin();


	for (; sIter != sPatches.end(); sIter ++)
	{

		HLabel* temp = components_[srFrame].hier_label_bucket[srLevel][*sIter];

		srVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
	}

	map<IndexType,HVertex*> tgVertex;

	set<IndexType>::iterator tIter = tPatches.begin();

	for (; tIter != tPatches.end(); tIter++)
	{
		HLabel* temp = components_[tgFrame].hier_label_bucket[tgLevel][*tIter];
		tgVertex.insert(temp->vertex_bucket.begin(), temp->vertex_bucket.end() );
	}

	map<IndexType,HVertex*> toComVtx; // ��¼��Щ�����ǵĶ�Ӧ�����ڶ�Ӧ����

	map<IndexType,HVertex*> backComVtx; // �໥��¼

	bool isTo = true;
	computer_common(srVertex,tgVertex,toComVtx,isTo);

	isTo = false;
	computer_common(srVertex, tgVertex,backComVtx,isTo);


	ScalarType srV_size = srVertex.size();
	ScalarType tgV_size = tgVertex.size();

	IndexType toCom_size = toComVtx.size();
	IndexType backCom_size = backComVtx.size();

	ScalarType corValue =  0.5 * (toCom_size/srV_size + backCom_size/tgV_size );

	ScalarType scaleShapeValue =  abs( srV_size - tgV_size) / max(srV_size, tgV_size );


	Logger<<" Source combine patches.\n";
	for (auto  iter = sPatches.begin(); iter != sPatches.end(); ++ iter)
	{
		Logger<<(*iter)<<"  ";
	}
	Logger<<".\n";



	Logger<<" Target combine patches.\n";
	for (auto  iter = tPatches.begin(); iter != tPatches.end(); ++ iter)
	{
		Logger<<(*iter)<<"  ";
	}
	Logger<<".\n";

	//ScalarType totDis = corValue + scaleShapeValue + deformationValue;

	Logger<<" Two condition type distances.\n";
	Logger<<" corValue = "<<corValue<<endl;
	Logger<<" scaleShapeValue = "<<scaleShapeValue<<endl;
	//Logger<<" deformationvalue = "<<deformationValue<<endl;

	if ( scaleShapeValue > 0.6 || corValue < 0.5 )
	{
		Logger<<" The distance is INF.\n";

		return INF;

	}else
	{
// 		ScalarType deformationValue  = 0.0;
// 
// 		if (toComVtx.size() > 5 && backComVtx.size() > 5)
// 		{
// 		    deformationValue  = deformationDis(toComVtx,backComVtx);
// 
// 		}else
// 		{
// 			Logger<<"��̫С�޷�������ξ���, �ͱ�����.\n";
// 		}

		//Logger<<" deformationvalue = "<<deformationValue<<endl;

		return  /* deformationValue *//*+ scaleShapeValue  + */ (1 - corValue);

	}

}

void GraphMatch::printNiceMatchingPatches()
{
	IndexType bestNum = 30;

	while (bestNum -- > 0 && equeue_.size() > 0)
	{
		PatchMatch bestPatch = equeue_.top();
		equeue_.pop();

		set<IndexType> srPatches = bestPatch.srPatches;
		set<IndexType> tgPateches = bestPatch.tgPatches;

		Logger<<bestNum<<"�����ŵ����.\n";

		if (srPatches.size() > 0)
		{
			Logger<<" Merge  srFrame  patches.\n";
			for (auto  iter = srPatches.begin(); iter != srPatches.end(); ++ iter)
			{
				Logger<<(*iter)<<"  ";
			}
			Logger<<".\n";
		}

		if (tgPateches.size() > 0)
		{
			Logger<<" Merge  tgFrame  patches.\n";
			for (auto  iter = tgPateches.begin(); iter != tgPateches.end(); ++ iter)
			{
				Logger<<(*iter)<<"  ";
			}
			Logger<<".\n";
		}
	}
}

ScalarType GraphMatch::toDeformationErrAllVtx(map<IndexType,HVertex*>& toVtx, map<IndexType,HVertex*>& backVtx)
{
	Sample& oriFrame = sample_set_[srFrame];
	Sample& tarFrame = sample_set_[tgFrame];

	IndexType toPsSize = toVtx.size();
	IndexType backPsSize = backVtx.size();

	Matrix3X s_coord, t_coord;
	s_coord.setZero(3, toPsSize);
	t_coord.setZero(3, backPsSize);

	IndexType i = 0;

	for (auto toIter = toVtx.begin(); toIter != toVtx.end(); ++toIter,++i)
	{
		IndexType sPsId = toIter->first;
		s_coord.col(i) = oriFrame.vertices_matrix().col(sPsId);
	}

	IndexType j = 0;
	for (auto backIter = backVtx.begin(); backIter != backVtx.end(); ++backIter,++j)
	{

		IndexType sPsId = backIter->first;
		t_coord.col(j) = tarFrame.vertices_matrix().col(sPsId);
	}

	//loacl ICP 
	SICP::Parameters pa(false,2,10,1.2,1e5,20,20,1,1e-5);
	MatrixXXi localCor;
	localCor.setZero(1,toPsSize);
	Matrix3X oriCoor = s_coord;

	Matrix3X corCoor;

	SICP::point_w_point(s_coord,t_coord,localCor,pa);//front rotmate

	getCorrCoor(t_coord,corCoor,localCor);


 	Matrix33 rot_mat;
 	MatrixXX tran_vec;
 
 	point2point(oriCoor, corCoor, rot_mat, tran_vec);
 
 	//����ÿ���㵽��Ӧ�����
 	MatrixXX traSCoor = rot_mat * oriCoor;
 
 	traSCoor.colwise() += tran_vec.col(0);
 
 	MatrixXX errMat = corCoor - traSCoor;
 
 	ScalarType totDis = errMat.colwise().norm().sum();

	return totDis;
}

ScalarType GraphMatch::backDeformationErrAllVtx(map<IndexType,HVertex*>& toVtx, map<IndexType,HVertex*>& backVtx)
{
	Sample& oriFrame = sample_set_[srFrame];
	Sample& tarFrame = sample_set_[tgFrame];

	IndexType toPsSize = toVtx.size();
	IndexType backPsSize = backVtx.size();

	Matrix3X s_coord, t_coord;
	s_coord.setZero(3, toPsSize);
	t_coord.setZero(3, backPsSize);

	IndexType i = 0;

	for (auto toIter = toVtx.begin(); toIter != toVtx.end(); ++toIter,++i)
	{

		IndexType sPsId = toIter->first;
		s_coord.col(i) = oriFrame.vertices_matrix().col(sPsId);
	}

	IndexType j = 0;
	for (auto backIter = backVtx.begin(); backIter != backVtx.end(); ++backIter,++j)
	{

		IndexType sPsId = backIter->first;
		t_coord.col(j) = tarFrame.vertices_matrix().col(sPsId);
	}

	//loacl ICP 
	SICP::Parameters pa(false,2,10,1.2,1e5,20,20,1,1e-5);
	MatrixXXi localCor;
	localCor.setZero(1,backPsSize);

	Matrix3X oriCoor = t_coord;

	Matrix3X corCoor;

	SICP::point_w_point(t_coord,s_coord,localCor,pa);//front rotmate

	getCorrCoor(s_coord,corCoor,localCor);


 	Matrix33 rot_mat;
 	MatrixXX tran_vec;
 
 	point2point(oriCoor, corCoor, rot_mat, tran_vec);
 
 	//����ÿ���㵽��Ӧ�����
 	MatrixXX traSCoor = rot_mat * oriCoor;
 
 	traSCoor.colwise() += tran_vec.col(0);
 
 	MatrixXX errMat = corCoor - traSCoor;
 
 	ScalarType totDis = errMat.colwise().norm().sum();
 
 	return totDis;
}
void GraphMatch::getCorrCoor(Matrix3X& tgCoor, Matrix3X& corrCoor, MatrixXXi& vtxMap)
{
	IndexType vtxSize = vtxMap.cols();

	corrCoor.setZero(3, vtxSize);

	for (IndexType i = 0; i < vtxSize; ++i)
	{
		corrCoor.col(i) = tgCoor.col( vtxMap(0,i) );
	}

}

void GraphMatch::mergeTinyPatches(IndexType frameId, IndexType outlierSize)
{

	HFrame& cur_frame = components_[frameId];

	IndexType gLevel = cur_frame.hier_graph.size();//�ϲ���߲��е�С��

	--gLevel;

	vector<HLabel*> vtxBucket = cur_frame.hier_label_bucket[gLevel];

	IndexType minSize = 1e5;

	for (auto vIter = vtxBucket.begin(); vIter != vtxBucket.end(); ++vIter)
	{
		IndexType vSize= (*vIter)->vertex_bucket.size();
		if (vSize < minSize)
		{
			minSize = vSize;
		}
	}

	while (minSize < outlierSize)
	{
		map<IndexType,IndexType>& indexMap = cur_frame.hier_label_vtxBucket_index[gLevel];

		LabelsGraph* curGraph = cur_frame.hier_graph[gLevel];

		VertexIterator vBeginIter,vEndIter;

		AdjacencyIterator vAdjBeginIter,vAdjEndIter;

		tie(vBeginIter,vEndIter) = boost::vertices(*curGraph);

		set<IndexType> mergeS;

		IndexType i = 0;

		IndexType tempSize = 1e5;

		vector<IndexType> minAdj;

		for (auto vIter = vtxBucket.begin(); vIter != vtxBucket.end(); ++vIter,++i)
		{
			IndexType vSize= (*vIter)->vertex_bucket.size();
			if( vSize > outlierSize) continue;
			mergeS.clear();

			VertexIterator nodeIter = (vBeginIter + i);
			VertexDescriptor nodeDec = *nodeIter;
			auto adjIter = boost::adjacent_vertices(nodeDec,*curGraph);

			mergeS.insert(*nodeIter);

			minAdj.clear();
			for (; adjIter.first != adjIter.second; ++ adjIter.first)
			{
				IndexType adjSize = vtxBucket[*(adjIter.first)]->vertex_bucket.size();	
				if (adjSize > vSize) //�ϲ���С�Ŀ���
				{
					if (minAdj.empty())
					{
						minAdj.push_back(*(adjIter.first) );
					}else
					{
						auto adjB =minAdj.begin();
						IndexType bSize = vtxBucket[*adjB]->vertex_bucket.size();
						if (adjSize < bSize)
						{
							minAdj.insert(adjB,*(adjIter.first) );
						}
					}
				}

			}

			if (!minAdj.empty())
			{
				mergeS.insert( *(minAdj.begin() ) );
			}

			if ( mergeS.size() > 1)
			{
				mergePatchesOri(gLevel,frameId,mergeS);
				break;
			}

		} //end for vtx bucket;


		//�ж��Ƿ����merge����

		gLevel = cur_frame.hier_graph.size();//�ϲ���߲��е�С��

		--gLevel;

		vtxBucket = cur_frame.hier_label_bucket[gLevel];

		minSize = 1e5;

		for (auto vIter = vtxBucket.begin(); vIter != vtxBucket.end(); ++vIter)
		{
			IndexType vSize= (*vIter)->vertex_bucket.size();
			if (vSize < minSize)
			{
				minSize = vSize;
			}
		}

	}

}
