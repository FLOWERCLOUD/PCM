#ifndef _GRAPH_MATCHING_H
#define _GRAPH_MATCHING_H

#include "pool_allocator.h"
#include "BoostGraph.h"
#include "sample_set.h"

#define  COUT_DEBUG 0

struct HVertex;
struct HLabel;
struct HFrame;


class GraphMatch
{

public:

	GraphMatch(SampleSet& _smpSet, map<IndexType,HFrame>& _componet,IndexType _srFrame,IndexType _tgFrame): sample_set_(_smpSet), components_(_componet)
	{
		srFrame = _srFrame;
		tgFrame = _tgFrame;
	}

public:
	void GraphMatch::computer_common(map<IndexType,HVertex*>& srPatchesVtx, map<IndexType,HVertex*>& tgPatchesVtx, 
		                             map<IndexType,HVertex*>& comVtx,bool isTo);

	ScalarType distance2PatchesLevel();

	void point2point(Matrix3X & srCloud,Matrix3X & tgCloud,Matrix33 & rotMat,MatrixXX & transVec);

	ScalarType deformationDis(map<IndexType,HVertex*>& toComVtx, map<IndexType,HVertex*>& backComVtx);

	ScalarType distBetween2Frame4Items(const set<IndexType>& sPatches, const set<IndexType>& tPatches,IndexType srLevel, IndexType tgLevel);

	void process_CFcmp_node( IndexType CFrameId ,IndexType NFrameId, set<IndexType>& CFcmp_node /* 当前帧的组合节点*/ ,set<IndexType>& NFcmp_node ,
		                     LabelsGraph* labelGraphOfNFrame,IndexType graLevel);

	void matchNextFrameOnePt(set<IndexType>& CFcmp_node,IndexType srGraLevel,IndexType tgGraLevel);

	void calculateSimilar2Frame();

	ScalarType distance2PatchesLevel(const set<IndexType>& sPatches,const set<IndexType>& tPatches,IndexType srLevel, IndexType tgLevel);
	                                 
	void matchNextFrameSet(set<IndexType>& curNodeSet,LabelsGraph& tgGraph,IndexType srGraLevel,IndexType tgGraLevel);

	void mergePatches(IndexType srLevel,IndexType tgLevel);

	void mergeSrPatches(IndexType srLevel,set<IndexType>& srBestCombine);

	void mergeTgPatches(IndexType tgLevel,set<IndexType>& tgBestCombine);

	void setMatchset(IndexType srLevel, IndexType tgLevel);

	IndexType mergePatchesOri(IndexType depth, IndexType frameId, set<IndexType>& patches);
	LabelsGraph* mergeIter(IndexType gLevel,IndexType frameId, LabelsGraph& _orIG ,vector<HLabel*>& old_label_bucket , 
		                   set<IndexType>& patches ,IndexType& new_label,map<IndexType,IndexType>& LabelIdIndex);
	vector<HLabel*>* DeepCopyHLableVec( const vector<HLabel*>* _p_oriHLabelvec );
	HLabel* DeepCopyHLabel( const HLabel* _p_oriHlabel );
	HVertex* DeepCopyVertex( const HVertex* _p_orivtx , HLabel* new_label_parent );
	void mergeTwoLabelAndPushBack(HLabel& _cLabel , HLabel& _bestNeiLabel_ );
	void processAndCreateNewGraphAndPushBack(LabelsGraph& _preGraph , VertexDescriptor& _targetVd ,vector<HLabel*>& _labelvec_ , map<IndexType ,IndexType>& labelofvtx_, LabelsGraph& newGraph_);
	void copyLabelBucket(vector<HLabel*>& leftLabels, const vector<HLabel*>& oriLabels);

	IndexType mergeAllNodes(IndexType frameId, set<IndexType>& patches);
	IndexType mergeAllSetOneTime(LabelsGraph& _orIG ,vector<HLabel*>& old_label_bucket , set<IndexType>& patches ,IndexType& new_label );

	ScalarType toDeformationErr(map<IndexType,HVertex*>& toComVtx);
	ScalarType backDeformationErr(map<IndexType,HVertex*>& backeComVtx);

	void buildPatchCorrespondenceByLabel();

	void printMatchingInformation();

	//20150831
	void findBestPatches(IndexType srLevel, IndexType tgLevel,set<IndexType>& srBestSet, set<IndexType>& tgBestSet); //same as setMatchset function

	void mergePatchesAfterCoSeg(IndexType srLevel, IndexType tgLevel,set<IndexType>& srBestSet, set<IndexType>& tgBestSet);

	ScalarType distBetween2Setpatches(const set<IndexType>& sPatches, const set<IndexType>& tPatches,IndexType srLevel, IndexType tgLevel);

	void printNiceMatchingPatches();

	ScalarType toDeformationErrAllVtx(map<IndexType,HVertex*>& toVtx, map<IndexType,HVertex*>& backVtx);

	ScalarType backDeformationErrAllVtx(map<IndexType,HVertex*>& toVtx, map<IndexType,HVertex*>& backVtx);

	void getCorrCoor(Matrix3X& tgCoor, Matrix3X& corrCoor, MatrixXXi& vtxMap);

	void mergeTinyPatches(IndexType frameId, IndexType outlierSize);

public:
	IndexType srFrame;
	IndexType tgFrame;

private:
	SampleSet& sample_set_;
	map<IndexType ,HFrame>& components_;

	set< set<IndexType> > CFrameSetSet_;
	set< set<IndexType> > NFrameSetSet_;

	std::map< IndexType , set< set<IndexType> > >  stotal_map_;
	std::map< IndexType , set< set<IndexType> > >  ttotal_map_;
	std::map< IndexType , set< IndexType > >  sadj_map_;
	std::map< IndexType , set< IndexType>  >  tadj_map_;

	typedef priority_queue< PatchMatch, vector<PatchMatch>, patchmatch_error_greater_than> ERQueue;
	ERQueue equeue_;
private:
	PoolAllocator allocator_;
};

#endif // !_GRAPH_MATCHING_H
