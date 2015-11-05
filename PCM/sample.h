#ifndef _EXAMPLE_H
#define _EXAMPLE_H
#include "selectable_item.h"
#include "vertex.h"
#include "pool_allocator.h"
#include "box.h"
#include "..\ICP\ICP.h"
#include "Eigen\Dense"
#include "basic_types.h"
#include <QMutex>

class Sample:public SelectableItem
{
public:
	Sample():vertices_(),allocator_(),kd_tree_(nullptr),
				kd_tree_should_rebuild_(true),
				mutex_(QMutex::NonRecursive){}
	~Sample();

	inline Vertex& operator[]( IndexType i ) const{ return *vertices_[i]; }
	

	Vertex* add_vertex( const PointType& pos,const NormalType& n,
		const ColorType& c);

	void delete_vertex_group( const std::vector<IndexType>& idx_grp );

	void draw(ColorMode::ObjectColorMode, const Vec3& bias = Vec3(0.,0.,0.));
	void draw(ColorMode::VertexColorMode,const Vec3& bias = Vec3(0.,0.,0.));
	void draw(ColorMode::LabelColorMode,const Vec3& bias = Vec3(0.,0.,0.));
	void draw_with_name();

	size_t num_vertices() const { return vertices_.size(); }

	typedef std::vector<Vertex*>::iterator	vtx_iterator;

	inline vtx_iterator begin(){ return vertices_.begin(); }
	inline vtx_iterator end(){ return vertices_.end(); }

	//Every time vertex change, the kdtree should rebuild
	void	build_kdtree();

	IndexType closest_vtx( const PointType& query_point ) const;
	bool		neighbours(const IndexType query_point_idx, const IndexType num_closet, IndexType* out_indices);
	/* 
		Get matrix for transforming world-sample space to 
		view-sample space , making sure all samples can be saw
		no matter what original coordinates it is
	*/
	inline Matrix44 matrix_to_scene_coord(  );

	/* Green channel to get all vertex position information */
	inline  Matrix3X&	vertices_matrix()
	{	
		if (kd_tree_should_rebuild_)
		{
			build_kdtree();
		}
		return vtx_matrix_; 
	}
	/*Update vertex position according vertex matrix*/
	void	update();

	inline void lock(){ mutex_.lock(); }
	inline void unlock(){ mutex_.unlock(); }

	bool Sample::neighbours(const IndexType query_point_idx, const IndexType num_closet,
		IndexType* out_indices,ScalarType* out_distances);


	const PointType box_center() const{ return box_.center(); }
	const ScalarType	box_diag() { return box_.diag(); }
	const PointType box_near_corner(){ return box_.low_corner(); }
	const PointType box_far_corner(){ return box_.high_corner(); }

	Box getBox(){return box_;}

public:
	inline Matrix3X& nor_matrix()
	{
		if (kd_tree_should_rebuild_)
		{
			build_kdtree();
		}
		return nor_matrix_;
	}
// public:
// 	QMutex										mutex_;

private:
	std::vector<Vertex*>	vertices_;
	PoolAllocator			allocator_;
	Box						box_;

	Matrix3X									vtx_matrix_; // 3 * NumOfVtx
	//
	Matrix3X                                    nor_matrix_;
	//
	//Attention: kdtree is just a adapter, it means it use the reference of its data source
	nanoflann::KDTreeAdaptor<Matrix3X, 3>*		kd_tree_;
	bool										kd_tree_should_rebuild_;
	QMutex										mutex_;
};

#endif