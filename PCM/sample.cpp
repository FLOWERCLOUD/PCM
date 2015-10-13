#include "sample.h"
#include "windows.h"
#include <gl/gl.h>
#include <gl/glu.h>
#include <fstream>
#include "globals.h"
#include <assert.h>

Sample::~Sample()
{ 
	vertices_.clear();
	allocator_.free_all(); 
	delete	kd_tree_;
}

Vertex* Sample::add_vertex(const PointType& pos = NULL_POINT,
						const NormalType& n = NULL_NORMAL, 
						const ColorType& c = NULL_COLOR)
{
	Vertex*	new_space = allocator_.allocate<Vertex>();
	Vertex* new_vtx = new(new_space)Vertex;
	if ( !new_vtx )
	{
		return nullptr;
	}
	int a = sizeof(Vertex);
	int b = sizeof(SelectableItem);
	new_vtx->set_position(pos);
	new_vtx->set_normal(n);
	new_vtx->set_color(c);

	vertices_.push_back(new_vtx);
	box_.expand( pos );
	kd_tree_should_rebuild_ = true;

	return new_vtx;
}

void Sample::draw(ColorMode::ObjectColorMode, const Vec3& bias )
{
	if (!visible_)
	{
		return;
	}
	
	glPointSize(Paint_Param::g_point_size);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);


	if ( selected_ )
	{
		ColorType c = HIGHTLIGHTED_COLOR;
		glColor4f(c(0), c(1), c(2),c(3));
	}
	else
		glColor4f( color_(0), color_(1), color_(2), color_(3) );
	
	Matrix44 mat = matrix_to_scene_coord();
	for ( IndexType i = 0; i < vertices_.size(); i++ )
	{
		vertices_[i]->draw_without_color(mat, bias);
	}
	glEnd();

}

void Sample::draw(ColorMode::VertexColorMode, const Vec3& bias)
{
	if (!visible_)
	{
		return;
	}

	glPointSize(Paint_Param::g_point_size);
	glEnable(GL_POINT_SMOOTH);
	glBegin(GL_POINTS);


	Matrix44 mat = matrix_to_scene_coord();
	for ( IndexType i = 0; i < vertices_.size(); i++ )
	{
		vertices_[i]->draw(mat,bias);
	}
	glEnd();

}

void Sample::draw(ColorMode::LabelColorMode, const Vec3& bias)
{
	if (!visible_)
	{
		return;
	}
	//if ( selected_ )
	{
		//glDisable(GL_LIGHTING);
		glPointSize(Paint_Param::g_point_size);
		glEnable(GL_POINT_SMOOTH);
		glBegin(GL_POINTS);

		Matrix44 mat = matrix_to_scene_coord();
		for ( IndexType i = 0; i < vertices_.size(); i++ )
		{
			vertices_[i]->draw_with_label(mat,bias);
		}
		glEnd();
		//glEnable(GL_LIGHTING);
	}

}

void Sample::draw_with_name()
{
	Matrix44 mat = matrix_to_scene_coord();
	for( unsigned int idx = 0; idx < vertices_.size(); idx++ )
	{
		vertices_[idx]->draw_with_name(idx,mat);
	}
}

void Sample::build_kdtree()
{
	if( !kd_tree_should_rebuild_  || vertices_.size() == 0 )
	{
		return;
	}


	if (kd_tree_)
	{
		delete kd_tree_;
	}

	size_t n_vtx = vertices_.size();
	vtx_matrix_ = Matrix3X( 3, n_vtx );

	//reconstruct vtx_matrix
	for ( IndexType v_idx = 0; v_idx < n_vtx; v_idx++ )
	{
		Vertex*	pv = vertices_[v_idx];
		vtx_matrix_.col(v_idx) << pv->x() , pv->y() , pv->z();
	}

	kd_tree_ = new nanoflann::KDTreeAdaptor<Matrix3X, 3>(vtx_matrix_);

	kd_tree_should_rebuild_ = false;


}

IndexType Sample::closest_vtx( const PointType& query_point ) const
{
	ScalarType	qp[3] = { query_point(0), query_point(1), query_point(2) };

	if (kd_tree_should_rebuild_)
	{
		return -1;
	}

	return kd_tree_->closest( qp );
	

}

Matrix44 Sample::matrix_to_scene_coord( )
{
	Matrix44 mat;

	const ScalarType scale_factor = 1./box_.diag();
	const PointType	 box_center = box_.center();

	mat << scale_factor, 0, 0, -box_center(0)*scale_factor,
			0, scale_factor, 0, -box_center(1)*scale_factor,
			0, 0, scale_factor, -box_center(2)*scale_factor,
			0, 0, 0, 1;

	return mat;
						
}

bool Sample::neighbours(const IndexType query_point_idx, const IndexType num_closet, IndexType* out_indices)
{
	if (kd_tree_should_rebuild_)
	{
		return false;
	}

	ScalarType* out_distances = new ScalarType[num_closet];
	ScalarType	qp[3] = {vertices_[query_point_idx]->x(), 
						vertices_[query_point_idx]->y(), 
						vertices_[query_point_idx]->z() };
	kd_tree_->query( qp, num_closet, out_indices, out_distances);

	delete out_distances;

	return true;
}

bool Sample::neighbours(const IndexType query_point_idx, const IndexType num_closet,
						IndexType* out_indices,ScalarType* out_distances)
{
	if (kd_tree_should_rebuild_)
	{
		return false;
	}

	ScalarType	qp[3] = {vertices_[query_point_idx]->x(), 
		vertices_[query_point_idx]->y(), 
		vertices_[query_point_idx]->z() };
	kd_tree_->query( qp, num_closet, out_indices, out_distances);

	return true;
}

void Sample::update()
{
	assert( vtx_matrix_.cols() == vertices_.size() );
	IndexType v_idx = 0;
	box_ = Box();
	for ( vtx_iterator v_iter = begin();
			v_iter != end(); v_iter++,v_idx++ )
	{
		PointType p( vtx_matrix_(0, v_idx),
					vtx_matrix_(1, v_idx),
					vtx_matrix_(2, v_idx));
		(*v_iter)->set_position( p );
		box_.expand( p );
	}
	kd_tree_should_rebuild_ = true;
	build_kdtree();
}

void Sample::delete_vertex_group(const std::vector<IndexType>& idx_grp )
{
	IndexType  i=0, j=0;
	IndexType size = idx_grp.size();
	if (size==0)
	{
		return;
	}
	for ( std::vector<Vertex*>::iterator iter = vertices_.begin();
		iter != vertices_.end(); i++)
	{
		if ( i == idx_grp[j] )
		{
			//This is the node to delete
			iter = vertices_.erase( iter );
			j++;
			if ( j>=size )
			{
				break;
			}
		}
		else
		{
			iter++;
		}
	}

	//kdtree dirty
	kd_tree_should_rebuild_ = true;
	build_kdtree();

}

