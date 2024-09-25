#include "Stencil.H"
#include <AMReX_Print.H>

namespace mycode
{

Stencil::Stencil()
{
}

Stencil::Stencil(unsigned int size)
{
    col_entry_.resize(size);
    source_ = 0.0;
}

void Stencil::define(unsigned int size)
{
    col_entry_.resize(size);
    source_ = 0.0;
}

std::vector<int> Stencil::getColIdVec()
{
    std::vector<int> col_id;
    col_id.reserve(col_entry_.size());
    for( const auto& col : col_entry_ )
    {
        col_id.emplace_back(col.first);
    }
    return col_id;
}

const std::vector<int> Stencil::getColIdVec() const
{
    std::vector<int> col_id;
    col_id.reserve(col_entry_.size());
    for( const auto& col : col_entry_ )
    {
        col_id.emplace_back(col.first);
    }
    return col_id;
}

std::vector<double> Stencil::getCoeffVec()
{
    std::vector<double> coeff;
    coeff.reserve(col_entry_.size());
    for( const auto& col : col_entry_ )
    {
        coeff.emplace_back(col.second);
    }
    return coeff;
}

const std::vector<double> Stencil::getCoeffVec() const
{
    std::vector<double> coeff;
    coeff.reserve(col_entry_.size());
    for( const auto& col : col_entry_ )
    {
        coeff.emplace_back(col.second);
    }
    return coeff;
}

void Stencil::setRowId(int i)
{
    row_id_ = i;
}

void Stencil::setSource(double val)
{
    source_ = val;
}

void Stencil::setColEntry(unsigned int i, int col_id, double val)
{
    assert(!col_entry_.empty());
    col_entry_[i] = std::pair<int, double>(col_id, val);
}

void Stencil::addToColEntry(unsigned int i, int col_id, double val)
{
    if( col_entry_.size() <= i )
    {
        push_back_colEntry(col_id, val);
    }
    else
    {
        col_entry_[i].second += val; 
    }
}

void Stencil::addToColEntry(int col_id, double val)
{
	auto it = std::find_if(col_entry_.begin(), col_entry_.end(), 
	                       [&]( const std::pair< int, double >& col )
	                       { return col.first == col_id; } 
	                      );
	
	if(it != col_entry_.end())
	{
		it->second += val;
	}
	else
	{
		push_back_colEntry(col_id, val);
	}
}

void Stencil::addToSource(double val)
{
    source_ += val;
}

void Stencil::push_back_colEntry(int col_id, double val)
{
    col_entry_.push_back( std::pair<int, double>(col_id, val) );
}

void Stencil::Print()
{
    amrex::Print() << "row id : " << row_id_ << "\n";
    amrex::Print() << "\t source : " << source_ << "\n";

    for (auto &&col : col_entry_)
    {
        amrex::Print() << "\tcol : " << col.first << "\tcoeff : " << col.second << "\n";
    }
}

void Stencil::Print() const
{
    amrex::Print() << "row id : " << row_id_ << "\n";
    amrex::Print() << "\t source : " << source_ << "\n";
    for (auto &&col : col_entry_)
    {
        amrex::Print() << "\tcol : " << col.first << "\tcoeff : " << col.second << "\n";
    }
}

Stencil::~Stencil()
{
    if(!col_entry_.empty())
    {
        col_entry_.clear();
        col_entry_.shrink_to_fit();
    }
}

Stencil operator + ( const Stencil& A, const Stencil& B )
{
    assert(A.getRowId() == B.getRowId());
    Stencil C = A;
    C += B;
    return C;
}

Stencil operator - ( const Stencil& A, const Stencil& B )
{
    assert(A.getRowId() == B.getRowId());
    Stencil C = A;
    C -= B;
    return C;
}

Stencil operator == ( const Stencil& A, const Stencil& B )
{
	return A-B;
}

} /*End namespace mycode */

