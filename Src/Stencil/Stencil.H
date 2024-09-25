#ifndef _Stencil_H_
#define _Stencil_H_


/************************************************************************

Author: Mahesh K Sawardekar
Date: Mon Dec 30 18:26:32 2019


************************************************************************/

#include <vector>
#include <AMReX_IntVect.H>
#include <AMReX_Array4.H>
#include <AMReX_Print.H>

namespace mycode
{

/**
 * @brief Class to hold the stencil data in FV/FD discretization
 *        needed for setting up the linear system 
 */
class Stencil
{
// public members
public:
	enum cell_level
	{
		current,
		coarse,
		fine
	};

	// default constructor:
	/**
	 * @brief Construct a new empty Stencil object
	 * 
	 */
	Stencil();

// constructor:
	/**
	 * @brief Construct a new Stencil object of given size
	 * 
	 * @param size : stencil size
	 */
	Stencil(unsigned int size);

// destructor:
	~Stencil();

/**************************** functions *********************************/
	
	/**
	 * @brief allocate memory to col_entry of size 'size'
	 * 
	 * @param size : stencil size
	 */
	void define(unsigned int size);

	/**
	 * @brief return no of values in stencil
	 * 
	 * @return unsigned int 
	 */
	unsigned int size() { return col_entry_.size(); }
	const unsigned int size() const { return col_entry_.size(); }


	/**
	 * @brief Get the row number corresponding to global System A
	 * 
	 * @return int : row index
	 */
	int getRowId() { return row_id_; }
	const int getRowId() const { return row_id_; }

	/**
	 * @brief Get the Source
	 * 
	 * @return double 
	 */
	double getSource() { return source_; }
	const double getSource() const { return source_; }


	/**
	 * @brief Get the ith Col Entry
	 * 
	 * @param i 
	 * @return std::pair< int, double >& : col index and corresponding value
	 */
	std::pair< int, double >& getColEntry( unsigned int i ) { return col_entry_[i]; }
	const std::pair< int, double >& getColEntry( unsigned int i ) const { return col_entry_[i]; }


	/**
	 * @brief reference to all col entries in stencil
	 * 
	 */
	std::vector< std::pair< int, double > >& getColEntry() { return col_entry_; }
	const std::vector< std::pair< int, double > >& getColEntry() const { return col_entry_; }
	
	std::vector<int> getColIdVec();
	const std::vector<int> getColIdVec() const;
	
	std::vector<double> getCoeffVec();
	const std::vector<double> getCoeffVec() const;

	/**
	 * @brief Set the row index
	 * 
	 * @param i 
	 */
	void setRowId(int i);

	void setSource(double val);

	/**
	 * @brief Set the ith col_id and value of stencil
	 *        The stencil should be preallocated
	 * 
	 * @param i      : ith entry of stencil
	 * @param col_id : col index
	 * @param val    : value corresponding to col_id
	 */
	void setColEntry(unsigned int i, int col_id, double val);

	/**
	 * @brief add val to ith colEntry value
	 * 
	 * @param i 
	 * @param col_id 
	 * @param val 
	 */
	void addToColEntry(unsigned int i, int col_id, double val);
	void addToColEntry(int col_id, double val);
	void addToSource(double val);

	/**
	 * @brief Set the col_id and value of stencil 
	 *        Adds to the end of already present stencil val
	 * 
	 * @param col_id : col index
	 * @param val    : value corresponding to col_id
	 */
	void push_back_colEntry(int col_id, double val);

	void Print();
	void Print() const;

/**************************** operators *********************************/

	Stencil& operator * (const double val)
	{
		for (unsigned int i = 0; i < size(); i++)
		{
			col_entry_[i].second *= val;
			source_ *= val;
		}
		return *this;
	}

	Stencil& operator / (const double val)
	{
		for (unsigned int i = 0; i < size(); i++)
		{
			col_entry_[i].second /= val;
			source_ /= val;
		}
		return *this;
	}

	Stencil& operator += (const Stencil& A)
	{
		assert(A.getRowId() == row_id_);
		for (size_t i = 0; i < A.size(); i++)
		{
			addToColEntry(A.getColEntry(i).first, A.getColEntry(i).second);		
		}
		return *this;
	}

	Stencil& operator -= (const Stencil& A)
	{
		assert(A.getRowId() == row_id_);
		for (size_t i = 0; i < A.size(); i++)
		{
			addToColEntry(A.getColEntry(i).first, -A.getColEntry(i).second);
		}
		return *this;
	}

	Stencil& operator = (const Stencil& A)
	{
		if(this == &A)
		{
			return *this;
		}
		else
		{
			row_id_ = A.row_id_;
			source_ = A.source_;
			col_entry_ = A.col_entry_;
			return *this;
		}
	}

/************************************************************************/

// private members
private:

/****************************   data    *********************************/
	
	/**
	 * @brief row id in global system
	 * 
	 */
	int row_id_;

	/**
	 * @brief source value i.e. value that goes in rhs or b in system Ax = b
	 * 
	 */
	double source_;

	/**
	 * @brief col id and coeff corresponding to col id in global system
	 * 
	 */
	std::vector< std::pair< int, double > > col_entry_;

	/************************************************************************/

}; /*End class Stencil */

Stencil operator + ( const Stencil& A, const Stencil& B );
Stencil operator - ( const Stencil& A, const Stencil& B );
Stencil operator == ( const Stencil& A, const Stencil& B );

} /*End namespace mycode */

#endif
