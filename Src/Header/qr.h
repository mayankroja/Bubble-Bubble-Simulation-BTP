/*
 * qr.h
 *      Author: sunder
 */

#ifndef QR_H_
#define QR_H_

#include <iostream>
#include <cmath>

/* Forward declarations */ 

template<typename T>
void allocate_mem_1d(T*&, int); 

template<typename T>
void allocate_mem_2d(T*&, int, int);

template<typename T>
void free_mem_1d(T*&);

template<typename T>
void free_mem_2d(T*&, int);

template<typename Real>
class QRdcmp; 

/* Allocate and deallocate memory for 1D and 2D arrays */

// Allocate memory for 1D array

template<typename T>
void allocate_mem_1d(T*& U, int m) {
	U = new T[m];
}


// Allocate memory for 2D array

template<typename T>
void allocate_mem_2d(T**& U, int m, int n) {
	U = new T*[m];

	for (int i = 0; i < m; ++i) {
		U[i] = new T[n];
	}
}

// Release memory for 1D array

template<typename T>
void free_mem_1d(T*& U) {
	delete[] U;
}

// Release memory for 2D array

template<typename T>
void free_mem_2d(T**& U, int m) {

	for (unsigned int i = 0; i < m; ++i) {
		delete[] U[i];
	}

	delete[] U;
}

/* QR decomposition class */ 

template <typename Real>
class QRdcmp {
	int m, n;
	Real** QR_;
	Real* Rdiag;
public:
	/* Constructors and destructor */

	QRdcmp();
	QRdcmp(Real**, int, int);
	QRdcmp(const QRdcmp&);
	~QRdcmp();

	/* Main methods of the class */

	bool is_full_rank() const;
	void get_Q(Real**) const;
	void get_R(Real**) const;
	void solve(const Real*, Real*) const;
	void Rsolve(const Real*, Real*) const;
};

//----------------------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------------------

template <typename Real>
QRdcmp<Real>::QRdcmp() :
	m(0),
	n(0),
	QR_(NULL),
	Rdiag(NULL)
{}

//----------------------------------------------------------------------------
// Constructor taking a C-style 2D array of size m x n as input
//----------------------------------------------------------------------------

template <typename Real>
QRdcmp<Real>::QRdcmp(Real** A, int mm, int nn) :
	m(mm),
	n(nn)
{

	// allocate memory to QR_ and Rdiag

	int i, j, k;

	QR_ = new Real*[m];

	for (int i = 0; i < m; ++i)
		QR_[i] = new Real[n];

	Rdiag = new Real[n];

	// Copy A into QR

	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			QR_[i][j] = A[i][j];

	i = 0; j = 0; k = 0;

    // Main loop.
    for (k = 0; k < n; k++) {
       // Compute 2-norm of k-th column without under/overflow.
       Real nrm = 0;
       for (i = k; i < m; i++) {
//          nrm = std::hypot(nrm,QR_[i][k]);
          nrm = sqrt(nrm*nrm + QR_[i][k]*QR_[i][k]);
       }

       if (nrm != 0.0) {
          // Form k-th Householder vector.
          if (QR_[k][k] < 0) {
             nrm = -nrm;
          }
          for (i = k; i < m; i++) {
             QR_[i][k] /= nrm;
          }
          QR_[k][k] += 1.0;

          // Apply transformation to remaining columns.
          for (j = k+1; j < n; j++) {
             Real s = 0.0;
             for (i = k; i < m; i++) {
                s += QR_[i][k]*QR_[i][j];
             }
             s = -s/QR_[k][k];
             for (i = k; i < m; i++) {
                QR_[i][j] += s*QR_[i][k];
             }
          }
       }
       Rdiag[k] = -nrm;
    }
}

//----------------------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------------------

template <typename Real>
QRdcmp<Real>::QRdcmp(const QRdcmp& QR2) :
	m(QR2.m),
	n(QR2.n)
{

	// allocate memory to QR_ and Rdiag

	int i, j, k;

	QR_ = new Real*[m];

	for (int i = 0; i < m; ++i)
		QR_[i] = new Real[n];

	Rdiag = new Real[n];

	// Copy A into QR

	for (i=0; i<m; i++)
		for (j=0; j<n; j++)
			QR_[i][j] = QR2.QR_[i][j];

	for (i=0; i<n; ++i)
		Rdiag[i] = QR2.Rdiag[i];
}

//----------------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------------

template <typename Real>
QRdcmp<Real>::~QRdcmp() {


	for (int i = 0; i < m; ++i)
		delete[] QR_[i];

	delete[] QR_;


	delete[] Rdiag;
}

//----------------------------------------------------------------------------
// Check if the matrix is full rank
//----------------------------------------------------------------------------

template <typename Real>
bool QRdcmp<Real>::is_full_rank() const {

	for (int j = 0; j < n; j++) {
       if (std::abs(Rdiag[j]) < 1.0e-14)
          return false;
    }
    return true;
}

//----------------------------------------------------------------------------
// Get the Q(m x n) matrix
//----------------------------------------------------------------------------

template <typename Real>
void QRdcmp<Real>::get_Q(Real** Q) const {

	int i=0, j=0, k=0;

    for (k = n-1; k >= 0; k--) {
       for (i = 0; i < m; i++) {
          Q[i][k] = 0.0;
       }
       Q[k][k] = 1.0;
       for (j = k; j < n; j++) {
          if (QR_[k][k] != 0) {
             Real s = 0.0;
             for (i = k; i < m; i++) {
                s += QR_[i][k]*Q[i][j];
             }
             s = -s/QR_[k][k];
             for (i = k; i < m; i++) {
                Q[i][j] += s*QR_[i][k];
             }
          }
       }
    }
}

//----------------------------------------------------------------------------
// Get the R(n x n) matrix
//----------------------------------------------------------------------------

template <typename Real>
void QRdcmp<Real>::get_R(Real** R) const {

    for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++) {
          if (i < j) {
             R[i][j] = QR_[i][j];
          } else if (i == j) {
             R[i][j] = Rdiag[i];
          } else {
             R[i][j] = 0.0;
          }
       }
    }
}

//----------------------------------------------------------------------------
// Solve system Ax = b; b -> m length array, x -> n length array
//----------------------------------------------------------------------------

template <typename Real>
void QRdcmp<Real>::solve(const Real* b, Real* x_) const {

	Real* x = new Real[m];

	// Copy b into x

	for (int i = 0; i < m; ++i)
		x[i] = b[i];

    // Compute Y = transpose(Q)*b
    for (int k = 0; k < n; k++)
	  {
          Real s = 0.0;
          for (int i = k; i < m; i++)
			{
             s += QR_[i][k]*x[i];
          }
          s = -s/QR_[k][k];
          for (int i = k; i < m; i++)
			{
             x[i] += s*QR_[i][k];
          }
    }

    // Solve R*X = Y;
    for (int k = n-1; k >= 0; k--)
	  {
       x[k] /= Rdiag[k];
       for (int i = 0; i < k; i++) {
             x[i] -= x[k]*QR_[i][k];
       }
    }


	// return n x nx portion of X

	for (int i=0; i<n; i++)
		x_[i] = x[i];


	  delete[] x;
}

//----------------------------------------------------------------------------
// Solve system Rx = b; b -> n length array, x -> n length array, R is the n X n matrix
//----------------------------------------------------------------------------

template <typename Real>
void QRdcmp<Real>::Rsolve(const Real* b, Real* x_) const {

	Real* x = new Real[n];

	// Copy b into x

	for (int i = 0; i < n; ++i)
		x[i] = b[i];


    // Solve R*X = Y;
    for (int k = n-1; k >= 0; k--)
	  {
       x[k] /= Rdiag[k];
       for (int i = 0; i < k; i++) {
             x[i] -= x[k]*QR_[i][k];
       }
    }


	// return n x nx portion of X

	for (int i=0; i<n; i++)
		x_[i] = x[i];


	  delete[] x;
}

#endif /* QR_H_ */
