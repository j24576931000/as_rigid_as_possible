#ifndef _ARAP_INTERACTOR_2_H__
#define _ARAP_INTERACTOR_2_H__

#include <map>
#include <set>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/numeric/bindings/umfpack/umfpack_inc.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
//using namespace cv;
using namespace std;

// Direct implementation of Igarashi's As-Rigid-As-Possible Shape Manipulation
// paper, without any tweaking or optimization
// 

// LU decomposition
template <class T, int N>
static inline bool ludcmp(T a[N][N], int indx[N], T *d = NULL)
{
	int i, j, k;
	T vv[N];

	if (d)
		*d = 1;
	for (i = 0; i < N; i++) {
		T big = 0.0;
		for (j = 0; j < N; j++) {
			T tmp = fabs(a[i][j]);
			if (tmp > big)
				big = tmp;
		}
		if (big == 0.0)
			return false;
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < N; j++) {
		for (i = 0; i < j; i++) {
			T sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		T big = 0.0;
		int imax = j;
		for (i = j; i < N; i++) {
			T sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			T tmp = vv[i] * fabs(sum);
			if (tmp > big) {
				big = tmp;
				imax = i;
			}
		}
		if (imax != j) {
			for (k = 0; k < N; k++)
				std::swap(a[imax][k], a[j][k]);
			if (d)
				*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0)
			return false;
		if (j != N-1) {
			T tmp = 1.0/(a[j][j]);
			for (i = j+1; i < N; i++)
				a[i][j] *= tmp;
		}
	}
	return true;
}


// Backsubstitution after ludcmp
template <class T, int N>
static inline void lubksb(T a[N][N], int indx[N], T b[N])
{
	int ii = -1, i, j;
	for (i = 0; i < N; i++) {
		int ip = indx[i];
		T sum = b[ip];
		b[ip] = b[i];
		if (ii != -1)
			for (j = ii; j < i; j++)
				sum -= a[i][j] * b[j];
		else if (sum)
			ii = i;
		b[i] = sum;
	}
	for (i = N-1; i >= 0; i--) {
		T sum = b[i];
		for (j = i+1; j < N; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}

// Matrix inverse
template <class T, int N>
static inline void inverse(T a[N][N]) {
	T y[N][N], d, col[N];
	int i, j, indx[N];

	ludcmp<double, N>(a, indx, &d);
	for (j = 0; j < N; j++) {
		for (i = 0;i < N; i++)
			col[i] = 0.0;
		col[j] = 1.0;
		lubksb<double, N>(a, indx, col);
		for(i = 0; i < N; i++)
			y[i][j] = col[i];
	}

	// copy back
	for ( int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			a[i][j] = y[i][j];
}

// Dot product of ith row of m1 and jth col of m2
template <class T, int M, int K, int N>
static T row_X_col(T m1[M][K], T m2[K][N], int i, int j)
{
	T r = 0;
	for (int k = 0; k < K; k++)
		r += m1[i][k] * m2[k][j];
	return r;
}

// Dot product of ith row of m1 and jth col of m2
template <class T, int M, int K, int N>
static void mul(T m1[M][K], T m2[K][N], T m[M][N])
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			m[i][j] = row_X_col<T, M, K, N>(m1, m2, i, j);
}

// Dot product of ith row of m1 and jth col of m2
template <class T, int M, int N>
static void mul(T A[M][N], T v[N], T b[M])
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			b[i] += A[i][j]*v[j];
}

// Transpose
template <class T, int M, int N>
static void transpose(T m[M][N], T mt[N][M])
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			mt[j][i] = m[i][j];
}


// Some utility function
template <class T>
T* new_vec(int n) {
	T* a = new T[n];
	memset(a, 0, sizeof(T)*n);
	return a;
}

template <class T>
T** new_matrix(int m, int n) {
	T**a = new T*[m];
	for (int i = 0; i < m; i++){
		a[i] = new T[n];
		memset(a[i], 0, sizeof(T)*n);
	}

	return a;
}

template <class T>
void delete_matrix(T** a, int m) {
	for (int i = 0; i < m; i++) {
		delete[] a[i];
		a[i] = NULL;
	}
	delete[] a;
}

template <class T>
T* copy(T* a, int n) {
	T *a1 = new T[n];
	memcpy(a1, a, sizeof(T)*n);
	return a1;
}

template <class T, int M, int N>
void copy(T a[M][N], T d[M][N]) {
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			d[i][j] = a[i][j];
}

template <class T>
T** copy(T** a, int m, int n) {
	T **a1 = new T*[m];
	for (int i = 0; i < m; i++) {
		a1[i] = new T[n];
		memcpy(a1[i], a[i], sizeof(T)*n);
	}
	return a1;
}

template <class T>
void print(FILE* f, const double v[], int n)
{
	for (int i = 0; i < n; i++)
		fprintf(f, "%g\n", v[i]);
}

template <class T, int M, int N>
void print(FILE* f, const double v[M][N])
{
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++)
			printf("%6.4f\t", v[i][j]);
		printf("\n");
	}
}

// Geometric primitive definitions
template <int D, class T = double>
class Vec {
protected:
	T v[D];

public:
	// Constructor for no arguments.  Everything initialized to 0.
	Vec() { for (int i = 0; i < D; i++) v[i] = T(0); }

	// Constructor for 1-4 arguments.  Ugly, but the compiler *should*
	// optimize away the ifs and unroll the for...  This one's explicit.
	explicit Vec(T x, T y = T(0), T z = T(0), T w = T(0))
	{ v[0] = x;  if (D > 1) v[1] = y;
	if (D > 2) v[2] = z;  if (D > 3) v[3] = w;
	for (int i = 4; i < D; i++) v[i] = T(0); }

	// Array reference and conversion to pointer - no bounds checking
	const T &operator [] (int i) const
	{ return v[i]; }
	T &operator [] (int i)
	{ return v[i]; }
	operator const T * () const
	{ return v; }
	operator const T * ()
	{ return v; }
	operator T * ()
	{ return v; }

	// Member operators
	Vec<D,T> &operator += (const Vec<D,T> &x)
	{ for (int i = 0; i < D; i++) v[i] += x[i];  return *this; }
	Vec<D,T> &operator -= (const Vec<D,T> &x)
	{ for (int i = 0; i < D; i++) v[i] -= x[i];  return *this; }
	Vec<D,T> &operator *= (const Vec<D,T> &x)
	{ for (int i = 0; i < D; i++) v[i] *= x[i];  return *this; }
	Vec<D,T> &operator *= (const T &x)
	{ for (int i = 0; i < D; i++) v[i] *= x;     return *this; }
	Vec<D,T> &operator /= (const Vec<D,T> &x)
	{ for (int i = 0; i < D; i++) v[i] /= x[i];  return *this; }
	Vec<D,T> &operator /= (const T &x)
	{ for (int i = 0; i < D; i++) v[i] /= x;     return *this; }

	// Outside of class: + - * / dist
};

// Nonmember operators that take two Vecs
template <int D, class T>
static inline const Vec<D,T> operator + (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return Vec<D,T>(v1) += v2;
}

template <int D, class T>
static inline const Vec<D,T> operator - (const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return Vec<D,T>(v1) -= v2;
}

// Vec/scalar operators
template <int D, class T>
static inline const Vec<D,T> operator * (const T &x, const Vec<D,T> &v)
{
	Vec<D,T> result(v);
	for (int i = 0; i < D; i++)
		result[i] = x * result[i];
	return result;
}

template <int D, class T>
static inline const Vec<D,T> operator * (const Vec<D,T> &v, const T &x)
{
	return Vec<D,T>(v) *= x;
}

template <int D, class T>
static inline const Vec<D,T> operator / (const T &x, const Vec<D,T> &v)
{
	Vec<D,T> result(v);
	for (int i = 0; i < D; i++)
		result[i] = x / result[i];
	return result;
}

template <int D, class T>
static inline const Vec<D,T> operator / (const Vec<D,T> &v, const T &x)
{
	return Vec<D,T>(v) /= x;
}

// Other functions

// Unitility functions first
template <class T>
static inline T sqr(const T &x)
{
	return x*x;
}

template <int D, class T>
static inline const T dist2(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	T d2 = sqr(v2[0]-v1[0]);
	for (int i = 1; i < D; i++)
		d2 += sqr(v2[i]-v1[i]);
	return d2;
}

// Euclidean distance between two points
template <int D, class T>
static inline const T dist(const Vec<D,T> &v1, const Vec<D,T> &v2)
{
	return sqrt(dist2(v1,v2));
}


typedef Vec<2> Point2D;
typedef Vec<2> Vec2D;
typedef Vec<3, int> Tri;


#define  GLUT_LEFT_BUTTON                   0x0000
#define  GLUT_MIDDLE_BUTTON                 0x0001
#define  GLUT_RIGHT_BUTTON                  0x0002
#define  GLUT_DOWN                          0x0000
#define  GLUT_UP                            0x0001
#define  GLUT_LEFT                          0x0000
#define  GLUT_ENTERED                       0x0001



typedef std::map<std::pair<int, int>, double > matrix_map;
typedef matrix_map::iterator mat_iter;
typedef matrix_map::const_iterator mat_citer;

typedef std::map<int, std::set<int> > iset_map;
typedef iset_map::iterator iset_iter;
typedef iset_map::const_iterator iset_citer;
typedef std::set<int>::iterator elt_iter;
typedef std::set<int>::const_iterator elt_citer;


namespace ublas = boost::numeric::ublas;
namespace umf = boost::numeric::bindings::umfpack;
class SparseMatrix
{
public:
	// m:	number of rows
	// n:	number of columns
    SparseMatrix(int m = 0, int n = 0);
    SparseMatrix(const SparseMatrix& B);
    ~SparseMatrix();

    void resetDim(int m, int n);


public:
    // Memeber reference
    double operator() (int i, int j) const;
    double& operator() (int i, int j);

    // Formated output, for debugging
    //void print(FILE* f) const;
    void print();
public:
    // algebra
    SparseMatrix transpose() const;
    SparseMatrix operator + (const SparseMatrix& B) const;

    // Multiplication, using own implementation
    vector<double> operator * (const vector<double>& x) const;

    // solve Ax = b for x
    vector<double> solve(const vector<double>& b);

    int m_rows;
    int n_cols;

    iset_map m_row_indices;
    iset_map m_col_indices;

    void PreSaveMatrix();
    double get( int i, int j ) const;

    ublas::compressed_matrix<double, 
	    ublas::column_major, 
	    0,
	    ublas::unbounded_array<int>, 
	    ublas::unbounded_array<double> > PreSave;
    matrix_map SMatrix_data;

};
// Another place holder that expects your implementation
class ShapeView
{
public:
	virtual Point2D window2world(const Point2D& win_pt)=0;
	virtual double window2world(double win_len)=0;
	virtual Point2D world2window(const Point2D& pos)=0;
	
	// Post an update-draw function
	virtual void refresh() =0;

public:
	enum ButtonState {
		LBUTTON_DOWN = (1<<GLUT_LEFT_BUTTON),
		MBUTTON_DOWN = (1<<GLUT_MIDDLE_BUTTON),
		RBUTTON_DOWN = (1<<GLUT_RIGHT_BUTTON)
	};
	ButtonState get_button_state() ;
	void button_state(int state);
private:
	int b_state;

};


// The trimesh in the 2D space - the representation that the ARAP used

class TriMesh2D
{
public:
	vector<Point2D> vertices;
	vector<vector<Point2D>> record_v;
	vector<Tri> tris;
	vector<int> normals;//1¥¿  -1­t

public:
	void compute_normal();
	int normal_detection();
	void draw(bool linemode);
};

class ArapInteractor
{
	// For user interaction (connection to your own code)
	ShapeView* shapeview;
	vector<Point2D> savelastFlagsPosition;
	vector<Point2D> boundary;
	TriMesh2D themesh;
	
	vector<Point2D> baseVertices;

	// Indicate the state of each vertices
	//	0:	clear
	//	1:	static
	vector<int> flags;

	// Indicates if it currently under dragging
	vector<int> ctrlPoints;
	bool beingDragged;

	// Algorithm
	void step1();
	void step2();

	// For step 1:
	// Helper functions
	void getGeachTri(const Tri& t, double G[][6]);
	// Precomputation right after the construction of the mesh
	void preCompG();
	// Precomputation before each dragging
	void preStep1();

	// For step 2:
	// Helper functions
	void getFeachTri(const Tri& t, double K[][4], double F[][4]);
	void getHeachTri(const Tri& t, double H[][6]);
	// Precomputation right after the construction of the mesh
	void preCompF();
	void preCompH();
	// Precomputation before each dragging
	void preStep2();

	// Do the deformation
	void deform();
	

public:
	ArapInteractor(ShapeView* view, const TriMesh2D& mesh);
	virtual ~ArapInteractor(void);

public:
	void recover_mesh();
	TriMesh2D tempmesh;
	vector<TriMesh2D>record_mesh;
	const char* get_string() { return "ARAP"; }
	void return_mesh(int);
	void save_mesh();
	void OnDraw(int vp = -1);
	void OnDisplay();
	void OnDrawInfo(int vp = -1);
	void OnMotion(int x, int y, int flag,bool mouse_down=0, int vp = -1);
	void OnMouse(int button, int button_state, int x, int y, int vp = -1);
	void OnKeyboard(unsigned char key, int x, int y);
	void OnSpecialKey(int key, int x, int y);

	//get the selected vertex position
	int getVertex(int x, int y);
};
void OpenGLinitial();
void PanelResize(int width , int height);
void Render_Init(int width , int height);
 #endif