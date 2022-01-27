#define _USE_MATH_DEFINES

#include <windows.h>
#include<iostream>
#include <vector>
#include <cmath>
#include <time.h>
#define GLUT_DISABLE_ATEXIT_HACK
#include <GL/gl.h>
#include <GL/glu.h>

using namespace std;

#include "ArapInteractor_2.h"

// This allows you to customize the deformer - if you don't understand,
// read that Igarashi's paper first :)

int step1_only = false;
int show_fitted = true;

#define ARAP_DO_VERIFY	0

#define ARAP_DO_VERIFY	0

void TriMesh2D::draw(bool linemode)
{	
	//glColor3f(1.0, 0.0, 0.0);
	if (linemode)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, GL_DOUBLE, sizeof(vertices[0]), &vertices[0][0]);

	glDrawElements(GL_TRIANGLES, tris.size()*3, GL_UNSIGNED_INT, (GLvoid*)tris[0]);
	

}
void TriMesh2D::compute_normal()
{
	double v1[2];
	double v2[2];
	for(int i=0;i<tris.size();i++)
	{
		v1[0]=vertices[tris[i][1]][0]-vertices[tris[i][0]][0];
		v1[1]=vertices[tris[i][1]][1]-vertices[tris[i][0]][1];
		v2[0]=vertices[tris[i][2]][0]-vertices[tris[i][0]][0];
		v2[1]=vertices[tris[i][2]][1]-vertices[tris[i][0]][1];
		if(((v1[0]*v2[1])-(v2[0]*v1[1]))>=0)
			normals.push_back(1);
		else 
			normals.push_back(-1);
		//cout<<normals[i]<<endl;
	}
}
int TriMesh2D::normal_detection()
{

	double v1[2];
	double v2[2];
	double P[3][2];
	double boundry=5.0;
	if(normals.size()==tris.size())
	{
		//cout<<normals.size()<<","<<tris.size()<<endl;
		for(int i=0;i<normals.size();i++)
		{
			P[0][0]=vertices[tris[i][0]][0];	P[0][1]=vertices[tris[i][0]][1];
			P[1][0]=vertices[tris[i][1]][0];	P[1][1]=vertices[tris[i][1]][1];
			P[2][0]=vertices[tris[i][2]][0];	P[2][1]=vertices[tris[i][2]][1];

			for(int j=0;j<3;j++)//三角形3點
			{
				for(int k=0;k<4;k++)//上下左右偵測將會反向
				{
					P[j][0]=vertices[tris[i][j]][0];
					P[j][1]=vertices[tris[i][j]][1];
					switch(k)
					{
					case 0://up
						P[j][1]=vertices[tris[i][j]][1]-boundry;
						break;
					case 1://down
						P[j][1]=vertices[tris[i][j]][1]+boundry;
						break;
					case 2://left
						P[j][0]=vertices[tris[i][j]][0]-boundry;
						break;
					case 3://right
						P[j][0]=vertices[tris[i][j]][0]+boundry;
						break;
					}
					v1[0]=P[1][0]-P[0][0];
					v1[1]=P[1][1]-P[0][1];
					v2[0]=P[2][0]-P[0][0];
					v2[1]=P[2][1]-P[0][1];
					if(((v1[0]*v2[1])-(v2[0]*v1[1]))>=0)
					{
						if(normals[i]!=1) return 1;
					}
					else
						if(normals[i]!=-1) return 1;
				}
			}
		}
	}
	return 0;
}
// Helper function - draw a simple widget for dragging.
// Change it to your own version to make it faster and prettier!
void draw_circle(const Point2D& center, double r,
	double err, double z, bool linemode)
{
	if (r > 1000)
		return;
	int segs = (int)(r / err);
	segs = min(segs, 720);

	glPolygonMode(GL_FRONT_AND_BACK, linemode? GL_LINE : GL_FILL);

	if (linemode)
		glBegin(GL_LINE_LOOP);
	else
		glBegin(GL_POLYGON);
	r=10;
	for (double theta = 0; theta < 2*M_PI; theta += 2*M_PI/segs) {
		Point2D p = center + r * Vec2D(cos(theta), sin(theta));
		glVertex3d(p[0], p[1], z);
	}
	glEnd();

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}
//=====================SparseMatrix=======================//

namespace ublas = boost::numeric::ublas;
namespace umf = boost::numeric::bindings::umfpack;

SparseMatrix::SparseMatrix(int m, int n)
{
	if(n_cols==0)m_rows=0;
	m_rows=m;
	n_cols=n;
	
}

SparseMatrix::SparseMatrix(const SparseMatrix& B)
{
	m_rows = B.m_rows;  
	n_cols = B.n_cols; 
	SMatrix_data = B.SMatrix_data;  
	m_row_indices = B.m_row_indices;  
	m_col_indices = B.m_col_indices;  
	PreSave.resize(m_rows,n_cols,false);
	
}
SparseMatrix::~SparseMatrix()
{

}
void SparseMatrix::resetDim(int m, int n)
{
	if(n_cols==0)m_rows=0;
	m_rows=m;
	n_cols=n;
	m_row_indices.clear();  
	m_col_indices.clear();  
	SMatrix_data.clear();  
}
double SparseMatrix::operator ()(int i ,int j)const
{
	return get(i, j);
}
double& SparseMatrix::operator ()(int i ,int j)
{
	if( i>=0 && i<m_rows && j>=0 && j<n_cols)
	{
		m_row_indices[i].insert( j );  
		m_col_indices[j].insert( i ); 
		return SMatrix_data[std::make_pair(i, j)]; 
	}
	else
		assert("no such data");
}
SparseMatrix SparseMatrix::transpose() const  
{  
	SparseMatrix T( n_cols, m_rows ); 
	for ( mat_citer i = SMatrix_data.begin(); i!= SMatrix_data.end(); i++ )  
	{  
		T( i->first.second, i->first.first ) = get( i->first.first, i->first.second ); ;
		
	}  


	return T;  
}
SparseMatrix SparseMatrix::operator+( const SparseMatrix& B ) const  
{  
	if ( m_rows != B.m_rows || n_cols != B.n_cols )  
		assert("two matrix are not same size!!");  
	SparseMatrix C( *this );  
	for ( iset_citer i = B.m_row_indices.begin(); i != B.m_row_indices.end(); i++ )  
	{  
		for ( elt_citer j = i->second.begin(); j != i->second.end(); j++ )  
		{  
			double elt = B( i->first, *j );  
			C( i->first, *j ) += elt;  
		}  
	}  
	return C;  
}  
vector<double> SparseMatrix::operator*(const vector<double>& x) const  //0.001s
{    
	
	if(x.size()!=n_cols)
		assert("x.size()!= cols"); 
	vector<double> R;
	R.resize(m_rows);
	for (int i=0 ; i<m_rows; i++ )  
	{	
		R[i]=0;
		for (int j=0; j<n_cols ; j++ )  
		{  
			double elt = get( i, j );
			R[i] += elt * x[j];  
		}  
	}  

	return R;  
}  
double SparseMatrix::get( int i, int j ) const  
{  
	if ( i >= m_rows || j >= n_cols )  
		assert(">=rows || >=cols"); 
	mat_citer iter = SMatrix_data.find( make_pair( i, j ) );  
	if ( iter == SMatrix_data.end() )  
		return 0.0;  
	return iter->second;  
} 

void SparseMatrix::PreSaveMatrix()
{
	PreSave.resize(m_rows,n_cols,false);
	for ( mat_citer i = SMatrix_data.begin(); i!= SMatrix_data.end(); i++ ) // 0.003s
	{  
		PreSave( i->first.first, i->first.second ) = get( i->first.first, i->first.second ); ;

	}
}
vector<double> SparseMatrix::solve(const vector<double>& b)
{
	if(b.size()!=n_cols)
		assert("size not the same!!");
	vector<double> x(m_rows);


	umf::symbolic_type<double> Symbolic;
	umf::numeric_type<double> Numeric;
	umf::symbolic (PreSave, Symbolic); 
	umf::numeric (PreSave, Symbolic, Numeric); 

	umf::solve (PreSave, x, b, Numeric);
	
	return x;
}
void SparseMatrix::print()
{
	//for(int i=0;i<m_rows;i++)
	//{
	//	for(int j=0;j<n_cols;j++)
	//	{
	//		//cout<<get(i,j)<<" ";
	//	}
	//	//cout<<endl;
	//}
}
//=======================ShapeView========================//
Point2D ShapeView::window2world(const Point2D& win_pt)
{
	return win_pt;
}
double ShapeView::window2world(double win_len)
{
	return win_len*1000;
}
Point2D ShapeView::world2window(const Point2D& pos)
{
	return pos;
}
void ShapeView::refresh()
{

}
ShapeView::ButtonState ShapeView::get_button_state()
{
	if(b_state==0)return LBUTTON_DOWN;
	else if(b_state==1)return MBUTTON_DOWN;
	else if(b_state==2)return RBUTTON_DOWN;
}
void ShapeView::button_state(int state)
{
	if(state==0)b_state=0;
	else if(state==1)b_state=1;
	else if(state==2)b_state=2;
}
//========================================================//
SparseMatrix BigG, G00, G01, G10, G11, Gprime, B;
SparseMatrix BigH, H00, H01, H10, H11, Hprime, D;
vector< double** > K, F, invF;
vector< vector<double> > C;
vector< vector<Point2D> > fittedVertices;
vector<vector< vector<Point2D> >> record;
ArapInteractor::ArapInteractor(ShapeView* view, const TriMesh2D& mesh)
{
	shapeview = view;

#if ARAP_DO_VERIFY
	// Build one-triangle mesh for testing
	boundary.push_back(Point2D(-0.5, -0.4));
	boundary.push_back(Point2D(0.5, -0.4));
	boundary.push_back(Point2D(0, 0.4));

	themesh.vertices = boundary;
	themesh.tris.push_back(Tri(0, 1, 2));
#else
	themesh = mesh;
	themesh.compute_normal();//!!
	boundary = themesh.vertices;
#endif
	// Save a copy
	baseVertices = themesh.vertices;

	// Reset the states
	flags.resize(themesh.vertices.size(), 0);
	savelastFlagsPosition.resize(themesh.vertices.size(),Point2D((double)0, (double)0));//!!

	beingDragged = false;
	preCompG();
	preCompF();
	preCompH();
	//cout<<"ArapInteractor set! vertices size : "<<themesh.vertices.size()<<endl;
}

ArapInteractor::~ArapInteractor(void)
{
	for (int t = 0; t < (int)F.size(); t++) {
		delete_matrix<double>(F[t], 4);
		delete_matrix<double>(invF[t], 4);
		delete_matrix<double>(K[t], 6);
	}
}
#if 0
void ArapInteractor::getGeachTri(const Tri& t, double G[][6])
{

	// G = A'A
	double A[2][6], At[6][2];
	for (int i = 0; i < 2; i++)
		for (int j = 0; j < 6; j++)
			A[i][j] = 0;

	for (int i = 0; i < 3; i++) {
		int i0 = i, i1 = (i+1)%3, i2 = (i+2)%3;

		const Point2D& v0 = themesh.vertices[t[i0]];
		const Point2D& v1 = themesh.vertices[t[i1]];
		const Point2D& v2 = themesh.vertices[t[i2]];

		// First compute x01 and y01;
		//	xdir = v1-v0, ydir = (-xdir[1], xdir[0]), v2dir = v2-v0
		//  xdir[0] * x01 + ydir[0] * y01 = v2dir[0]
		//  xdir[1] * x01 + ydir[1] * y01 = v2dir[1]
		//

		Vec2D xdir = v1-v0;
		Vec2D ydir(-xdir[1], xdir[0]);
		Vec2D v2dir = v2-v0;
		double det = xdir[0]*ydir[1] - xdir[1]*ydir[0];
		double x01 = (v2dir[0]*ydir[1] - v2dir[1]*ydir[0])/det;
		double y01 = (xdir[0]*v2dir[1] - xdir[1]*v2dir[0])/det;

#if ARAP_DO_VERIFY
		// verify:
		Point2D v22 = v0 + x01*xdir + y01*ydir;
		double err = len(v22-v2);
#endif
		// Now ready for A

		// for v0'
		A[0][i0*2]	 += 1-x01;
		A[1][i0*2]	 += -y01;
		A[0][i0*2+1] += y01;
		A[1][i0*2+1] += 1-x01;

		// for v1'
		A[0][i1*2]	 += x01;
		A[1][i1*2]	 += y01;
		A[0][i1*2+1] += -y01;
		A[1][i1*2+1] += x01;

		// for v2'
		A[0][i2*2]	 += -1;
		A[1][i2*2]	 += 0;
		A[0][i2*2+1] += 0;
		A[1][i2*2+1] += -1;
	}

	transpose<double, 2, 6>(A, At);
	mul<double, 6, 2, 6>(At, A, G);

#if ARAP_DO_VERIFY
	print<double, 2, 6>(stdout, A);
	print<double, 6, 2>(stdout, At);

	const Point2D& v0 = themesh.vertices[t[0]];
	const Point2D& v1 = themesh.vertices[t[1]];
	const Point2D& v2 = themesh.vertices[t[2]];

	// Verify again
	double vt[] = { v0[0], v0[1], v1[0], v1[1], v2[0], v2[1] };
	double v[] =  { v0[0], v0[1], v1[0], v1[1], v2[0], v2[1] };
	double b[]  =  { 0, 0, 0, 0, 0, 0 };;
	mul<double, 6, 6>(G, v, b);
	double err = 0;
	for (int i = 0; i < 6; i++)
		err += vt[i]*b[i];

#endif

}
#else
void ArapInteractor::getGeachTri(const Tri& t, double G[][6])
{
	//std::cout << "getgeachtri" << std::endl;
	// Reset G
	for (int ii = 0; ii < 6; ii++)
		for (int jj = 0; jj < 6; jj++)
			G[ii][jj] = 0;

	// G = A'A
	double A[2][6], At[6][2], tmpG[6][6];

	for (int i = 0; i < 3; i++) {
		// Reset A
		for (int ii = 0; ii < 2; ii++)
			for (int jj = 0; jj < 6; jj++)
				A[ii][jj] = 0;

		int i0 = i, i1 = (i+1)%3, i2 = (i+2)%3;

		const Point2D& v0 = themesh.vertices[t[i0]];
		const Point2D& v1 = themesh.vertices[t[i1]];
		const Point2D& v2 = themesh.vertices[t[i2]];

		// First compute x01 and y01;
		//	xdir = v1-v0, ydir = (-xdir[1], xdir[0]), v2dir = v2-v0
		//  xdir[0] * x01 + ydir[0] * y01 = v2dir[0]
		//  xdir[1] * x01 + ydir[1] * y01 = v2dir[1]
		//

		Vec2D xdir = v1-v0;
		Vec2D ydir(-xdir[1], xdir[0]);
		Vec2D v2dir = v2-v0;
		double det = xdir[0]*ydir[1] - xdir[1]*ydir[0];
		double x01 = (v2dir[0]*ydir[1] - v2dir[1]*ydir[0])/det;
		double y01 = (xdir[0]*v2dir[1] - xdir[1]*v2dir[0])/det;

#if ARAP_DO_VERIFY
		// verify:
		Point2D v22 = v0 + x01*xdir + y01*ydir;
		double err = len(v22-v2);
#endif
		// Now ready for A

		// for v0'
		A[0][i0*2]	 = 1-x01;
		A[1][i0*2]	 = -y01;
		A[0][i0*2+1] = y01;
		A[1][i0*2+1] = 1-x01;

		// for v1'
		A[0][i1*2]	 = x01;
		A[1][i1*2]	 = y01;
		A[0][i1*2+1] = -y01;
		A[1][i1*2+1] = x01;

		// for v2'
		A[0][i2*2]	 = -1;
		A[1][i2*2]	 = 0;
		A[0][i2*2+1] = 0;
		A[1][i2*2+1] = -1;

		transpose<double, 2, 6>(A, At);
		mul<double, 6, 2, 6>(At, A, tmpG);

		// Add to G
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 6; j++)
				G[i][j] += tmpG[i][j];

#if ARAP_DO_VERIFY
		// print<double, 2, 6>(stdout, A);
		// print<double, 6, 2>(stdout, At);

		const Point2D& _v0 = themesh.vertices[t[0]];
		const Point2D& _v1 = themesh.vertices[t[1]];
		const Point2D& _v2 = themesh.vertices[t[2]];

		// Verify again
		double vt[] = { _v0[0], _v0[1], _v1[0], _v1[1], _v2[0], _v2[1] };
		double v[] =  { _v0[0], _v0[1], _v1[0], _v1[1], _v2[0], _v2[1] };
		double b[]  =  { 0, 0, 0, 0, 0, 0 };
		mul<double, 6, 6>(G, v, b);
		double err2 = 0;
		for (int i = 0; i < 6; i++)
			err2 += vt[i]*b[i];
		printf("%g\n", err2);
#endif

	}

// 	//=============================print===========================
// 	cout<<"\n\nprint G[6][6]\n";
// 	for (int i = 0; i < 6; i++)
// 	{
// 		for (int j = 0; j < 6; j++)
// 		{
// 			cout<<G[i][j]<<" ";
// 		}
// 		cout<<endl;
// 	}
// 	//=============================================================
}
#endif

void ArapInteractor::preCompG()
{
	//std::cout << "precompg" << std::endl;
	int nv = flags.size();
	BigG.resetDim(nv*2, nv*2);

	// For each triangle
	int nt = (int)themesh.tris.size();
	for (int i = 0; i < nt; i++) {
		double G[6][6];
		const Tri& t = themesh.tris[i];
		getGeachTri(t, G);

		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++) {
				BigG(t[j]*2,   t[k]*2)   += G[j*2  ][k*2  ];
				BigG(t[j]*2+1, t[k]*2)   += G[j*2+1][k*2  ];
				BigG(t[j]*2,   t[k]*2+1) += G[j*2  ][k*2+1];
				BigG(t[j]*2+1, t[k]*2+1) += G[j*2+1][k*2+1];
			}
	}
}

void ArapInteractor::preStep1()
{
	//std::cout << "prestep1" << std::endl;
	// The map from vertex index to index in the matrix (then multiply by 2)
	int nv = flags.size(), cur_free = 0, cur_ctrl = 0;
	vector<int> vert_map(nv, 0);
	for (int i = 0; i < nv; i++) {
		if (flags[i] == 0)
			vert_map[i] = cur_free++;
		else
			vert_map[i] = cur_ctrl++;
	}

	if (cur_ctrl == 0 || cur_free == 0)
		return;

	G00.resetDim(cur_free*2, cur_free*2);
	G01.resetDim(cur_free*2, cur_ctrl*2);
	G10.resetDim(cur_ctrl*2, cur_free*2);
	G11.resetDim(cur_ctrl*2, cur_ctrl*2);

	// For each triangle
	int nt = (int)themesh.tris.size();
	for (int i = 0; i < nt; i++) {
		const Tri& t = themesh.tris[i];
		Tri m(vert_map[t[0]], vert_map[t[1]], vert_map[t[2]]);

		// assign G's to G00, G01, G10, G11 accordingly
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++) {
				SparseMatrix* GG = NULL;
				// Non-zeros could only be the diagonal entries 
				if (!flags[t[j]] && !flags[t[k]])
					GG = &G00;
				else if (!flags[t[j]] && flags[t[k]])
					GG = &G01;
				else if (flags[t[j]] && !flags[t[k]])
					GG = &G10;
				else // if (flags[t[j]] && flags[t[k]])
					GG = &G11;

				(*GG)(m[j]*2,   m[k]*2)   = BigG(t[j]*2,   t[k]*2);
				(*GG)(m[j]*2+1, m[k]*2)   = BigG(t[j]*2+1, t[k]*2);
				(*GG)(m[j]*2,   m[k]*2+1) = BigG(t[j]*2,   t[k]*2+1);
				(*GG)(m[j]*2+1, m[k]*2+1) = BigG(t[j]*2+1, t[k]*2+1);
			}
	}

	SparseMatrix G00_T = G00.transpose();
	SparseMatrix G10_T = G10.transpose();

	Gprime = G00+G00_T;
	B = G01 + G10_T;
	Gprime.PreSaveMatrix();

#if ARAP_DO_VERIFY
	// Verifying
	BigG.print(stdout);
	printf("\n");

	G00.print(stdout);
	printf("\n");
	G00_T.print(stdout);
	printf("\n");
	Gprime.print(stdout);

	vector<double> v;
	vector<double> q, u;

	for (int i = 0; i < nv; i++) {
		if (!flags[i]) {
			u.push_back(themesh.vertices[i][0]);
			u.push_back(themesh.vertices[i][1]);
		} else {
			q.push_back(themesh.vertices[i][0]);
			q.push_back(themesh.vertices[i][1]);
		}
		v.push_back(themesh.vertices[i][0]);
		v.push_back(themesh.vertices[i][1]);
	}

	vector<double> BigGv = BigG * v;
	double err = 0;
	for (int i = 0; i < nv*2; i++)
		err += BigGv[i]*v[i];

	vector<double> Gpu = Gprime * u;
	vector<double> Bq = B * q;
#endif
}

void ArapInteractor::step1()
{
	//std::cout << "step1" << std::endl;

	// Compute q

	vector<double> q, Bq;

	int nv = themesh.vertices.size();;
	for (int i = 0; i < nv; i++) {
		if (!flags[i])
			continue;
		q.push_back(themesh.vertices[i][0]);
		q.push_back(themesh.vertices[i][1]);
	}
	Bq = B*q;
	for (int i = 0; i < (int)Bq.size(); i++)
		Bq[i] = -Bq[i];

	vector<double> u = Gprime.solve(Bq);//0.004s

	int ui = 0;
	for (int i = 0; i < nv; i++) {
		if (flags[i])
			continue;
		themesh.vertices[i][0] = u[ui++];
		themesh.vertices[i][1] = u[ui++];
	}
}

void ArapInteractor::getFeachTri(const Tri& t, double K[][4], double F[][4])
{
	//std::cout << "getfeachtri" << std::endl;
	// F = K'K
	double Kt[4][6];

	// Reset F
	for (int ii = 0; ii < 6; ii++)
		for (int jj = 0; jj < 4; jj++)
			K[ii][jj] = 0;
	K[0][0] = K[1][1] = K[2][2] = K[3][3] = 1;

	const Point2D& v0 = themesh.vertices[t[0]];
	const Point2D& v1 = themesh.vertices[t[1]];
	const Point2D& v2 = themesh.vertices[t[2]];

	// First compute x01 and y01;
	//	xdir = v1-v0, ydir = (-xdir[1], xdir[0]), v2dir = v2-v0
	//  xdir[0] * x01 + ydir[0] * y01 = v2dir[0]
	//  xdir[1] * x01 + ydir[1] * y01 = v2dir[1]
	//

	Vec2D xdir = v1-v0;
	Vec2D ydir(-xdir[1], xdir[0]);
	Vec2D v2dir = v2-v0;
	double det = xdir[0]*ydir[1] - xdir[1]*ydir[0];
	double x01 = (v2dir[0]*ydir[1] - v2dir[1]*ydir[0])/det;
	double y01 = (xdir[0]*v2dir[1] - xdir[1]*v2dir[0])/det;

#if ARAP_DO_VERIFY
	// verify:
	Point2D v22 = v0 + x01*xdir + y01*ydir;
	double err = len(v22-v2);
#endif
	// Now ready for F

	// for v0'
	K[4][0] = 1-x01;
	K[5][0] = -y01;
	K[4][1] = y01;
	K[5][1] = 1-x01;

	// for v1'
	K[4][2] = x01;
	K[5][2] = y01;
	K[4][3] = -y01;
	K[5][3] = x01;

	// Do transpost
	transpose<double, 6, 4>(K, Kt);
	mul<double, 4, 6, 4>(Kt, K, F);
}

void ArapInteractor::preCompF()
{
	//std::cout << "precompf" << std::endl;
	// For each triangle
	int nt = (int)themesh.tris.size();
	F.resize(nt);
	invF.resize(nt);
	K.resize(nt);

	for (int i = 0; i < nt; i++) {
		double tmpF[4][4], tmpInvF[4][4], tmpK[6][4];
		const Tri& t = themesh.tris[i];
		getFeachTri(t, tmpK, tmpF);
		copy<double, 4, 4>(tmpF, tmpInvF);
		inverse<double, 4>(tmpInvF);

		F[i] = new_matrix<double>(4, 4);
		invF[i] = new_matrix<double>(4, 4);
		K[i] = new_matrix<double>(6, 4);

		for (int ii = 0; ii < 4; ii++) {
			for (int jj = 0; jj < 4; jj++) {
				F[i][ii][jj] = tmpF[ii][jj];
				invF[i][ii][jj] = tmpInvF[ii][jj];
			}
		}

		for (int ii = 0; ii < 6; ii++)
			for (int jj = 0; jj < 4; jj++) {
				K[i][ii][jj] = tmpK[ii][jj];
			}
	}
}

// Here ti is actually not useful, just for verification purpose
void ArapInteractor::getHeachTri(const Tri& t, double H[][6])
{
	//std::cout << "getheachtri" << std::endl;
	// Reset H
	for (int ii = 0; ii < 6; ii++)
		for (int jj = 0; jj < 6; jj++)
			H[ii][jj] = 0;

	// G = tI'tI
	double tI[2][6], tIt[6][2], tmpH[6][6];

	for (int i = 0; i < 3; i++) {
		// Reset tI
		for (int ii = 0; ii < 2; ii++)
			for (int jj = 0; jj < 6; jj++)
				tI[ii][jj] = 0;

		int i0 = i, i1 = (i+1)%3, i2 = (i+2)%3;

		// for v0'
		tI[0][i0*2]	  = -1;
		tI[1][i0*2+1] = -1;

		// for v1'
		tI[0][i1*2]	  = 1;
		tI[1][i1*2+1] = 1;

		transpose<double, 2, 6>(tI, tIt);

		mul<double, 6, 2, 6>(tIt, tI, tmpH);

		// Add to H
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 6; j++)
				H[i][j] += tmpH[i][j];
	}
	//=============================print===========================
// 	cout<<"\n\nprint H[6][6]\n";
// 	for (int i = 0; i < 6; i++)
// 	{
// 		for (int j = 0; j < 6; j++)
// 		{
// 			cout<<H[i][j]<<" ";
// 		}
// 		cout<<endl;
// 	}
	//=============================================================

#if ARAP_DO_VERIFY
#endif
}

void ArapInteractor::preCompH()
{
	std::cout << "precomph" << std::endl;
	int nv = flags.size();
	BigH.resetDim(nv*2, nv*2);

	// For each triangle
	int nt = (int)themesh.tris.size();
	for (int ti = 0; ti < nt; ti++) {
		double H[6][6];
		const Tri& t = themesh.tris[ti];
		getHeachTri(t, H);
#if ARAP_DO_VERIFY
		vector<Point2D> fitted(3);
		fitted[0] = themesh.vertices[t[0]];
		fitted[1] = themesh.vertices[t[1]];
		fitted[2] = themesh.vertices[t[2]];

		double vfitted[6], v[6], f[6];
		memset(vfitted, 0, sizeof(vfitted));
		memset(v, 0, sizeof(v));
		memset(f, 0, sizeof(f));

		for (int i = 0; i < 3; i++) {
			int j = (i+1) % 3;
			Vec2D vij_f = fitted[j]-fitted[i];
			vfitted[2*i] = vij_f[0];
			vfitted[2*i+1] = vij_f[1];

			v[2*i] = themesh.vertices[t[i]][0];
			v[2*i+1] = themesh.vertices[t[i]][1];

			f[2*i] += vij_f[0];
			f[2*i+1] += vij_f[1];
			f[2*j] += -vij_f[0];
			f[2*j+1] += -vij_f[1];
		}

		double tI[6][6];
		for (int i = 0; i < 6; i++)
			for (int j = 0; j < 6; j++)
				tI[i][j] = 0;
		for (int i = 0; i < 3; i++) {
			int j = (i+1) % 3;
			tI[i*2][i*2] = -1;
			tI[i*2+1][i*2+1] = -1;
			tI[i*2][j*2] = 1;
			tI[i*2+1][j*2+1] = 1;
		}
		double vfitted2[6] = {0, 0, 0, 0, 0, 0};
		mul<double, 6, 6>(tI, v, vfitted2);

		double Hv[6];
		memset(Hv, 0, sizeof(Hv));
		mul<double, 6, 6>(H, v, Hv);
		double err = 0;
		for (int i = 0; i < 6; i++)
			err += Hv[i]*v[i];
		for (int i = 0; i < 6; i++)
			err += 2*f[i]*v[i];
		for (int i = 0; i < 6; i++)
			err += vfitted[i]*vfitted[i];
#endif
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++) {
				if (H[j*2][k*2] != 0)
					BigH(t[j]*2,   t[k]*2)   += H[j*2  ][k*2  ];
				if (H[j*2+1][k*2] != 0)
					BigH(t[j]*2+1, t[k]*2)   += H[j*2+1][k*2  ];
				if (H[j*2][k*2+1] != 0)
					BigH(t[j]*2,   t[k]*2+1) += H[j*2  ][k*2+1];
				if (H[j*2+1][k*2+1] != 0)
					BigH(t[j]*2+1, t[k]*2+1) += H[j*2+1][k*2+1];
			}
	}
}

void ArapInteractor::preStep2()//0.05
{
	//std::cout << "prestep2" << std::endl;
	// Preparing for the memory
	int nt = (int)themesh.tris.size();
	C.resize(nt);
	for (int i = 0; i < nt; i++) {
		C[i].resize(4);
	}
	fittedVertices.resize(nt);
	for (int i = 0; i < nt; i++) {
		fittedVertices[i].resize(3);
	}

	// Compute Hprime and D

	// The map from vertex index to index in the matrix (then multiply by 2)
	int nv = flags.size(), cur_free = 0, cur_ctrl = 0;
	vector<int> vert_map(nv, 0);
	for (int i = 0; i < nv; i++) {
		if (flags[i] == 0)
			vert_map[i] = cur_free++;
		else
			vert_map[i] = cur_ctrl++;
	}

	if (cur_ctrl == 0 || cur_free == 0)
		return;

	H00.resetDim(cur_free*2, cur_free*2);
	H01.resetDim(cur_free*2, cur_ctrl*2);
	H10.resetDim(cur_ctrl*2, cur_free*2);
	H11.resetDim(cur_ctrl*2, cur_ctrl*2);

	// For each triangle
	for (int i = 0; i < nt; i++) {
		const Tri& t = themesh.tris[i];
		Tri m(vert_map[t[0]], vert_map[t[1]], vert_map[t[2]]);

		// assign H's to H00, H01, H10, H11 accordingly
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++) {
				SparseMatrix* HH = NULL;
				// Non-zeros could only be the diagonal entries 
				if (!flags[t[j]] && !flags[t[k]])
					HH = &H00;
				else if (!flags[t[j]] && flags[t[k]])
					HH = &H01;
				else if (flags[t[j]] && !flags[t[k]])
					HH = &H10;
				else // if (flags[t[j]] && flags[t[k]])
					HH = &H11;

				(*HH)(m[j]*2,   m[k]*2)   = BigH(t[j]*2,   t[k]*2);
				(*HH)(m[j]*2+1, m[k]*2)   = BigH(t[j]*2+1, t[k]*2);
				(*HH)(m[j]*2,   m[k]*2+1) = BigH(t[j]*2,   t[k]*2+1);
				(*HH)(m[j]*2+1, m[k]*2+1) = BigH(t[j]*2+1, t[k]*2+1);
			}
	}

	SparseMatrix H00_T = H00.transpose();
	SparseMatrix H10_T = H10.transpose();

	Hprime = H00+H00_T;
	D = H01 + H10_T;
	Hprime.PreSaveMatrix();
}

void ArapInteractor::step2()//0.6
{
	//std::cout << "step2" << std::endl;
	// Compute C

	// For each triangle
	int nt = (int)themesh.tris.size();

	for (int i = 0; i < nt; i++) {
		double _K[6][4], _Kt[4][6];
		for (int ii = 0; ii < 6; ii++) {
			for (int jj = 0; jj < 4; jj++)
				_K[ii][jj] = K[i][ii][jj];
		}

		transpose<double, 6, 4>(_K, _Kt);

		const Tri& t = themesh.tris[i];
		// Get intermediate result
		const Point2D& v0 = themesh.vertices[t[0]];
		const Point2D& v1 = themesh.vertices[t[1]];
		const Point2D& v2 = themesh.vertices[t[2]];

		double vprime[] = { v0[0], v0[1], v1[0], v1[1], v2[0], v2[1] };
		double _C_tmp[4] = {0, 0, 0, 0};
		mul<double, 4, 6>(_Kt, vprime, _C_tmp);
		for (int j = 0; j < 4; j++)
			C[i][j] = _C_tmp[j];

		// Compute fitted triangle
		double _invF[4][4];
		for (int ii = 0; ii < 4; ii++) {
			for (int jj = 0; jj < 4; jj++)
				_invF[ii][jj] = invF[i][ii][jj];
		}
		memset(vprime, 0, sizeof(vprime));
		mul<double, 4, 4>(_invF, _C_tmp, vprime);

		// Also compute the other vertex as well
		double vfit[] = { 0, 0, 0, 0, 0, 0 };
		mul<double, 6, 4>(_K, vprime, vfit);

		fittedVertices[i][0][0] = vfit[0];
		fittedVertices[i][0][1] = vfit[1];
		fittedVertices[i][1][0] = vfit[2];
		fittedVertices[i][1][1] = vfit[3];
		fittedVertices[i][2][0] = vfit[4];
		fittedVertices[i][2][1] = vfit[5];

		// Compute scale and scale back to get congruent triangles
		// The paper doesn't say anything about how to scale, so I assue
		// it's around the gravity center
		Point2D center = fittedVertices[i][0] +
			fittedVertices[i][1] + fittedVertices[i][2];
		center /= 3.0;

		double scale = dist(fittedVertices[i][0], fittedVertices[i][1]) +
			dist(fittedVertices[i][1], fittedVertices[i][2]) + 
			dist(fittedVertices[i][2], fittedVertices[i][0]);
		scale /= dist(baseVertices[t[0]], baseVertices[t[1]]) + 
			dist(baseVertices[t[1]], baseVertices[t[2]]) + 
			dist(baseVertices[t[2]], baseVertices[t[0]]);

		for (int j = 0; j < 3; j++) {
			Vec2D v = fittedVertices[i][j]-center;
			v /= scale;
			fittedVertices[i][j] = center + v;
		}
	}

	if (step1_only)
		return;
	// The second sub-step of step2

	// Compute q
	vector<double> q, Dq_plus_f0;
	int nv = themesh.vertices.size(), cur_free = 0, cur_ctrl = 0;

	vector<int> vert_map(nv, 0);
	for (int i = 0; i < nv; i++) {
		if (flags[i] == 0)
			vert_map[i] = cur_free++;
		else {
			vert_map[i] = cur_ctrl++;
			q.push_back(themesh.vertices[i][0]);
			q.push_back(themesh.vertices[i][1]);
		}
	}

	if (cur_ctrl == 0 || cur_free == 0)
		return;

	// This is Dq

 	Dq_plus_f0 = D*q;

	// Compute f0 and f1
 	vector<double> f(nv*2, 0), f0(cur_free*2), f1(cur_ctrl*2);
 
	for (int ti = 0; ti < nt; ti++) {
		const Tri& t = themesh.tris[ti];
		vector<Point2D>& fitted = fittedVertices[ti];
		for (int i = 0; i < 3; i++) {
			int j = (i+1) % 3;
			Vec2D vij_f = fitted[j]-fitted[i];
			f[2*t[i]] += -2*vij_f[0];
			f[2*t[i]+1] += -2*vij_f[1];
			f[2*t[j]] += 2*vij_f[0];
			f[2*t[j]+1] += 2*vij_f[1];
		}
	}
	// Map them into f0, f1
	for (int i = 0; i < nv; i++) {
		if (flags[i] == 0) {
			f0[2*vert_map[i]] = f[2*i];
			f0[2*vert_map[i]+1] = f[2*i+1];
		} else {
			f1[2*vert_map[i]] = f[2*i];
			f1[2*vert_map[i]+1] = f[2*i+1];
		}
	}
	// Compute Dq_plus_f0, and negate it
	for (int i = 0; i < (int)Dq_plus_f0.size(); i++) {
		Dq_plus_f0[i] -= f0[i];
		Dq_plus_f0[i] = -Dq_plus_f0[i];
	}

	vector<double> u = Hprime.solve(Dq_plus_f0);

	int ui = 0;
	for (int i = 0; i < nv; i++) {
		if (flags[i])
			continue;
		themesh.vertices[i][0] = u[ui++];
		themesh.vertices[i][1] = u[ui++];
	}
	
	themesh.record_v.push_back(themesh.vertices);
	for (int i = 0; i < themesh.record_v.size(); i++)
	{
		for (int j = 0; j < themesh.record_v[i].size(); j++)
		{
			//std::cout << "control: " << (themesh.record_v)[i][j][0] << std::endl;
			//std::cout << "control: " << (themesh.record_v)[i][j][1] << std::endl;
		}
	}
}

void ArapInteractor::deform()
{
	//std::cout << "deform"<<std::endl;
	// The map from vertex index to index in the matrix (then multiply by 2)
	int nv = flags.size();
	vector<int> ctrl_verts;
	for (int i = 0; i < nv; i++) {
		if (flags[i])
			ctrl_verts.push_back(i);
	}

	if (ctrl_verts.size() >= 2 && (int)ctrl_verts.size() < nv) {

		step1();//0.004s
		step2();//0.004s
		return;
	}

	if (ctrl_verts.size() == nv) 
	{
	} 
	else if (ctrl_verts.size() == 1) 
	{
		Vec2D v = themesh.vertices[ctrl_verts[0]]-baseVertices[ctrl_verts[0]];
		for (int i = 0; i < nv; i++) 
		{
			themesh.vertices[i] = baseVertices[i] + v;
			//std::cout << "themesh.vertices deform: " << themesh.vertices[i] << std::endl;
		}
	} else 
	{
		themesh.vertices = baseVertices;
		
	}
	step2();//0.004

}
void ArapInteractor::recover_mesh()
{
	themesh = tempmesh;
}
void ArapInteractor::save_mesh()
{
	
	record_mesh.push_back(themesh);
}
void ArapInteractor::return_mesh(int i)
{
	tempmesh = record_mesh[0];
	if((record_mesh.size() - 1 )>i)
		themesh = record_mesh[i];
	/*cv::VideoCapture capture(0);
	cv::VideoWriter writer("VideoTest.avi", CV_FOURCC('M', 'J', 'P', 'G'), 25.0, cv::Size(640, 480));
	cv::Mat frame;

	while (capture.isOpened())
	{
		capture >> frame;
		writer << frame;
		imshow("video", frame);
		if (cvWaitKey(20) == 27) {
			break;
		}
	}*/
	/*cv::VideoCapture capture(0);
	if (!capture.isOpened()) {
		cout<<"false"<<endl;
	}
	cv::Size videoSize = cv::Size((int)capture.get(CV_CAP_PROP_FRAME_WIDTH), (int)capture.get(CV_CAP_PROP_FRAME_HEIGHT));
	cv::VideoWriter writer;
	writer.open("VideoTest.avi", CV_FOURCC('M', 'J', 'P', 'G'), 30, videoSize);
	cv::namedWindow("show image", 0);

	while (true) {
		cv::Mat frame;
		capture >> frame;
		if (!frame.empty()) {
			writer.write(frame);
			imshow("show image", frame);
			if (cv::waitKey(33) == 27) {
				break;
			}
		}
	}*/
}
void ArapInteractor::OnDraw(int vp)
{
	glTranslated(0, 0, .1);

	// Draw the mesh
	glLineWidth(1);

	glColor3f(.9f, .7f, .5f);

	glTranslated(0, 0, -.01);
	themesh.draw(false);
	glTranslated(0, 0, .01);

	glColor3f(.1f, .1f, .1f);
	themesh.draw(true);
	
	if (show_fitted && step1_only) {
		// Draw fitted triangles
		glTranslated(0, 0, .01);
		//std::cout <<"fuck" <<std::endl;
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
		glColor4f(.5f, .5f, .5f, 0.7f); 

		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

		int nt = fittedVertices.size();

		glBegin(GL_TRIANGLES);
		for (int i = 0; i < nt; i++) {
			glVertex2dv(fittedVertices[i][0]);
			glVertex2dv(fittedVertices[i][1]);
			glVertex2dv(fittedVertices[i][2]);
		}
		glEnd();

		glDisable(GL_BLEND);
		glTranslated(0, 0, -.01);
	}

	// Draw the control points
	glTranslated(0, 0, .1);
	glLineWidth(2);
	glColor3f(.0f, .0f, .8f);
	glEnable(GL_LINE_SMOOTH);
	for (int i = 0; i < (int)themesh.vertices.size(); i++){
		if (flags[i] > 0) {
			draw_circle(themesh.vertices[i], 0.02, 0.002, 0, true);
		}
	}
	glDisable(GL_LINE_SMOOTH);
	glTranslated(0, 0, -.1);

	glTranslated(0, 0, -.1);

}

void ArapInteractor::OnDisplay()
{
}

void ArapInteractor::OnDrawInfo(int vp)
{
}

void ArapInteractor::OnMotion(int x, int y, int flag,bool mouse_down, int vp)
{
	//std::cout << "onmotion" << std::endl;
	Point2D pos = Point2D((double)x, (double)y);
	
	double thresh = 15.0;
	int detection;//偵測normal反向
	if (mouse_down) {
		beingDragged = true;
		// setup constraint
		themesh.vertices[flag] = pos;
		deform();

		if(!themesh.normal_detection())//normal不會反向
		{
			for(int i=0;i<themesh.vertices.size();i++)
			{
				savelastFlagsPosition[i]=themesh.vertices[i];
				//std::cout << "savelastFlagsPosition: "<< savelastFlagsPosition [i]<<std::endl;
			}			
		}
		else
		{
			for(int i=0;i<themesh.vertices.size();i++)
			{
				themesh.vertices[i]=savelastFlagsPosition[i];
				//std::cout << "themesh.vertices[i]: " << themesh.vertices[i] << std::endl;
			}
		}
	}
}
int ArapInteractor::getVertex(int x, int y)
{	
	Point2D pos = Point2D((double)x, (double)y);
	double thresh = 15.0;
	for (int i = 0; i < (int)themesh.vertices.size(); i++){
		if (dist(themesh.vertices[i], pos) > thresh)
			continue;

		if (flags[i] == 1)
			return i;
	}
	return -1;
}
void ArapInteractor::OnMouse(int button, int button_state, int x, int y, int vp)
{
	Point2D pos = Point2D((double)x, (double)y);
 	double thresh = 15.0;
 
	
	if (button_state == GLUT_UP && button ==  GLUT_LEFT_BUTTON) //0.15s
	{   
		// left click to pick the approximated vertex
		
		for (int i = 0; i < (int)themesh.vertices.size(); i++){
			if (dist(themesh.vertices[i], pos) > thresh)
				continue;

			if (flags[i] == 0)
				flags[i] = 1;
			else
				break;
				
			preStep1();
 			preStep2();

			break;
		}
		ctrlPoints.clear();
		beingDragged = false;

		// Prepare for dragging

		for (int i = 0; i < (int)themesh.vertices.size(); i++){
			if (dist(themesh.vertices[i], pos) > thresh || !flags[i])
				continue;

			ctrlPoints.push_back(i);
		}
	}
	if(button_state == GLUT_UP && button ==  GLUT_RIGHT_BUTTON)
	{
		// right click to delete the control point

		for (int i = 0; i < (int)themesh.vertices.size(); i++){
			if (dist(themesh.vertices[i], pos) > thresh)
				continue;

			if (flags[i] == 1)// && !beingDragged)
				flags[i] = 0;
			else
				break;

			preStep1();
			preStep2();

			break;
		}
		ctrlPoints.clear();
		beingDragged = false;
	}
}

void ArapInteractor::OnKeyboard(unsigned char key, int x, int y)
{
	if (key == '1')
		step1_only = !step1_only;
	if (key == 'f')
		show_fitted = !show_fitted;
}

void ArapInteractor::OnSpecialKey(int key, int x, int y)
{
}
void OpenGLinitial()
{
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	glDisable(GL_LIGHTING);
	glEnable(GL_TEXTURE_2D);
	glEnable(GL_COLOR_MATERIAL);

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);		// Set Texture Max Filter
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);	    // Set Texture Min Filter

	GLfloat  ambientLight[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat  diffuseLight[] = { 0.8f, 0.5f, 0.0f, 1.0f };
	GLfloat  specular[] = { 0.6f, 0.6f, 0.6f, 1.0f};
	GLfloat  Lposition[] = {0.0f, 10.0f, 0.0f , 1.0f};

	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glLightfv(GL_LIGHT0,GL_AMBIENT,ambientLight);
	glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseLight);
	glLightfv(GL_LIGHT0,GL_SPECULAR,specular);
	glLightfv(GL_LIGHT0,GL_POSITION,Lposition);
	glEnable(GL_LIGHT0);

	glClearColor(1.0f , 1.0f , 1.0f , 1.0f);
}
void PanelResize(int width , int height)
{
	if (height == 0)										    // Prevent A Divide By Zero By
		height = 1;	    									    // Making Height Equal One

	glViewport(0, 0, width, height);								    // Reset The Current Viewport

	glMatrixMode(GL_PROJECTION);									    // Select The Projection Matrix
	glLoadIdentity();									            // Reset The Projection Matrix

	gluOrtho2D(0, width, height, 0);

	glMatrixMode(GL_MODELVIEW);							                    // Select The Modelview Matrix
	glLoadIdentity();
}
void Render_Init(int width , int height)
{

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();

	glColor3f(1.0, 1.0, 1.0);
	glPushMatrix();


	glDisable( GL_DEPTH_TEST ) ;
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glEnable(GL_COLOR_MATERIAL);
          
	glClear(GL_COLOR_BUFFER_BIT);

	//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//glEnable(GL_DOUBLEBUFFER);
	//glDisable(GL_BLEND);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, width, height, 0);					// 平行六面體即View Volume =（左下x座標，右上x座標, 左下y座標,右上y座標)
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glViewport(0, 0, width, height);					// 窗口大小 =（左下x座標，左下y座標, 寬, 高)
	glTranslated(50,50,0);

}