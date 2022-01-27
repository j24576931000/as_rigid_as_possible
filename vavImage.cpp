#include "stdafx.h"
#include <auto_link_opencv.hpp>
#include <algorithm>
#include "vavImage.h"

#include <windows.h>
#define GLUT_DISABLE_ATEXIT_HACK
#include <GL/gl.h>
#include <GL/glu.h>

#include <iostream>
using namespace std;

vavImage::vavImage(void)
{
}

vavImage::vavImage( const cv::Mat& im )
{
	m_Image = im;
}

vavImage::~vavImage(void)
{
}

bool vavImage::ReadImage( std::string path )
{
	m_Image = cv::imread(path.c_str());
	if (m_Image.cols > 0 && m_Image.rows > 0)
		return true;
	else 
		return false;
}

cv::Mat vavImage::CannyEdge( double threshold1/*=0*/, double threshold2/*=30*/, int apertureSize/*=3*/, bool L2gradient/*=false*/ )
{
	cv::Mat edges;
	cvtColor(m_Image, edges, CV_BGR2GRAY); //convert from RGB color space to GRAY
	Canny(edges,edges,
		threshold1,
		threshold2,
		apertureSize, L2gradient);
	cvtColor(edges, edges, CV_GRAY2BGR);
	bwmorph_clean();

	return edges;
}

bool vavImage::Vaild()
{
	if (m_Image.cols >0 && m_Image.rows > 0)
		return true;
	return false;
}

int vavImage::IsinsidePoint(int x1, int y1,int x2,int y2,int x3,int y3)
{
	int x;
	int y;
	x=(x1+x2)/2;
	y=(y1+y2)/2;
	x=(x+x3)/2;
	y=(y+y3)/2;
	//std::cout << x << ", " << y << std::endl;

 	cv::Vec3b intensity = m_Image.at<cv::Vec3b>((int)y,(int) x);

	if (intensity[0] == 255 && intensity[1] == 255 && intensity[2] == 255)
	{
		return 1;
	}
	return 0;
	
}
Vector2s vavImage::GetContour()
{
	Vector2s out;
	cv::Mat m_Image_temp;
	vector<vector<cv::Point>> contour;
	floodFill(m_Image,cv::Point2i(m_Image.cols/2, m_Image.rows/2),cv::Scalar(255,255,255));

	//**find Contours**//

	cvtColor(m_Image, m_Image_temp, CV_BGR2GRAY);
	findContours(m_Image_temp, contour, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0,0));

	for(int i=0;i<contour[0].size();i++)
	{
		out.push_back(Vector2(contour[0][i].x,contour[0][i].y));
	}
	return out;
}

Lines vavImage::ComputeEdgeLine(const cv::Mat& image)
{
	cv::Mat tImage = image;
	Lines	res;
	for (int j=0;j < tImage.cols;++j)
	{
		for (int i=0;i < tImage.rows;++i)
		{
			cv::Vec3b& intensity = tImage.at<cv::Vec3b>(i, j);
			if (intensity[0] != 0 || intensity[1] != 0 || intensity[2] != 0)
			{
				Line line;
				line.push_back(Vector2(i, j));
				intensity[0] = 0;
				intensity[1] = 0;
				intensity[2] = 0;
				EdgeLink(tImage, line);
				res.push_back(line);
			}
		}
	}
	return res;
}

Weights wm_init;
Weights wm_init2;
struct _tmp_function
{
	_tmp_function()
	{
		wm_init.push_back(WeightData(Vector2(-1,-1), 1));
		wm_init.push_back(WeightData(Vector2( 0,-1), 1));
		wm_init.push_back(WeightData(Vector2( 1,-1), 1));
		wm_init.push_back(WeightData(Vector2(-1, 0), 1));
		wm_init.push_back(WeightData(Vector2( 1, 0), 1));
		wm_init.push_back(WeightData(Vector2(-1, 1), 1));
		wm_init.push_back(WeightData(Vector2( 0, 1), 1));
		wm_init.push_back(WeightData(Vector2( 1, 1), 1));

		wm_init2.push_back(WeightData(Vector2(-2,-2), 1));
		wm_init2.push_back(WeightData(Vector2(-1,-2), 1));
		wm_init2.push_back(WeightData(Vector2( 0,-2), 1));
		wm_init2.push_back(WeightData(Vector2( 1,-2), 1));
		wm_init2.push_back(WeightData(Vector2( 2,-2), 1));
		
		wm_init2.push_back(WeightData(Vector2(-2,-1), 1));
		wm_init2.push_back(WeightData(Vector2( 2,-1), 1));

		wm_init2.push_back(WeightData(Vector2(-2, 0), 1));
		wm_init2.push_back(WeightData(Vector2( 2, 0), 1));

		wm_init2.push_back(WeightData(Vector2(-2, 1), 1));
		wm_init2.push_back(WeightData(Vector2( 2, 1), 1));

		wm_init2.push_back(WeightData(Vector2(-2, 2), 1));
		wm_init2.push_back(WeightData(Vector2(-1, 2), 1));
		wm_init2.push_back(WeightData(Vector2( 0, 2), 1));
		wm_init2.push_back(WeightData(Vector2( 1, 2), 1));
		wm_init2.push_back(WeightData(Vector2( 2, 2), 1));
	}
}__tmp_function;

void vavImage::EdgeLink( cv::Mat& image, Line& now_line )
{
	bool	edgefail = false;
	for (;!edgefail;)
	{
		edgefail = true;
		Weights wm = wm_init;
		if (now_line.size() > 1)
		{
			int x, y;
			Vector2 move = now_line.back() - *(now_line.end()-2);
			for (int i=0;i<wm.size();i++)
			{
				if (move.y != 0 && move.y == wm[i].pos.y)
				{
					wm[i].weight++;
				}
				if (move.x != 0 && move.x == wm[i].pos.x)
				{
					wm[i].weight++;
				}
				if (wm[i].pos == move)
				{
					wm[i].weight++;
				}
			}
		}
		std::sort(wm.begin(), wm.end());
		const Vector2& v = now_line.back();
		for (int i=0;i<wm.size();i++)
		{
			int x = v.x+wm[i].pos.x;
			int y = v.y+wm[i].pos.y;
			if (x<0)x=0;
			if (x>=image.rows)x=image.rows-1;
			if (y<0)y=0;
			if (y>=image.cols)y=image.cols-1;
			cv::Vec3b& intensity = image.at<cv::Vec3b>(x, y);
			if (intensity[0] != 0 || intensity[1] != 0 || intensity[2] != 0)
			{
				now_line.push_back(Vector2(x, y));
				intensity[0] = 0;
				intensity[1] = 0;
				intensity[2] = 0;
				edgefail = false;
				break;
			}
		}
		if (edgefail)
		{
			for (int i=0;i<wm_init2.size();i++)
			{
				int x = v.x+wm_init2[i].pos.x;
				int y = v.y+wm_init2[i].pos.y;
				if (x<0)x=0;
				if (x>=image.rows)x=image.rows-1;
				if (y<0)y=0;
				if (y>=image.cols)y=image.cols-1;
				cv::Vec3b& intensity = image.at<cv::Vec3b>(x, y);
				if (intensity[0] != 0 || intensity[1] != 0 || intensity[2] != 0)
				{
					//std::cout<<"x:"<<x<<", y:"<<y<<std::endl;
					now_line.push_back(Vector2(x, y));
					intensity[0] = 0;
					intensity[1] = 0;
					intensity[2] = 0;
					edgefail = false;
					break;
				}
			}
		}
	}
	std::reverse(now_line.begin(), now_line.end());
	for (;!edgefail;)
	{
		edgefail = true;
		Weights wm = wm_init;
		if (now_line.size() > 1)
		{
			int x, y;
			Vector2 move = now_line.back() - *(now_line.end()-2);
			for (int i=0;i<wm.size();i++)
			{
				if (move.y != 0 && move.y == wm[i].pos.y)
				{
					wm[i].weight++;
				}
				if (move.x != 0 && move.x == wm[i].pos.x)
				{
					wm[i].weight++;
				}
				if (wm[i].pos == move)
				{
					wm[i].weight++;
				}
			}
		}
		std::sort(wm.begin(), wm.end());
		const Vector2& v = now_line.back();
		for (int i=0;i<wm.size();i++)
		{
			int x = v.x+wm[i].pos.x;
			int y = v.y+wm[i].pos.y;
			if (x<0)x=0;
			if (x>=image.rows)x=image.rows-1;
			if (y<0)y=0;
			if (y>=image.cols)y=image.cols-1;
			cv::Vec3b& intensity = image.at<cv::Vec3b>(x, y);
			if (intensity[0] != 0 || intensity[1] != 0 || intensity[2] != 0)
			{
				now_line.push_back(Vector2(x, y));
				intensity[0] = 0;
				intensity[1] = 0;
				intensity[2] = 0;
				edgefail = false;
				break;
			}
		}
		if (edgefail)
		{
			for (int i=0;i<wm_init2.size();i++)
			{
				int x = v.x+wm_init2[i].pos.x;
				int y = v.y+wm_init2[i].pos.y;
				if (x<0)x=0;
				if (x>=image.rows)x=image.rows-1;
				if (y<0)y=0;
				if (y>=image.cols)y=image.cols-1;
				cv::Vec3b& intensity = image.at<cv::Vec3b>(x, y);
				if (intensity[0] != 0 || intensity[1] != 0 || intensity[2] != 0)
				{
					//std::cout<<"x:"<<x<<", y:"<<y<<std::endl;
					now_line.push_back(Vector2(x, y));
					intensity[0] = 0;
					intensity[1] = 0;
					intensity[2] = 0;
					edgefail = false;
					break;
				}
			}
		}
	}
}

const cv::Vec3b& vavImage::GetColor( int x, int y ) const
{
	if (x<0) x=0;
	if (x>=m_Image.rows) x=m_Image.rows-1;
	if (y<0) y=0;
	if (y>=m_Image.cols) y=m_Image.cols-1;
	//cv::Vec3b a;a=m_Image.at<cv::Vec3b>(x, y);
	//cout<<(int)a[0]<<endl;
	return m_Image.at<cv::Vec3b>(x, y);
}

void vavImage::GetFeatureEdge( Lines& lines )
{
	return;
	for (Lines::iterator it = lines.begin(); it!=lines.end();++it)
	{
		Line& li = *it;
		if (li.size() == 3)
		{
			li.erase(li.begin()+1);
		}
		else if (li.size()>3)
		{
			Vector2 lastmove = li[1] - li[0];
			for (int i=2;i<li.size()-1;++i)
			{
				Vector2 move = li[i+1] - li[i];
				if (move == lastmove)
				{
					li[i-1].x = -999;
				}
				lastmove = move;
			}
			for (int i=1;i<li.size();++i)
			{
				if (li[i].x == -999)
				{
					li.erase(li.begin()+i);
					i--;
				}
			}
		}
		
	}
}

void vavImage::bwmorph_clean()
{
	for (int j=1;j < m_Image.cols-1;++j)
	{
		for (int i=1;i < m_Image.rows-1;++i)
		{
			cv::Vec3b& intensity = m_Image.at<cv::Vec3b>(i, j);
			if (intensity[0] != 0 && intensity[1] != 0 && intensity[2] != 0)
			{
				bool zero = true;
				zero &= IsZero(i-1, j-1);
				zero &= IsZero(i-1, j  );
				zero &= IsZero(i-1, j+1);
				zero &= IsZero(i  , j-1);
				zero &= IsZero(i  , j+1);
				zero &= IsZero(i+1, j-1);
				zero &= IsZero(i+1, j  );
				zero &= IsZero(i+1, j+1);
				if (zero)
				{
					intensity[0] = 0;
					intensity[1] = 0;
					intensity[2] = 0;
				}
			}
		}
	}
}

bool vavImage::IsZero( int x, int y )
{
	cv::Vec3b& intensity = m_Image.at<cv::Vec3b>(x, y);
	if (intensity[0] == 0 && intensity[1] == 0 && intensity[2] == 0)
		return true;
	return false;
}

void vavImage::Skel()
{
	cv::Mat skel(m_Image.size(), m_Image.type(), cv::Scalar(0));
	cv::Mat temp;
	cv::Mat eroded;
	cv::Mat element = cv::getStructuringElement(cv::MORPH_CROSS, cv::Size(3, 3));

	bool done;
	do
	{
		cv::erode(m_Image, eroded, element);
		cv::dilate(eroded, temp, element); // temp = open(img)
		cv::subtract(m_Image, temp, temp);
		cv::bitwise_or(skel, temp, skel);
		eroded.copyTo(m_Image);

		done = (cv::norm(m_Image) == 0);
	} while (!done);
	m_Image = skel;
}

void vavImage::Threshold( int v )
{
	v *= 3;
	for (int j=1;j < m_Image.cols-1;++j)
	{
		for (int i=1;i < m_Image.rows-1;++i)
		{
			cv::Vec3b& intensity = m_Image.at<cv::Vec3b>(i, j);
			if ((intensity[0] + intensity[1] + intensity[2]) >= v)
			{
				intensity[0] = 255;
				intensity[1] = 255;
				intensity[2] = 255;
			}
			else
			{
				intensity[0] = 0;
				intensity[1] = 0;
				intensity[2] = 0;
			}
		}
	}
}

cv::Mat vavImage::Laplace( int aperture_size/*=3*/ )
{
	cv::Mat edges;
	cvtColor(m_Image, edges, CV_BGR2GRAY); //convert from RGB color space to GRAY
	Laplacian(edges, edges, CV_8U, aperture_size, 1);
	cvtColor(edges, edges, CV_GRAY2BGR);
	return edges;
}

void vavImage::ShowEdgeLine()
{
	Lines lines = ComputeEdgeLine(m_Image);
	for (Lines::iterator it = lines.begin(); it!=lines.end();++it)
	{
		for (Line::iterator it2 = it->begin(); it2!=it->end();++it2)
		{
			cv::Vec3b& intensity = m_Image.at<cv::Vec3b>(it2->x, it2->y);
			intensity[0] = 0;
			intensity[1] = 0;
			intensity[2] = 0;
		}
	}
	GetFeatureEdge(lines);
	for (Lines::iterator it = lines.begin(); it!=lines.end();++it)
	{
		int R = rand()%256, B = rand()%256, G = rand()%256;
		for (Line::iterator it2 = it->begin(); it2!=it->end();++it2)
		{
			cv::Vec3b& intensity = m_Image.at<cv::Vec3b>(it2->x, it2->y);
			intensity[0] = R;
			intensity[1] = G;
			intensity[2] = B;
		}
	}
}

void vavImage::Resize( int x, int y, int method )
{
	cv::Mat tmp = m_Image;
	resize(tmp, m_Image, cv::Size(x, y), 0, 0, method);
}
void vavImage::drawImage()
{	
	glPointSize(1.0);
	glBegin(GL_POINTS);
	for (int j=1;j < m_Image.cols-1;++j)
	{
		for (int i=1;i < m_Image.rows-1;++i)
		{
			cv::Vec3b& intensity = m_Image.at<cv::Vec3b>(i, j);
			glColor3d(intensity[2], intensity[1], intensity[0]);
			glVertex2d(j,i);
		}
	}
	glEnd();
}
