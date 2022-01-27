#pragma once
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "math/Vector2.h"
#include "stdafx.h"

struct WeightData
{
	WeightData(Vector2& p, int w):pos(p), weight(w)	{}
	Vector2 pos;
	int  weight;
	bool operator<(const WeightData& wd)
	{
		return weight < wd.weight;
	}
};
typedef std::vector<WeightData> Weights;

typedef Vector2s Line;
typedef std::vector<Line> Lines;

class vavImage
{
public:
	vavImage(const cv::Mat& im);
	vavImage(void);
	~vavImage(void);
	bool	ReadImage(std::string path);
	cv::Mat	CannyEdge(double threshold1=0, double threshold2=30, int apertureSize=3, bool L2gradient=false);
	bool	Vaild();
	int	GetWidth() {return m_Image.rows;}
	int	GetHeight() {return m_Image.cols;}
	Lines	ComputeEdgeLine() { return ComputeEdgeLine(m_Image);}
	Lines	ComputeEdgeLine(const cv::Mat& image);
	void	EdgeLink(cv::Mat& image, Line& now_line);
	const cv::Vec3b& GetColor(int x, int y) const;
	void	GetFeatureEdge(Lines& lines);
	void	TPSFromFeatureEdge(const Lines& lines);
	void	bwmorph_clean();
	bool	IsZero(int x, int y);
	void	Skel();
	void	Threshold(int v);
	cv::Mat	Laplace(int aperture_size=3);
	void	ShowEdgeLine();
	void	Resize(int x, int y, int method = cv::INTER_LINEAR);
	void	drawImage();

	int IsinsidePoint(int x1, int y1,int x2,int y2,int x3,int y3);
	Vector2s GetContour();
private:
	cv::Mat	m_Image;
	cv::Mat save_Image;
};
