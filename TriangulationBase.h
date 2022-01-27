#pragma once
#include "math/Vector2.h"
#include <vector>

struct Triangle
{
	Vector2	m_Points[3];
};
typedef std::vector<Triangle> Triangles;

class TriangulationBase
{
public:
	virtual ~TriangulationBase(void){}
	virtual void AddPoint(float x, float y)=0;
	virtual void Compute()=0;
	virtual void DelaunayMesher2()=0;
	Triangles& GetTriangles()
	{
		return m_Triangles;
	}
protected:
	Triangles	m_Triangles;
};

