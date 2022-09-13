#include <iostream>
#include "tgaimage.h"
#include <cmath>
#include <string>
#include <vector>



#define Width 2048
#define Height 2048
using namespace std;
const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
const double PI = 4 * atan(1);

void DrawLine(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color)
{
	if (x0 > x1) std::swap(x0, x1);
	if (y0 > y1) std::swap(y0, y1);
	for (int i = x0; i <= x1; i++) image.set(i, y0, color);
	for (int i = y0; i <= y1; i++) image.set(x0, i, color);
}
void line_Bresenham(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color)
{
	double dx = abs(x1 - x0), dy = abs(y1 - y0);
	bool sw = false;
	if (dx == 0 || dy == 0)
	{
		DrawLine(x0, y0, x1, y1, image, color);
		return;
	}
	if (dy > dx)
	{
		std::swap(x0, y0);
		std::swap(x1, y1);
		sw = true;
	}
	if (x0 > x1)
	{
		std::swap(x0, x1);
		std::swap(y0, y1);
	}
	for (int x = x0; x <= x1; x++)
	{
		float t = (x - x0) / (float)(x1 - x0);
		int y = round(y0 * (1. - t) + y1 * t);
		if (sw) image.set(y, x, color);
		else image.set(x, y, color);
	}
}


struct Matrix
{
	float T[4][4];
};



struct O_Vertex
{
	float x;
	float y;
	float z;
	float w;

	O_Vertex(float X = 0, float Y = 0, float Z = 0, float W = 1)
	{
		x = X;
		y = Y;
		z = Z;
		w = W;
	}
	O_Vertex operator +(const O_Vertex& v) const { return O_Vertex(x + v.x, y + v.y, z + v.z); }
	O_Vertex operator -(const O_Vertex& v) const { return O_Vertex(x - v.x, y - v.y, z - v.z); }
	O_Vertex operator ^(const O_Vertex& v) const { return O_Vertex(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
	float operator *(const O_Vertex& v) const { return x * v.x + y * v.y + z * v.z; }
	O_Vertex operator *(float f) const { return O_Vertex(x * f, y * f, z * f); }
	float norm() const { return std::sqrt(x * x + y * y + z * z); }
	O_Vertex& normalize() { *this = (*this) * (1 / (*this).norm()); return *this; }
};
TGAColor Diffuse(TGAImage& image, O_Vertex uv)
{
	return image.get(uv.x, uv.y);
}
O_Vertex Normal(O_Vertex v0, O_Vertex v1, O_Vertex v2)
{
	O_Vertex rez, t, v;

	t.x = v2.x - v0.x;
	t.y = v2.y - v0.y;
	t.z = v2.z - v0.z;

	v.x = v1.x - v0.x;
	v.y = v1.y - v0.y;
	v.z = v1.z - v0.z;

	rez.x = t.y * v.z - t.z * v.y;
	rez.y = t.z * v.x - t.x * v.z;
	rez.z = t.x * v.y - t.y * v.x;

	float norm = sqrt(rez.x * rez.x + rez.y * rez.y + rez.z * rez.z);
	O_Vertex r = rez;
	rez.x = r.x / norm;
	rez.y = r.y / norm;
	rez.y = r.y / norm;

	return rez;
}

float Scalar(O_Vertex t, O_Vertex v)
{
	return t.x * v.x + t.y * v.y + t.z * v.z;
}



void triangle_Z(O_Vertex t0, O_Vertex t1, O_Vertex t2, O_Vertex uv0, O_Vertex uv1, O_Vertex uv2, TGAImage& image, TGAImage& diffuse, float intensity, int* zbuffer)
{
	if (t0.y == t1.y && t0.y == t2.y) return;
	O_Vertex v0, v1, uva, uvb;
	if (t0.y > t1.y) swap(t0, t1);
	if (t0.y > t2.y) swap(t0, t2);
	if (t1.y > t2.y) swap(t1, t2);
	int height = t2.y - t0.y;
	for (int i = 0; i < height; i++)
	{
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
		v0 = t0 + (t2 - t0) * alpha;
		v1 = second_half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
		uva = uv0 + (uv2 - uv0) * alpha;
		uvb = second_half ? uv1 + (uv2 - uv1) * beta : uv0 + (uv1 - uv0) * beta;

		if (v0.x > v1.x) { swap(v0, v1); swap(uva, uvb); }
		
		for (int j = v0.x; j <= v1.x; j++) {
			float phi = v1.x == v0.x ? 1. : (float)(j - v0.x) / (float)(v1.x - v0.x);
			O_Vertex p, uvP;
			p.x = (int)(v0.x + (v1.x - v0.x) * phi);
			p.y = (int)(v0.y + (v1.y - v0.y) * phi);
			p.z = (int)(v0.z + (v1.z - v0.z) * phi);
			uvP.x = (int)(uva.x + (uvb.x - uva.x) * phi);
			uvP.y = (int)(uva.y + (uvb.y - uva.y) * phi);
			
			int idx = p.x + p.y * Width;
			if (zbuffer[idx] < p.z) {
				zbuffer[idx] = p.z;
				TGAColor color = Diffuse(diffuse, uvP);
				image.set(p.x, p.y, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity, 255));
				
			}
		}
	}
}void triangle_Z(O_Vertex t0, O_Vertex t1, O_Vertex t2, TGAImage& image, TGAColor color, int* zbuffer, float &z_max, float &z_min)
{
	if (t0.y == t1.y && t0.y == t2.y) return;
	O_Vertex v0, v1;
	if (t0.y > t1.y) swap(t0, t1);
	if (t0.y > t2.y) swap(t0, t2);
	if (t1.y > t2.y) swap(t1, t2);
	int height = t2.y - t0.y;
	for (int i = 0; i < height; i++)
	{
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height;
		
		v0.x = t0.x + float(t2.x - t0.x) * alpha;
		v0.y = t0.y + float(t2.y - t0.y) * alpha;
		v0.z = t0.z+ float(t2.z - t0.z) * alpha;
		v1.x = second_half ? t1.x + float(t2.x - t1.x) * beta : t0.x + float(t1.x - t0.x) * beta;
		v1.y = second_half ? t1.y + float(t2.y - t1.y) * beta : t0.y + float(t1.y - t0.y) * beta;
		v1.z = second_half ? t1.z + float(t2.z - t1.z) * beta : t0.z + float(t1.z - t0.z) * beta;

		if (v0.x > v1.x) swap(v0, v1);
		
		for (int j = v0.x; j <= v1.x; j++) {
			float phi = v1.x == v0.x ? 1. : (float)(j - v0.x) / (float)(v1.x - v0.x);
			
			O_Vertex p;
			p.x = (int)(v0.x + (v1.x - v0.x) * phi);
			p.y = (int)(v0.y + (v1.y - v0.y) * phi);
			p.z = (int)(v0.z + (v1.z - v0.z) * phi);
			int idx = p.x + p.y * Width;
			if (zbuffer[idx] < p.z) {
				zbuffer[idx] = p.z;
				image.set(p.x, p.y, color);
				if (p.z > z_max)
					z_max = p.z;
				if (p.z < z_min)
					z_min = p.z;
			}
		}
	}
}

void Multiplication(O_Vertex& v, float T[4][4])
{
	float x0 = v.x, y0 = v.y, z0 = v.z, w0 = v.w;
	v.w = x0 * T[0][3] + y0 * T[1][3] + z0 * T[2][3] + w0 * T[3][3];
	v.x = round((x0 * T[0][0] + y0 * T[1][0] + z0 * T[2][0] + w0 * T[3][0]) / v.w);
	v.y = round((x0 * T[0][1] + y0 * T[1][1] + z0 * T[2][1] + w0 * T[3][1]) / v.w);
	v.z = round((x0 * T[0][2] + y0 * T[1][2] + z0 * T[2][2] + w0 * T[3][2]) / v.w);
}
void Parallel(Matrix& t, float tx, float ty, float tz)
{
	Matrix t0 = t;

	t.T[0][0] = t0.T[0][0] + t0.T[0][3] * tx;
	t.T[0][1] = t0.T[0][1] + t0.T[0][3] * ty;
	t.T[0][2] = t0.T[0][2] + t0.T[0][3] * tz;

	t.T[1][0] = t0.T[1][0] + t0.T[1][3] * tx;
	t.T[1][1] = t0.T[1][1] + t0.T[1][3] * ty;
	t.T[1][2] = t0.T[1][2] + t0.T[1][3] * tz;

	t.T[2][0] = t0.T[2][0] + t0.T[2][3] * tx;
	t.T[2][1] = t0.T[2][1] + t0.T[2][3] * ty;
	t.T[2][2] = t0.T[2][2] + t0.T[2][3] * tz;

	t.T[3][0] = t0.T[3][0] + t0.T[3][3] * tx;
	t.T[3][1] = t0.T[3][1] + t0.T[3][3] * ty;
	t.T[3][2] = t0.T[3][2] + t0.T[3][3] * tz;
}


int main(int argc, char** argv) {
	ifstream in("african_head.obj");
	if (!in.is_open())
	{
		std::cout << "There is no african head\n";
		
	}
	TGAImage image(Width, Height, TGAImage::RGB);


	Matrix t;
	t.T[0][0] = 1; t.T[0][1] = 0; t.T[0][2] = 0; t.T[0][3] = 0;
	t.T[1][0] = 0; t.T[1][1] = 1; t.T[1][2] = 0; t.T[1][3] = 0;
	t.T[2][0] = 0; t.T[2][1] = 0; t.T[2][2] = 1; t.T[2][3] = 0;
	t.T[3][0] = 0; t.T[3][1] = 0; t.T[3][2] = 0; t.T[3][3] = 1;


	int* zbuffer = new int[Width * Height];
	for (int i = 0; i < Width * Height; i++) {
		zbuffer[i] = numeric_limits<int>::min();
	}


	Parallel(t, 1024, 1024, 0);

	string a;
	vector <O_Vertex> v, vn, vt;
	O_Vertex var;
	O_Vertex l(0, 0, 1);
	O_Vertex vision(0, 0, 1);
	vision.normalize();
	l.normalize();
	float z_max = numeric_limits<int>::min();
	float z_min = numeric_limits<int>::max();
	while (in >> a)
	{
		if (a == "v")
		{
			in >> var.x >> var.y >> var.z;
			v.push_back(var);
		}

		else
			if (a == "vt") //текстура
			{
				in >> var.x >> var.y >> var.z;
				vt.push_back(var);
			}
			else
				if (a == "vn") //нормаль
				{
					in >> var.x >> var.y >> var.z;
					vn.push_back(var);
				}
				else
					if (a == "f")
					{
						string b;
						string::size_type sz;
						O_Vertex world_coords[3], v_t[3], v_n[3], screen_coords[3];
						float v_[3], vt_[3], vn_[3];
						
						getline(in, b, '/');
						v_[0] = stof(b, &sz);
						getline(in, b, '/');
						vt_[0] = stof(b, &sz);
						getline(in, b, ' ');
						vn_[0] = stof(b, &sz);

						getline(in, b, '/');
						v_[1] = stof(b, &sz);
						getline(in, b, '/');
						vt_[1] = stof(b, &sz);
						getline(in, b, ' ');
						vn_[1] = stof(b, &sz);

						getline(in, b, '/');
						v_[2] = stof(b, &sz);
						getline(in, b, '/');
						vt_[2] = stof(b, &sz);
						getline(in, b);
						vn_[2] = stof(b, &sz);
						for (int i = 0; i < 3; i++)
						{
							screen_coords[i].x = 500 * v[(int)v_[i] - 1].x;
							screen_coords[i].y = 500 * v[(int)v_[i] - 1].y;
							screen_coords[i].z = 500 * v[(int)v_[i] - 1].z;
							Multiplication(screen_coords[i], t.T);
							world_coords[i] = v[(int)v_[i] - 1];
							v_t[i] = vt[(int)vt_[i] - 1];
							v_n[i] = vn[(int)vn_[i] - 1];
						}
						
						O_Vertex n = (world_coords[1] - world_coords[0]) ^ (world_coords[2] - world_coords[0]);
						
						n.normalize();
						float intensity = n * l;
						if (intensity > 0)
						{
							triangle_Z(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(255 * intensity, 255 * intensity, 255 * intensity, 255), zbuffer, z_max, z_min);
							v_n[0].normalize(); v_n[1].normalize(); v_n[2].normalize();
							
						}
						

					}
					else
						getline(in, a);

	}

	image.flip_vertically();
	image.write_tga_file("output.tga");
	

	for (int i = 0; i < Width; ++i)
		for (int j = 0; j < Height; ++j)
			zbuffer[i + j * Height] += -z_min;
	float z_total = -z_min + z_max;
	TGAImage z_buffer(Width, Height, TGAImage::RGB);
	for (int i = 0; i < Width; i++) 
	{
		for (int j = 0; j < Height; j++) 
		{
			
			if(zbuffer[i + j * 2048] > 0)
				z_buffer.set(i, j, TGAColor(255 * (zbuffer[i + j * 2048]) / z_total, 255 * (zbuffer[i + j * 2048]) / z_total, 
					255 * (zbuffer[i + j * 2048]) / z_total, 255));
		}
	}

	z_buffer.flip_vertically();
	z_buffer.write_tga_file("zbuffer.tga");

	
	return 0;
}
