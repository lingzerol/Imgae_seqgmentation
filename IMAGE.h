#pragma once
#ifndef IMAGE_DEBUG
#define IMAGE_DEBUG
#include <core.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <queue>
#define RGB_PIXEL false
#define LAB_PIXEL true

const int direction[8][2] = { {1,0}, {0,1},{-1,0},{0,-1},{1,1},{1,-1},{-1,1},{-1,-1} };

class pixel {
public:
	friend double dc2(const pixel&, const pixel&);
	friend double ds2(const pixel&, const pixel&);
	friend class IMAGE;
	friend class compare;
	pixel() :_p{ 0,0,0 },_k(-1),_r(-1), _c(-1),_d(DBL_MAX) {}
	pixel(double l, double a, double b,long long r,long long c) :_p{ l,a,b },_k(-1),_r(r),_c(c),_d(DBL_MAX) {}
	void upgrade(const pixel&r) {
		_p[0] += r._p[0];
		_p[1] += r._p[1];
		_p[2] += r._p[2];
		_r += r._r;
		_c += r._c;
	}
	void divide(int num) {
		_p[0] /= num;
		_p[1] /= num;
		_p[2] /= num;
		_r /= num;
		_c /= num;
	}
private:
	double _p[3];
	long long _k;
	long long _r,_c;
	double _d;
};

class compare {
public:
	bool operator ()(const pixel&a, const pixel&b) {
		return a._d > b._d;
	}
};

class IMAGE {
public:
	IMAGE(const cv::Mat&source);
	IMAGE(const int const*, const int const*, const int const*,int,int);
	int width() { return _width; }
	int height() { return _height; }
	pixel operator[](int i) { return _img[i]; }
	void rgb_to_lab();
	double D2(const pixel&a, const pixel&b,double m=1) {
		return dc2(a,b) + ds2(a,b) * (m*m / (_S*_S)); }
	std::vector<int> SLIC(int *&,int &, int, int, int, double, double);
	int index(int r, int c) { return r * _width + c; }
private:
	void find_seed(int&,std::vector<pixel>&);
	std::vector<pixel>_img;
	const int _width;
	const int _height;
	double _S;
};

IMAGE::IMAGE(const cv::Mat &source):_width(source.cols),_height(source.rows),_S(source.cols*source.rows) {
	for (int i = 0; i < source.rows; ++i) {
		for (int j = 0; j < source.cols; ++j) {
			_img.push_back(pixel(source.ptr(i)[j*3], source.ptr(i)[j*3 + 1], source.ptr(i)[j*3 + 2], i, j));
		}
	}
}
IMAGE::IMAGE(const int const*p0, const int const*p1, const int const*p2,int height,int width):_width(width),_height(height),_S(width*height) {
	int sz = height * width;
	for (int i = 0; i < sz; ++i) {
		_img.push_back(pixel(p0[i], p1[i], p2[i], i / width, i%width));
	}
}
double dc2(const pixel &a, const pixel&b) {
	return (a._p[0] - b._p[0])*(a._p[0] - b._p[0]) +
		(a._p[1] - b._p[1])*(a._p[1] - b._p[1]) +
		(a._p[2] - b._p[2])*(a._p[2] - b._p[2]);
}
double ds2(const pixel &a, const pixel &b) {
	return (a._r - b._r)*(a._r - b._r) + (a._c - b._c)*(a._c - b._c);
}

void IMAGE::rgb_to_lab() {
	int num = _img.size();
	for (int i = 0; i < num; ++i) {
		double R = _img[i]._p[0] / 255.0;
		double G = _img[i]._p[1] / 255.0;
		double B = _img[i]._p[2] / 255.0;

		if (R > 0.04045) R = pow((R + 0.055) / 1.055, 2.4);
		else R = R / 12.92;
		if (G > 0.04045) G = pow((G + 0.055) / 1.055, 2.4);
		else G = G / 12.92;
		if (B > 0.04045) B = pow((B + 0.055) / 1.055, 2.4);
		else B = B / 12.92;

		double X = R * 0.4124 + G * 0.3576 + B * 0.1805;
		double Y = R * 0.2126 + G * 0.7152 + B * 0.0722;
		double Z = R * 0.0193 + G * 0.1192 + B * 0.9505;

		X = X / 95.047;
		Y = Y / 100.0;
		Z = Z / 108.883;
		
		double bounce = pow(6.0 / 29.0, 3);

		if (X > bounce)X = pow(X, 1.0/3.0);
		else X = 1.0 / 3.0*(29.0 / 6.0)*(29.0 / 6.0)*X + 4.0 / 29.0;
		if (Y > bounce)Y = pow(Y, 1.0 / 3.0);
		else Y = 1.0 / 3.0*(29.0 / 6.0)*(29.0 / 6.0)*Y + 4.0 / 29.0;
		if (Z > bounce)Z = pow(Z, 1.0 / 3.0);
		else Z = 1.0 / 3.0*(29.0 / 6.0)*(29.0 / 6.0)*Z + 4.0 / 29.0;
		/*
		const double epsilon = 0.008856;	//actual CIE standard
		const double kappa = 903.3;		//actual CIE standard

		const double Xr = 0.950456;	//reference white
		const double Yr = 1.0;		//reference white
		const double Zr = 1.088754;	//reference white
		X = X / Xr;
		Y = Y / Yr;
		Z = Z / Zr;

		if (X > epsilon)	X = pow(X, 1.0 / 3.0);
		else				X = (kappa*X + 16.0) / 116.0;
		if (Y > epsilon)	Y = pow(Y, 1.0 / 3.0);
		else				Y = (kappa*Y + 16.0) / 116.0;
		if (Z > epsilon)	Z = pow(Z, 1.0 / 3.0);
		else				Z = (kappa*Z + 16.0) / 116.0;
		*/
		

		_img[i]._p[0] = 116 * Y - 16;
		_img[i]._p[1] = 500 * (X - Y);
		_img[i]._p[2] = 200 * (Y - Z);
		//std::cout << _img[i]._p[0] << " " << _img[i]._p[1] << " " << _img[i]._p[2] << " " << std::endl;
	}
}
//label is to log the pixel's label, 
//lenght is to tell how many pixel there is, 
//numk is to tell how many superpixels there are, 
//compactness is to tell whether the cluster is compactness, 
//least is to tell the mini size of the cluster, 
//the combine is to tell the threshold of the simiarity between pixel, 
//range is to tell the farthest pixel of the cluster can touch.
std::vector<int> IMAGE::SLIC(int *&label,int &length, int numk=100, int compactness=10, int least=100, double combine=2, double scope=1) {
	std::priority_queue<pixel,std::vector<pixel>,::compare> q;
	std::vector<pixel> cluster_center;
	std::vector<int> number;
	std::vector<int> seq;
	find_seed(numk, cluster_center);
	label = new int[_width*_height];
	length = _width * _height;
	number.resize(cluster_center.size());
	memset(label, -1, sizeof(int)*_width*_height);
	int pixel_count = 0;
	for (int i = 0; i < cluster_center.size(); ++i) {
		pixel t = cluster_center[i];
		cluster_center[i]._p[0] = cluster_center[i]._p[0]/combine;
		cluster_center[i]._p[1] = cluster_center[i]._p[1]/combine;
		cluster_center[i]._p[2] = cluster_center[i]._p[2]/combine;
		cluster_center[i]._r = cluster_center[i]._r / scope;
		cluster_center[i]._c = cluster_center[i]._c / scope;
		q.push(t);
		number[i] = 0;
	}
	//std::cout << cluster_center.size() << std::endl;
	const int CONNECTIVITY = 4;
	while (!q.empty()) {
		pixel temp = q.top();
		q.pop();
		int k = temp._k; 
		//cluster_center[k].upgrade(temp);
		if (label[index(temp._r,temp._c)] < 0) {
			label[index(temp._r, temp._c)] = k;
			seq.push_back(index(temp._r, temp._c));
			++pixel_count;
			cluster_center[k].upgrade(temp);
			++number[k];
			for (int i = 0; i < CONNECTIVITY; ++i) {
				int r, c;
				r =	temp._r+direction[i][0];
				c = temp._c+direction[i][1];
				if (!(r < 0 || r >= _height || c < 0 || c >= _width)) {
					int p = index(r,c);
					if (label[p] < 0) {
						pixel t = _img[p];
						pixel tt = cluster_center[k];
						tt.divide(number[k]);
						t._d = D2(t, tt, compactness);
						t._k = k;
						q.push(t);
					}
				}
			}
		}
	}

		for (int i = 0; i < _height; ++i) {
			for (int j = 0; j < _width; ++j) {
				int p = index(i, j);
				if (label[p] < 0) {
					if (p > 0 && label[p - 1]>=0)label[p] = label[p - 1];
					else if (p < _width*_height - 1 && label[p + 1] >= 0)label[p] = label[p + 1];
				}
				else if (number[label[p]] <= least) {
					--number[label[p]];
					if (p > 0 && label[p - 1] != label[p]) { 
						label[p] = label[p - 1]; 
						++number[label[p-1]];
					}
					else if (p < _width*_height - 1 && label[p + 1] != label[p]) { 
						label[p] = label[p + 1]; 
						++number[label[p+1]];
					}
				}
			}
		}
		return seq;
}
void IMAGE::find_seed(int& numk, std::vector<pixel> &q) {
	const int sz = _width * _height;
	int gridstep = sqrt(double(sz) / double(numk)) + 0.5;
	int halfstep = gridstep / 2;
	double h = _height; double w = _width;
	
	int xsteps = int(_width / gridstep);
	int ysteps = int(_height / gridstep);
	int err1 = abs(xsteps*ysteps - numk);
	int err2 = abs(int(_width / (gridstep - 1))*int(_height / (gridstep - 1)) - numk);
	if (err2 < err1)
	{
		gridstep -= 1.0;
		xsteps = _width / (gridstep);
		ysteps = _height / (gridstep);
	}
	//_S = gridstep*gridstep;
	int k = 0;
	for (int i = halfstep; i < _height; i += gridstep) {
		for (int j = halfstep; j < _width; j += gridstep) {
			pixel temp = _img[index(i, j)];
			temp._d = 0;
			temp._k = k++;
			q.push_back(temp);
		}
	}
	numk = q.size();
	//_S = (double)_width * _height / numk;
}
#endif