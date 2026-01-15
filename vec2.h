#pragma once

#include <cmath>

struct vec2 {
	double x = 0., y = 0.;

	vec2() = default;
	
	vec2(const vec2& vec)
		: x(vec.x), y(vec.y) {
	};

	vec2(double _x, double _y)
		: x(_x), y(_y) {
	};

	double operator*(const vec2& vec) const {
		return x * vec.x + y * vec.y;
	}

	vec2 operator*(double k) const {
		return { x * k, y * k };
	}

	vec2 operator/(double k) const {
		return { x / k, y / k };
	}

	vec2 operator+(const vec2& vec) const {
		return { x + vec.x, y + vec.y };
	}

	vec2 operator-() const {
		return { -x, -y };
	}

	vec2 operator-(const vec2& vec) const {
		return { x - vec.x, y - vec.y };
	}

	vec2 operator+=(const vec2& vec) {
		x += vec.x;
		y += vec.y;
		return *this;
	}

	vec2 operator-=(const vec2& vec) {
		x -= vec.x;
		y -= vec.y;
		return *this;
	}

	vec2 operator*=(double k) {
		x *= k;
		y *= k;
		return *this;
	}

	vec2 operator/=(double k) {
		x /= k;
		y /= k;
		return *this;
	}

	//Компонента z векторного произведения
	double crossZ(const vec2& vec) {
		return x * vec.y - y * vec.x;
	}

	double norm() const {
		return sqrt(x * x + y * y);
	}

	vec2 normalize() const {
		return *this * (1. / norm());
	}
};

//template<typename fp>
//vec2<fp> operator*(fp k, const vec2<fp>& v) {
//	return v * k;
//}

vec2 operator*(double k, const vec2& v);