#include <math.h>

typedef struct vec {
	double x, y, z;	
} vec;

void vec_add(vec *v1, vec v2) {
	v1->x += v2.x;
	v1->y += v2.y;
	v1->z += v2.z;
}

void vec_sub(vec *v1, vec v2) {
	v1->x -= v2.x;
	v1->y -= v2.y;
	v1->z -= v2.z;
}

void vec_scale(vec *v, double s) {
	v->x *= s;
	v->y *= s;
	v->z *= s;
}

double vec_dot(vec v1, vec v2) {
	return (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
}

void vec_unit(vec *v) {
	vec_scale(v, vec_dot(*v ,*v));
}

