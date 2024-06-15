#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vec.h"
#include "tigr.h"

#define CHUNK_SIZE 1024

typedef struct face {
	size_t v[3], n[3], t[3];
} face;

typedef struct matrix {
	double m[4][4];
} matrix;

typedef struct trig {
	vec p[3];
	vec n[3];
	vec t[3];
} trig;

typedef struct mesh {
	size_t vert_len, vert_cap;
	vec *verts;
	size_t norm_len, norm_cap;
	vec *norms;
	size_t uv_len, uv_cap;
	vec *uvs;
	size_t face_len, face_cap;
	face *faces;
} mesh;

typedef struct camera {
	double near, far, fov, aspect_ratio, theta, phi;
	vec origin;
} camera;

typedef struct scene {
	mesh mesh;
	camera camera;
	vec light;
} scene;

typedef struct screen {
	size_t width, height;
	Tigr *scr;
	double *z;		
} screen;

void add_vert(mesh *m, vec *v) {
	if (m->vert_cap == m->vert_len)
		m->verts = realloc(m->verts, (m->vert_cap += CHUNK_SIZE) * sizeof(vec));
	memcpy(&m->verts[m->vert_len++], v, sizeof(vec));
}

void add_face(mesh *m, face *f) {
	if (m->face_cap == m->face_len)
		m->faces = realloc(m->faces, (m->face_cap += CHUNK_SIZE) * sizeof(face));
	memcpy(&m->faces[m->face_len++], f, sizeof(face));
}

void add_norm(mesh *m, vec *n) {
	if (m->norm_cap == m->norm_len)
		m->norms = realloc(m->norms, (m->norm_cap += CHUNK_SIZE) * sizeof(vec));
	memcpy(&m->norms[m->norm_len++], n, sizeof(vec));
}

void add_uv(mesh *m, vec *uv) {	
	if (m->uv_cap == m->uv_len)
		m->uvs = realloc(m->uvs, (m->uv_cap += CHUNK_SIZE) * sizeof(vec));
	memcpy(&m->uvs[m->uv_len++], uv, sizeof(vec));
}

void matmul(vec i, matrix m, vec *o) {
	o->x = i.x * m.m[0][0] +  i.y * m.m[1][0] + i.z * m.m[2][0] + m.m[3][0];
	o->y = i.x * m.m[0][1] +  i.y * m.m[1][1] + i.z * m.m[2][1] + m.m[3][1];
	o->z = i.x * m.m[0][2] +  i.y * m.m[1][2] + i.z * m.m[2][2] + m.m[3][2];

	double w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + m.m[3][3];

	vec_scale(o, 1 / w);
}


void draw_trig(screen *scr, trig t, vec light) {
	if (t.p[0].x > 1 && t.p[1].x > 1 && t.p[2].x > 1) return;
	if (t.p[0].x < -1 && t.p[1].x < -1 && t.p[2].x < -1) return;
	if (t.p[0].y > 1 && t.p[1].y > 1 && t.p[2].y > 1) return;
	if (t.p[0].y < -1 && t.p[1].y < -1 && t.p[2].y < -1) return;

	for (int i = 0; i < 3; i++) {
		t.p[i].x = (t.p[i].x + 1) * scr->width / 2;
		t.p[i].y = (t.p[i].y + 1) * scr->height / 2;
		// printf("%lf %lf\n", t.p[i].x, t.p[i].y);
	}

	vec v1, v2, v3, v4;
	vec n1, n2, n3, n4;
	if (t.p[0].y < t.p[1].y && t.p[0].y < t.p[2].y) {
		v1 = t.p[0];
		n1 = t.n[0];
		v2 = t.p[1];
		n2 = t.n[1];
		v3 = t.p[2];
		n3 = t.n[2];
	} else if (t.p[1].y < t.p[2].y && t.p[1].y < t.p[0].y) {
		v1 = t.p[1];
		n1 = t.n[1];
		v2 = t.p[2];
		n2 = t.n[2];
		v3 = t.p[0];
		n3 = t.n[0];
	} else {
		v1 = t.p[2];
		n1 = t.n[2];
		v2 = t.p[0];
		n2 = t.n[0];
		v3 = t.p[1];
		n3 = t.n[1];
	}
	if (v3.y < v2.y) {
		v4 = v3;
		n4 = n3;
		v3 = v2;
		n3 = n2;		
		v2 = v4;
		n2 = n4;
	}	
	v4.y = v2.y;
	v4.x = ((v1.y - v2.y) * v3.x + (v2.y - v3.y) * v1.x) / (v1.y - v3.y);
	v4.z = ((v1.y - v2.y) * v3.z + (v2.y - v3.y) * v1.z) / (v1.y - v3.y);
	
	if (v2.x > v4.x) {
		double tmp = v2.x;
		v2.x = v4.x;
		v4.x = tmp;
	}

	{
		vec tmp = n3;
		n4 = n1;
		vec_scale(&tmp, v1.y - v2.y);
		vec_scale(&n4, v2.y - v3.y);
		vec_add(&n4, tmp);
		vec_scale(&n4, 1 / (v1.y - v3.y));
	}

	v1.x = floor(v1.x);
	v1.y = floor(v1.y);
	v2.x = floor(v2.x);
	v2.y = floor(v2.y);
	v3.x = floor(v3.x);
	v3.y = ceil(v3.y);
	v4.x = ceil(v4.x);
	v4.y = floor(v4.y);

	{
		double xl = v1.x,
			   xr = v1.x,
			   dxl = (v2.x - v1.x) / (v2.y - v1.y),
			   dxr = (v4.x - v1.x) / (v2.y - v1.y),
			   zl = v1.z,
			   zr = v1.z,
			   dzl = (v2.z - v1.z) / (v2.y - v1.y),
			   dzr = (v4.z - v1.z) / (v2.y - v1.y);
		vec nl = n1, nr = n1, dnl = n2, dnr = n4;
		vec_sub(&dnl, n1);
		vec_scale(&dnl, 1 / (v2.y - v1.y));
		vec_sub(&dnr, n1);
		vec_scale(&dnr, 1 / (v2.y - v1.y));
		for (double y = v1.y; y <= v2.y; y++) {
			double z = zl, dz = (zr - zl) / (xr - xl);
			vec n = nl, dn = nr;
			vec_sub(&dn, nl);
			vec_scale(&dn, 1 / (xr - xl));
			for (int x = xl; x <= xr; x++, z+=dz, vec_add(&n, dn)) {
				if (y < 0 || y >= scr->height) continue;
				if (x < 0 || x >= scr->width) continue;
				if (z < scr->z[(int) y * scr->width + x]) {
					double dot = vec_dot(n, light) / vec_dot(n, n);
					double brightness = (1 - dot) / 2;
					TPixel c = tigrRGB(255 * brightness, 255 * brightness, 255 * brightness);

					tigrPlot(scr->scr, x, y, c);
					scr->z[(int) y * scr->width + x] = z;
				}
			}
			xl += dxl;
			xr += dxr;
			zl += dzl;
			zr += dzr;
			vec_add(&nl, dnl);
			vec_add(&nr, dnr);
		}
	}
	{
		double xl = v3.x,
			   xr = v3.x,
			   dxl = (v2.x - v3.x) / (v3.y - v2.y),
			   dxr = (v4.x - v3.x) / (v3.y - v2.y),
			   zl = v3.z,
			   zr = v3.z,
			   dzl = (v2.z - v3.z) / (v3.y - v2.y),
			   dzr = (v4.z - v3.z) / (v3.y - v2.y);
		vec nl = n1, nr = n1, dnl = n2, dnr = n4;
		vec_sub(&dnl, n1);
		vec_scale(&dnl, 1 / (v2.y - v1.y));
		vec_sub(&dnr, n1);
		vec_scale(&dnr, 1 / (v2.y - v1.y));
		for (double y = v3.y; y >= v2.y; y--) {
			double z = zl, dz = (zr - zl) / (xr - xl);
			vec n = nl, dn = nr;
			vec_sub(&dn, nl);
			vec_scale(&dn, 1 / (xr - xl));
			for (int x = xl; x <= xr; x++, z+=dz, vec_add(&n, dn)) {
				if (y < 0 || y >= scr->height) continue;
				if (x < 0 || x >= scr->width) continue;
				if (z < scr->z[(int) y * scr->width + x]) {
					double dot = vec_dot(n, light) / vec_dot(n, n);
					double brightness = (1 - dot) / 2;
					TPixel c = tigrRGB(255 * brightness, 255 * brightness, 255 * brightness);

					tigrPlot(scr->scr, x, y, c);
					scr->z[(int) y * scr->width + x] = z;
				}
			}
			xl += dxl;
			xr += dxr;
			zl += dzl;
			zr += dzr;
		}
	}

	if (v2.y < 0 || v2.y >= scr->height) return;

	{
		double z = v2.z,
			   dz = (v4.z - v2.z) / (v4.x - v2.x);
		vec n = n2, dn = n4;
		vec_sub(&dn, n);
		vec_scale(&dn, 1 / (v4.x - v2.x));
		for (int x = v2.x; x < v4.x; x++, z+=dz, vec_add(&n, dn)) {
			if (x < 0 || x >= scr->width) continue;
			if (z < scr->z[(int) v2.y * scr->width + x]) {
				double dot = vec_dot(n, light) / vec_dot(n, n);
				double brightness = (1 - dot) / 2;
				TPixel c = tigrRGB(255 * brightness, 255 * brightness, 255 * brightness);
				
				tigrPlot(scr->scr, x, v2.y, c);
				scr->z[(int) v2.y * scr->width + x] = z;
			}

		}
	}
	
	/*
	tigrLine(scr->scr, t.p[0].x, t.p[0].y, t.p[1].x, t.p[1].y, c);
	tigrLine(scr->scr, t.p[1].x, t.p[1].y, t.p[2].x, t.p[2].y, c);
	tigrLine(scr->scr, t.p[2].x, t.p[2].y, t.p[0].x, t.p[0].y, c);
	*/
}

void draw(screen *scr, scene *scn) {
	// Resets the z buffer	
	for (int i = 0; i < scr->width * scr->height; i++) {
		scr->z[i] = INFINITY;
	}
	
	double hfov = tan(scn->camera.fov * M_PI / 360),
			vfov = hfov / scn->camera.aspect_ratio,
			proj_A = (scn->camera.near + scn->camera.far) / (scn->camera.far - scn->camera.near),
			proj_B = 2 * scn->camera.near * scn->camera.far / (scn->camera.far - scn->camera.near);
	matrix proj = {{
		{ 1 / hfov,        0,      0,      0 },
		{        0, 1 / vfov,      0,      0 },
		{        0,        0, proj_A,      1 },
		{        0,        0, proj_B,      0 }
	}};

	double ct = cos(scn->camera.theta),
		   st = sin(scn->camera.theta);
	matrix rot_y = {{	
		{ ct, 0, -st, 0 },
		{  0, 1,   0, 0 },
		{ st, 0,  ct, 0 },
		{  0, 0,   0, 1 } 
	}};

	double cp = cos(scn->camera.phi),
		   sp = sin(scn->camera.phi);
	matrix rot_z = {{
		{ 1,   0,  0, 0 },
		{ 0,  cp, sp, 0 },
		{ 0, -sp, cp, 0 },
		{ 0,   0,  0, 1 }
	}};

	for (int i = 0; i < scn->mesh.face_len; i++) {
		face f = scn->mesh.faces[i];

		vec v1, v2, v3, n1, n2, n3, tmp;		
		matmul(scn->mesh.verts[f.v[0] - 1], rot_y, &tmp);
		matmul(tmp, rot_z, &v1);
		matmul(scn->mesh.verts[f.v[1] - 1], rot_y, &tmp);
		matmul(tmp, rot_z, &v2);
		matmul(scn->mesh.verts[f.v[2] - 1], rot_y, &tmp);
		matmul(tmp, rot_z, &v3);
		
		vec_sub(&v1, scn->camera.origin);
		vec_sub(&v2, scn->camera.origin);
		vec_sub(&v3, scn->camera.origin);

		trig t;	
		matmul(v1, proj, &t.p[0]);
		matmul(v2, proj, &t.p[1]);
		matmul(v3, proj, &t.p[2]);
		
		matmul(scn->mesh.norms[f.n[0] - 1], rot_y, &tmp);
		matmul(tmp, rot_z, &t.n[0]);
		matmul(scn->mesh.norms[f.n[1] - 1], rot_y, &tmp);
		matmul(tmp, rot_z, &t.n[1]);
		matmul(scn->mesh.norms[f.n[2] - 1], rot_y, &tmp);
		matmul(tmp, rot_z, &t.n[2]);

		if (fabs(t.p[0].z) > 1 || fabs(t.p[1].z) > 1 || fabs(t.p[2].z) > 1) continue;
		draw_trig(scr, t, scn->light);	

	}
}


void load_obj(FILE *f, mesh *m) {
	char s[16] = {0};
	while(fscanf(f, "%15s", s) != EOF) {
		if (strcmp(s, "v") == 0) {
			vec v;
			fscanf(f, " %lf %lf %lf", &v.x, &v.y, &v.z);
			add_vert(m, &v);
		} else if (strcmp(s, "vn") == 0) {
			vec v;
			fscanf(f, " %lf %lf %lf", &v.x, &v.y, &v.z);
			add_norm(m, &v);
		} else if (strcmp(s, "vt") == 0) {
			vec v;
			fscanf(f, " %lf %lf %lf", &v.x, &v.y, &v.z);
			add_uv(m, &v);
		} else if (strcmp(s, "f") == 0) {
			face b;
			for (int i = 0; i < 3; i++) {
				fscanf(f, "%zu", &b.v[i]);
				if (fgetc(f) != '/') continue;
				fscanf(f, "%zu", &b.t[i]);
				if (fgetc(f) != '/') continue;
				fscanf(f, "%zu", &b.n[i]);
			}
			add_face(m, &b);
		}
	}
};

#define WIDTH 1920
#define HEIGHT 1080

int main(int argc, char *argv[]) {	
	if (argc != 2) {
		printf("usage: %s file.obj\n", argv[0]);
		return 1;
	}

	FILE *f = fopen(argv[1], "r");
	if (!f) {
		printf("Failed to load model\n");
		return 1;
	}
	mesh m = {0};
	load_obj(f, &m);
	fclose(f);

	double zmax = 0;
	for (int i = 0; i < m.vert_len; i++) {
		zmax = m.verts[i].z > zmax? m.verts[i].z: zmax;	
	}

	vec light = { 0.56, -0.56, -0.56 };
	camera c = {0.1, 10000.0, 90.0, (float) WIDTH / HEIGHT, 0, 0, {0, 0, 3 * zmax}}; 
	Tigr *screen = tigrWindow(WIDTH, HEIGHT, "c3d", 0);
	tigrUpdate(screen);	

	double *z = malloc((WIDTH * HEIGHT + 2) * sizeof(double));
	
	int t = 0;
	while (!tigrClosed(screen)) {
		if (tigrKeyHeld(screen, TK_LEFT)) c.theta += 1./60;
		if (tigrKeyHeld(screen, TK_RIGHT)) c.theta -= 1./60;
		if (tigrKeyHeld(screen, TK_UP)) c.phi += 1./60;
		if (tigrKeyHeld(screen, TK_DOWN)) c.phi -= 1./60;	
		if (tigrKeyHeld(screen, 'I')) c.origin.z *= 60./61;	
		if (tigrKeyHeld(screen, 'O')) c.origin.z *= 61./60;	
		if (tigrKeyHeld(screen, 'W')) c.origin.y += c.origin.z / 100;	
		if (tigrKeyHeld(screen, 'S')) c.origin.y -= c.origin.z / 100;	
		if (tigrKeyHeld(screen, 'A')) c.origin.x += c.origin.z / 100;	
		if (tigrKeyHeld(screen, 'D')) c.origin.x -= c.origin.z / 100;	
		if (tigrKeyHeld(screen, 'Z')) c.fov *= 61./60;	
		if (tigrKeyHeld(screen, 'X')) c.fov *= 60./61;	

		scene scn = { m, c, light };
		struct screen scr = { WIDTH, HEIGHT, screen, z };

		tigrClear(screen, tigrRGB(0, 0, 0));
		draw(&scr, &scn);
		tigrUpdate(screen);	
		printf("%d\n", t++);
	}
	
	free(z);
	return 0;
}
