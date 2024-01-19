#ifndef flatcont_h
#define flatcont_h

#include "raster.h"

// "at the border" qualifiers
#define BORDER_SIZE 3
#define AT_BORDER_PCT 15.
#define AT_BORDER_IMPORT_PCT 20.

Boundary boundary_from_path(std::vector<Point> & path);
double out_of_boundary(Boundary & bnd, const std::vector<Point> & path, double *p_len);
void simplify_path(std::vector<Point> &path);

bool detect_particle_contour(Raster8 & msk, Contour& cont, Particle &ptc, unsigned char bordc, bool simplify=true);

std::vector<int> particles_to_contours(Raster8 &msk, std::vector<Particle> & particles);

struct CellContour
{
	int cell_id;
	int z;
	int cont_idx;
	Contour cont;
	bool valid;
	Slice ptc;
};

std::vector<CellContour> unflatten_cell_contours(const std::vector<int>& flat_cells);


#endif
