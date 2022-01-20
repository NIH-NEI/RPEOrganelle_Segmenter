#ifndef rpeutil_h
#define rpeutil_h

const int FIX_BORDERS_DNA = 1;
const int FIX_BORDERS_ACTIN = 2;

std::vector<int> import_contours(int w, int h, const char *in_csv);
std::vector<int> import_reshape_contours(unsigned char *mask, int hm, int wm, const char *rs_csv);
std::vector<int> import_3d_contours(int w, int h, int d, const char *in_csv);

std::vector<int> cell_for_contour(int w, int h, int d, int z0,
		const std::vector<int>& flat_cells, const std::vector<int>& flat_cont);
std::vector<int> validate_contour(int w, int h,
		const std::vector<int>& flat_slices, const std::vector<int>& flat_cont);
std::vector<int> simplify_contour(const std::vector<int>& flat_cont);

std::vector<int> intersecting_slices(const std::vector<int>& flat_slices, const std::vector<int>& flat_cont);
std::vector<int> cut_slices(int w, int h,
		const std::vector<int>& flat_slices, const std::vector<int>& flat_cont);
std::vector<int> join_slices(int w, int h,
		const std::vector<int>& flat_slices, const std::vector<int>& flat_cont);

void export_3d_contours(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const std::vector<int> & flat_cells, const char *out_csv, bool validate, bool separate);
void export_2d_contours(unsigned char *mask, int hm, int wm,
		const std::vector<int>& flat_slices, const char *out_csv);

std::vector<int> cells_at_border(int w, int h, const std::vector<int> & flat_cells);
std::vector<int> contours_at_border(int w, int h, const std::vector<int> & flat_slices);

std::vector<int> mask_rcnn_to_particles(unsigned char *ptmask, int npts, int hptm, int wptm,
		int *rois, int nrois, int roisz, int x_orig=0, int y_orig=0);

std::vector<int> recombine_flat_particles(unsigned char *mask, int hm, int wm,
		std::vector<int> flat_particles, std::vector<double> scores, int fix_borders);

void mask2to3(unsigned char *mask, int hm, int wm);

void test_fix_borders_actin(unsigned char *mask, int hm, int wm);

#endif

