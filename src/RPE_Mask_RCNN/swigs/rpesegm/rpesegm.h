#ifndef rpesegm_h
#define rpesegm_h

const int POSTPROC_NONE = 0;
const int POSTPROC_DNA = 1;
const int POSTPROC_ACTIN = 2;

void assemble_basic(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile, int postproc=0);

void assemble_ml(std::vector<std::vector<int>> all_flat_particles,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d, const char *csvfile, int postproc=POSTPROC_NONE);

void segment_dna(int w, int h,
	unsigned short *data,
	unsigned char *angle,
	unsigned char *mask,
	int otsu1, int otsu2);
void assemble_dna(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile, bool validate=true);

void segment_actin(int w, int h,
	unsigned short *data,
	unsigned char *angle,
	unsigned char *mask,
	int otsu1, int otsu2);
void segment_actin_z01(int w, int h,
	unsigned short *data,
	unsigned char *angle,
	unsigned char *mask,
	unsigned char *mask2,
	int otsu1, int otsu2);
void assemble_actin(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
	const char *csvfile, const char *dna_csvfile, const char *z01_csvfile, bool validate=true);
void detect_top_bottom(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
	const char *actin_csvfile, const char *dna_csvfile);

void extract_z01(unsigned short *data3d, int zd3d, int hd3d, int wd3d);

// Full-frame (old) Z01 extraction algorithm
void combo_z01(unsigned short *data3d, int zd3d, int hd3d, int wd3d);

void normalize_frame(unsigned short *sdata, int shd, int swd, int kernel_size);
int adjusted_threshold(unsigned short *data, int hd, int wd, int thresh, double pct);
void export_2d_segmentation(unsigned char *mask, int hm, int wm, const char *csvfile);
int import_2d_segmentation(unsigned char *mask, int hm, int wm, const char *csvfile);

int import_3d_segmentation(unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile, bool border=true);

void compare_with_reshape(int w, int h,
	unsigned char *mask,
	unsigned char *rs_mask,
	const char *rs_csvfile,
	const char *out_csvfile);
void colorize_reshape_comparison(int w, int h,
	unsigned char *mask3d,
	unsigned short *data,
	unsigned char *mask2);
void read_reshape(unsigned char *mask, int hm, int wm, const char *csvfile);

void match_z01(unsigned char *mask, int hm, int wm, unsigned char *mask2, int hm2, int wm2);

int export_mask_id_3d(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		const char *csvfile, int num_dilations=0);

int import_mask_id_3d(unsigned short *data3d, int zd3d, int hd3d, int wd3d,
		unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
		const char *csvfile);

#endif
