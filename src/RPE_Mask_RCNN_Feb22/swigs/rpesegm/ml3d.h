#ifndef ml3d_h
#define ml3d_h

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif

#include <algorithm>
#include "raster.h"

struct LinkML
{
	int z, idx;
	LinkML() : z(-1), idx(-1) {}
	LinkML(int _z, int _idx) : z(_z), idx(_idx) {}
};

class SliceML : public Slice
{
public:
	int z, idx;
	int ptid;
	double sc;
	SliceML() : Slice(), z(-1), idx(-1), ptid(-1), sc(0.) {}
	bool valid() { return area>0; }
	bool taken() { return ptid>=0; }
};

struct SortScoreML
{
	int ptid;
	double sc;
	SortScoreML() : ptid(-1), sc(0.) {}
	SortScoreML(int _ptid, double _sc) : ptid(_ptid), sc(_sc) {}
};

class AssemblerML
{
public:
	double MIN_DUP_OVL = 0.5;
	double GOOD_IOU = 0.75;
	double ACCEPTABLE_IOU = 0.65;
	int MAX_GAP = 3;
	int MIN_SEED = 3;
	double DNA_PHANTOM_FILL = 0.6;

	int w, h, d;
	Raster3D mstack;
	std::vector<std::vector<SliceML>> all_slices;
	std::vector<Particle3D> cells;
	std::vector<SortScoreML> sorted_cells;
	
	AssemblerML(unsigned char *mask3d, int zm3d, int hm3d, int wm3d) :
		w(wm3d), h(hm3d), d(zm3d),
		mstack(wm3d, hm3d, zm3d, mask3d),
		all_slices(zm3d)
	{
		
	}

	void load_flat_slices(std::vector<std::vector<int>> & all_flat_particles);
	void de_dupe();
	void initial_linkage();
	void fill_gaps();
	void smooth_cells();
	void fix_borders_actin();
	void fix_borders_dna();
	int filter_cells(int min_height=4, double min_sc=5.);
	
protected:
	SliceML& lnsl(LinkML& ln) { return all_slices[ln.z][ln.idx]; }
	std::vector<int> find_dups(SliceML& ptc);
	int find_below(SliceML& ptc, int z, double good_iou=-1.);
	std::vector<int> find_above(SliceML& ptc, int z);
	double chain_score(std::vector<LinkML>& chain);
	double fwd_chain_score(SliceML& ptc, std::vector<LinkML>& chain, int z1);
	double nearby_chain_score(SliceML& ptc);
	void interpolate_frame(Particle3D& cell, int zmiss);
	
};

#endif
