/*
 * sparse_ints.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#ifndef SPARSE_INTS_H_
#define SPARSE_INTS_H_

// Include some useful preprocessor macros
//   I realize that this sort of practice, in general, makes
//   code harder to read.  Thus, I've tried to keep the custom
//   macro usage to a minimum and documented the few that I do
//   use thoroughly in the my_macros.h file.
#include "my_macros.h"

#define WRITE_AS_FLOAT 0

// typedefs
#if WRITE_AS_FLOAT
typedef float value_t;
#else
typedef double value_t;
#endif

// defined constants
#define MASTER 0
#define NO_MORE_SHELLS -1
#define NEED_WORK 987
#define EMPTY 0
#define ZERO 1e-14
#define BUFF_SIZE 1048576



// mpqc includes
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/basis/distshpair.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/mbptr12/r12technology.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <util/misc/regtime.h>
#include <chemistry/qc/mbptr12/r12wfnworld.h>
#include <util/misc/consumableresources.h>
#include <util/group/pregtime.h>
#include <util/group/thread.h>
#include <math/scmat/repl.h>
#include <math/scmat/local.h>

// system library includes
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <arpa/inet.h>


typedef int16_t idx_t;
typedef std::pair<sc::RefSymmSCMatrix, sc::RefSymmSCMatrix> SymmSCMatrixPair;
typedef std::map<std::string, SymmSCMatrixPair > DensityMap;
typedef std::map<std::string, SymmSCMatrixPair >::iterator  DensityMapIterator;

typedef struct {
	sc::RefSymmSCMatrix first;
	sc::Ref<sc::GaussianBasisSet> first_basis;
	sc::RefSymmSCMatrix second;
	sc::Ref<sc::GaussianBasisSet> second_basis;
} PairDensityMatrices;


namespace sparse_ints {

typedef enum BinfileContent {
	AllInts = 0,
	MaxAbs = 2,
	Average = 4,
	StdDev = 8,
	Median = 16,
	AllIntsUntransformed = 32,
	DensityMatrix = 64
} binfile_type;

class MultiTimer {
	std::vector<sc::Ref<sc::RegionTimer> > rtimers_;
	std::vector<sc::Timer> timers_;
	sc::Ref<sc::ThreadGrp> thr_;
	sc::Ref<sc::MessageGrp> msg_;
	const char *name_;


public:

	MultiTimer(){ }

	MultiTimer(const char *name, sc::Ref<sc::MessageGrp> msg, sc::Ref<sc::ThreadGrp> thr){
		msg_ = msg;
		thr_ = thr;
		for_each(ithr, thr_->nthread()){
			name_ = name;
			sc::Ref<sc::RegionTimer> tim = new sc::ParallelRegionTimer(msg, name, 1, 1);
			rtimers_.push_back(tim);
			sc::Timer timerthr = sc::Timer(tim);
			timers_.push_back(timerthr);
		}
	}

	~MultiTimer() {}


	//void enter(const char *region, int threadnum = 0){
	//	timers_[threadnum].enter(region);
	//}
	void enter(const char *region, int threadnum = -1){
		if(threadnum == -1){
			for_each(ithr, thr_->nthread()){
				timers_[ithr].enter(region);
			}
		} else
			timers_[threadnum].enter(region);
	}

	//void exit(const char *region, int threadnum = 0){
	//	timers_[threadnum].exit(region);
	//}
	void exit(const char *region, int threadnum = -1){
		if(threadnum == -1){
			for_each(ithr, thr_->nthread()){
				timers_[ithr].exit(region);
			}
		} else
			timers_[threadnum].exit(region);
	}

	void print(std::ostream& o = sc::ExEnv::out0()){
		sc::Ref<sc::RegionTimer> main_timer = new sc::ParallelRegionTimer(msg_, name_, 1, 1);
		for(int ithr = 0; ithr < thr_->nthread(); ithr++){
			//rtimers_[ithr]->update_top();
			main_timer->merge(rtimers_[ithr]);
		}
		rtimers_[0]->print(o);
	}

};


typedef struct {
	bool debug;
	bool quiet;
	bool verbose;
	bool max_only;
	bool dynamic;
	bool do_tar;  // Should we make tar files instead of making one, big, bin file?
	int use_fake_ints;  // for debugging purposes
	binfile_type out_type;
} SparseIntOptions;

using namespace sc;

// Global variables
extern Ref<MessageGrp> msg;
extern Ref<ThreadGrp> thr;
extern MultiTimer timer;
extern SparseIntOptions opts;
extern Ref<ThreadLock> print_lock;


} // end namespace sparse_ints

#endif /* SPARSE_INTS_H_ */

