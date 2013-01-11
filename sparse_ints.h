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

// defined constants
#define MASTER 0
#define NO_MORE_SHELLS -1
#define NEED_WORK 987
#define EMPTY 0
#define ZERO 1e-14
#define BUFF_SIZE 1048576

// macros
#define DBG_MSG(mymsg) \
	if(opts.debug){\
		lock_->lock();\
		cout << "          " << msg->me() << "." << threadnum_ << mymsg << endl;\
		lock_->unlock();\
	}


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


// typedefs

typedef int16_t idx_t;
typedef float value_t;
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

typedef struct {
	bool debug = false;
	bool quiet = false;
	bool verbose = false;
	bool max_only = true;
	bool dynamic = true;
} SparseIntOptions;

using namespace sc;

// Global variables
extern Ref<MessageGrp> msg;
extern Ref<ThreadGrp> thr;
extern Timer timer;
extern SparseIntOptions opts;

class Globals {
public:
};


} // end namespace sparse_ints

#endif /* SPARSE_INTS_H_ */
