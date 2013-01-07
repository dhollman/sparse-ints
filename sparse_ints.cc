//
// ao_trans.cc
//
// Copyright (C) 2013 David Hollman
//
// Author: David S. Hollman <dhollman@uga.edu>
// Maintainer: DSH
//
// License
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
////////////////////////////////////////////////////////////////////////////////////

// Include some useful preprocessor macros
//   I realize that this sort of practice, in general, makes
//   code harder to read.  Thus, I've tried to keep the custom
//   macro usage to a minimum and documented the few that I do
//   use thoroughly in the my_macros.h file.
#include "my_macros.h"

#define MASTER 0
#define NO_MORE_SHELLS -1
#define NEED_WORK 987
#define EMPTY 0
#define ZERO 1e-14

#define DBG_MSG(mymsg) \
	if(debug){\
		lock_->lock();\
		cout << "          " << msg->me() << "." << threadnum_ << mymsg << endl;\
		lock_->unlock();\
	}

#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
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
#include <arpa/inet.h>

using namespace sc;
using namespace std;

typedef int16_t idx_t;
typedef float value_t;
typedef std::pair<RefSymmSCMatrix, RefSymmSCMatrix> SymmSCMatrixPair;
typedef std::map<string, SymmSCMatrixPair >  DensityMap;
typedef std::map<string, SymmSCMatrixPair >::iterator  DensityMapIterator;

Timer timer;
bool low_comm_mode = false;
bool dynamic;
int quartets_per_task;
int max_cutoff;
bool print_distribution;
int num_masters = 1;
int n;
bool debug = false;
bool quiet = false;
bool verbose = false;
bool max_only = true;
bool max_of_all = true;
const double cutoff = pow(10.0, -double(max_cutoff));
double* cutoffs;
Ref<ThreadGrp> thr;
Ref<MessageGrp> msg;
int buff_size = 1024*1024;

template <class T=value_t>
inline T
max_abs(RefSCMatrix& mat){
	int nrows = mat->nrow();
	int ncols = mat->ncol();
	T maxval = 0.0;
	for_each(row,nrows, col,ncols){
		T tmpval = (T)fabs(mat->get_element(row, col));
		if(tmpval > maxval){
			maxval = tmpval;
		}
	}
	return maxval;
}

typedef struct {
	RefSymmSCMatrix first;
	Ref<GaussianBasisSet> first_basis;
	RefSymmSCMatrix second;
	Ref<GaussianBasisSet> second_basis;
} PairDensityMatrices;

class ComputeThread : public Thread {


    int threadnum_;
    Ref<ThreadLock> lock_;
    Ref<GaussianBasisSet> basis1_;
    Ref<GaussianBasisSet> basis2_;
    Ref<GaussianBasisSet> basis3_;
    Ref<GaussianBasisSet> basis4_;
    Ref<TwoBodyInt> inteval_;
    TwoBodyOper::type* otypes_;
    std::map<string, vector<string> > prefixes_;
    int num_types_;
    DensityMap dens_pairs_;
    DensityMapIterator mapit_;
    Ref<SCMatrixKit> kit_;
    RefSCDimension dim1_;
    RefSCDimension dim2_;
    int* quartets_processed_;



public:

    ComputeThread(
            int num,
            const Ref<TwoBodyIntDescr>& intdescr,
            const Ref<ThreadLock>& lock,
            const Ref<GaussianBasisSet>& bs1,
            const Ref<GaussianBasisSet>& bs2,
            const Ref<GaussianBasisSet>& bs3,
            const Ref<GaussianBasisSet>& bs4,
            TwoBodyOper::type* otypes,
            std::map<string, vector<string> > prefixes,
            int num_types,
            DensityMap& dens_pairs,
            Ref<SCMatrixKit>& kit,
            int* quartets_processed
    )
    {
        basis1_ = bs1;
        basis2_ = bs2;
        basis3_ = bs3;
        basis4_ = bs4;
        quartets_processed_ = quartets_processed;

        lock_ = lock;
        threadnum_ = num;

        timer.enter("get inteval");
        inteval_ = intdescr->inteval();
        timer.exit("get inteval");

        otypes_ = otypes;
        prefixes_ = prefixes;
        num_types_ = num_types;
        dens_pairs_ = dens_pairs;
        kit_ = kit;
        mapit_ = dens_pairs.begin();
        int maxn1 = 0, maxn2 = 0;
        SymmSCMatrixPair first_pair = (*mapit_).second;
        dim1_ = first_pair.first.dim();
        dim2_ = first_pair.second.dim();
        (*quartets_processed) = 0;
        /*
        for(; mapit_ != dens_pairs_.end(); mapit_++){
        	RefSCDimension dim = (*mapit_).second.first->dim();
        	assert(dim.nonnull());
        	int n1 = (*mapit_).second.first.dim().n();
        	if(n1 > maxn1){
        		maxn1 = n1;
        		dim1_ = (*mapit_).second.first.dim();
        	}
        	int n2 = (*mapit_).second.second.dim().n();
        	if(n2 > maxn2){
        		maxn2 = n2;
        		dim2_ = (*mapit_).second.second.dim();
        	}
        }
        */
    }

    ~ComputeThread()
    {
        inteval_ = 0;
    }


    void run()
    {
    	//=========================================================//
    	// open the output file(s)
    	//ofstream* o = new ofstream[num_types_];
    	std::map<string, ofstream* > o;
    	for(mapit_ = dens_pairs_.begin(); mapit_ != dens_pairs_.end(); mapit_++){
    		string pair_name = (*mapit_).first;
    		o[pair_name] = new ofstream[num_types_];
			for_each(ity, num_types_){
				stringstream sstr;
				sstr << prefixes_[pair_name][ity] << msg->me() << "_" << threadnum_ << ".bin";
				o[pair_name][ity].open(sstr.str().c_str(), ios::binary | ios::out);
			}
    	}

    	//=========================================================//
    	// Create the DistShellPair object
        DistShellPair shellpairs(msg, thr->nthread(), threadnum_,
                lock_, basis1_, basis3_, dynamic);
        // Setup the print frequency of the DistShellPair object
        shellpairs.set_print_percent(10);
        if(quiet) shellpairs.set_print_percent(1000);
        if(debug) shellpairs.set_debug(1);

    	//=========================================================//
        // Compute the integrals assigned to this thread

        int sh1 = 0, sh3 = 0;
        int sh2, sh4, nsh2, nsh4;
        idx_t identifier[4];
        const double* buffers[num_types_];
        for_each(ity, num_types_){
        	buffers[ity] = inteval_->buffer(otypes_[ity]);
        }
        while(shellpairs.get_task(sh1, sh3)) {
        	if (debug) {
        		lock_->lock();
        		cout << "      Computing shell pair (" << sh1
        			 << ", " << sh3 << ") on node "  << msg->me()
        			 << ", thread " << threadnum_ << "." << endl;
        		lock_->unlock();
        	}
			nsh2 = basis2_->nshell();
			nsh4 = basis4_->nshell();
        	(*quartets_processed_) += nsh2*nsh4;
        	identifier[0] = (idx_t)sh1;
        	identifier[2] = (idx_t)sh3;
        	const int nbf1 = basis1_->shell(sh1).nfunction();
        	const int nbf3 = basis3_->shell(sh3).nfunction();

        	// Set up some stuff for each pair
        	idx_t max_idents[4*num_types_];
        	int nbfpairs = nbf1*nbf3;
        	vector<vector<RefSCMatrix> > halft(num_types_);
        	for_each(ity, num_types_){
        		vector<RefSCMatrix> tmpv(nbfpairs);
        		halft[ity] = tmpv;
        		for_each(ipair, nbfpairs){
        			halft[ity][ipair] = kit_->matrix(dim1_, dim2_);
					halft[ity][ipair].assign(0.0);
        		}
        	}
        	for(sh2 = 0; sh2 < nsh2; ++sh2){
        		for(sh4 = 0; sh4 < nsh4; ++sh4){
        			identifier[1] = (idx_t)sh2;
        			identifier[3] = (idx_t)sh4;
        			const int nbf2 = basis2_->shell(sh2).nfunction();
					const int bfoff2 = basis2_->shell_to_function(sh2);
					const int nbf4 = basis4_->shell(sh4).nfunction();
					const int bfoff4 = basis4_->shell_to_function(sh4);

        			// Compute the shell
        			if (debug) {
        				lock_->lock();
        				cout << "        Computing shell quartet (" << sh1
        					 << "," << sh2 << "|" << sh3 << "," << sh4
        					 << ") on node "  << msg->me() << ", thread "
        					 << threadnum_ << "." << endl;
        				lock_->unlock();
        			}
        			inteval_->compute_shell(sh1, sh2, sh3, sh4);


        			// Loop over basis functions and store
        			int nfunc = nbf1*nbf2*nbf3*nbf4;
        			for_each(ity, num_types_){
        				const double* buff = buffers[ity];
						int bf1234 = 0;
						DBG_MSG("::nbf1="<<nbf1<<", nbf3="<<nbf3)
						for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
							int bfpair = bf1*nbf3 + bf3;
							DBG_MSG("::  sh1="<<sh1<<", sh3="<<sh3<<", bfpair=" << bfpair << ", bf2="<<bf2<<", bf4="<<bf4<<", value="<<buff[bf1234]);
							halft[ity][bfpair].set_element(bfoff2+bf2, bfoff4+bf4, buff[bf1234]);
							bf1234++;
						}
        			} // end loop over types

        		} // end loop over index 4
        	} // end loop over index 2
        	// Now transform and find maxima
        	iterate(mapit_, dens_pairs_){
        		string pairname = (*mapit_).first;
        		RefSymmSCMatrix P1 = (*mapit_).second.first;
        		RefSymmSCMatrix P2 = (*mapit_).second.second;
				for_each(ity, num_types_){
					value_t tmpval[nbfpairs];
					for_each(ipair, nbfpairs){
						DBG_MSG("Transforming ints for bf pair " << ipair
									<< "/" << nbfpairs << " of shell pair (" << sh1 << ", " << sh3 << ")")

						// First quarter transform
						RefSCMatrix transtmp = P1 * halft[ity][ipair];
						DBG_MSG(":: First quarter transform done for " << pairname)

						// Second quarter transform
						RefSCMatrix myhalf = transtmp * P2;
						DBG_MSG(":: Second quarter transform done for " << pairname)

						// get the maximum absolute value
						// SCMatrix->maxabs doesn't work for some reason
						tmpval[ipair] = max_abs(myhalf);

						DBG_MSG(":: Found maxabs for " << pairname)
					}
					if (debug) {
						lock_->lock();
						cout << "        Writing maximum entries for " << pairname
							 << " for shell pair (" << sh1 << ", " << sh3
							 << ") on node " << msg->me() << ", thread "
							 << threadnum_ << "." << endl;
						lock_->unlock();
					}
					o[pairname][ity].write((char*)&identifier, 4*sizeof(idx_t));
					o[pairname][ity].write((char*)&tmpval, nbfpairs*sizeof(value_t));
				}
			} // End loop over density matrix pairs
        } // End while get task

    	//=========================================================//
        // close the output
    	for(mapit_ = dens_pairs_.begin(); mapit_ != dens_pairs_.end(); mapit_++){
    		string pair_name = (*mapit_).first;
			for_each(ity, num_types_){
				o[pair_name][ity].close();
			}
			delete[] o[pair_name];
    	}

    }

};


inline off_t
file_size(string file){
	ifstream i(file.c_str(), ios::in | ios::binary | ios::ate);
	off_t length = (off_t)i.tellg();
	i.close();
	return length;
}

string
file_size_string(string file){
	double length = (double)file_size(file);
	stringstream sstr;
	if(length < 8.e2)
		sstr << length << " B";
	else if(length < 8.e5)
		sstr << setprecision(2) << fixed << length/1024.0 << " KB";
	else if(length < 8.e8)
		sstr << setprecision(2) << fixed << length/1024.0/1024.0 << " MB";
	else
		sstr << setprecision(2) << fixed << length/1024.0/1024.0/1024.0 << " GB";
	string str = sstr.str();
	return str;
}

inline void
copy_buffer(ifstream& i, ofstream& o){

	assert(i.is_open());
	int timeout = 10;
	i.seekg(0, ios::end);
	long length = i.tellg();
	int count = 0;
	while(length < 0 && count < timeout) {
		count++;
		cout << "  Failed to get size of file...retrying..." << endl;
		sleep(1);
		i.seekg(0, ios::end);
		length = i.tellg();
	}
	assert(length >= 0);
	i.seekg(0, ios::beg);
	char* buffer = new char[buff_size];
	while(length > buff_size){
		i.read(buffer, buff_size);
		o.write(buffer, buff_size);
		length -= buff_size;
	}
	i.read(buffer, length);
	o.write(buffer, length);
	delete[] buffer;

}

void
write_header(ofstream& o, const Ref<GaussianBasisSet>& basis)
{
	// Write the int_size, the double size, endianness, and number of atoms
    int8_t big_endian;
    if ( htonl(47) == 47 ) {
      big_endian = (int8_t)1;
    } else {
      big_endian = (int8_t)0;
    }
    o.write((char*)&big_endian, sizeof(int8_t));
    int8_t int_size = (int8_t)sizeof(int);
    o.write((char*)&int_size, sizeof(int8_t));
    int8_t value_size;
    value_size = (int8_t)sizeof(value_t);
	o.write((char*)&value_size, sizeof(int8_t));
    int8_t max8 = (int8_t)max_only;
    o.write((char*)&max8, sizeof(int8_t));

    // Write a description of the basis
    int natoms = basis->ncenter();
    o.write((char*)&natoms, sizeof(int));
    int nshell = basis->nshell();
    o.write((char*)&nshell, sizeof(int));
    // note that at this point, we are 64-bit aligned

    // Write the info about each shell
    for_each(ish, nshell){
    	GaussianShell sh = basis->shell(ish);
    	// Write the center number
    	int tmp = basis->shell_to_center(ish);
    	o.write((char*)&tmp, sizeof(int));
    	// Write the number of functions
    	tmp = sh.nfunction();
    	o.write((char*)&tmp, sizeof(int));
    	// Write the maximum and minimum angular momentum
    	tmp = sh.max_am();
    	o.write((char*)&tmp, sizeof(int));
    	tmp = sh.min_am();
    	o.write((char*)&tmp, sizeof(int));
    	// TODO Write some sort of "average contraction exponent" here...
    }
}

void
compute_ints_threaded(
        const Ref<TwoBodyIntDescr>& intdescr,
        TwoBodyOper::type* otypes, string* descs, int num_types,
        string prefix, string tmpdir,
        const Ref<GaussianBasisSet>& bs13,
        const Ref<GaussianBasisSet>& bs24,
    	DensityMap dens_pairs,
    	const Ref<SCMatrixKit>& inkit = (Ref<SCMatrixKit>)0
)
{
    //============================================================
    // Setup

	timer.enter("compute threads setup");
    int me = msg->me();
    int nthr = thr->nthread();
    int* quartets_processed[msg->n()];
    for_each(iproc, msg->n()){
    	quartets_processed[iproc] = new int[nthr];
    }

    if(!quiet && me == MASTER){
    	cout << "Computing integrals ";
    	for_each(ity, num_types-1)
    		cout << descs[ity] << ", ";
    	cout << descs[num_types-1] << "..." << endl;
    }

    // Setup the worker threads
    Ref<ThreadLock> lock = thr->new_lock();
    std::map<string, vector<string> > tmp_prefixes;
    DensityMapIterator mapit = dens_pairs.begin();
    Ref<SCMatrixKit> kit;
    if(inkit.nonnull())
    	kit = inkit;
    else
    	kit = (*mapit).second.first->kit();
    for(; mapit != dens_pairs.end(); mapit++) {
    	string pair_name = (*mapit).first;
		for_each(ity, num_types) {
			tmp_prefixes[pair_name].push_back(tmpdir + "/" + pair_name + descs[ity] + "_");
		}
    }
    for_each(ithr,nthr){
    	timer.enter("create compute threads");
        ComputeThread* thread = new ComputeThread(
        		ithr, intdescr, lock,
        		bs13, bs24, bs13, bs24,
        		otypes, tmp_prefixes, num_types,
        		dens_pairs, kit, &(quartets_processed[me][ithr])
		);
    	timer.exit("create compute threads");
        thr->add_thread(ithr, thread);
    }
	timer.exit("compute threads setup");

    //============================================================
    // Run the threads
	timer.enter("run threads");
    thr->start_threads();
    thr->wait_threads();
	timer.exit("run threads");

    //============================================================
    // Gather the results
    // First gather on each node...
    timer.enter("gather results");
    timer.enter("node gather");
	if(!quiet && me == MASTER)
		cout << "  Gathering results on a per-node basis..." << endl;
    for(mapit = dens_pairs.begin(); mapit != dens_pairs.end(); mapit++) {
    	string pair_name = (*mapit).first;
		for_each(ity, num_types){
			stringstream sstr;
			sstr << tmp_prefixes[pair_name][ity] << msg->me() << ".bin";
			string outfile = sstr.str();
			//const char* filename = sstr.str().c_str();
			ofstream o(outfile.c_str(), ios::out | ios::binary);
			// Loop over the temporary files created by the threads
			ifstream i;
			for_each(ithr, nthr){
				stringstream sstr;
				sstr << tmp_prefixes[pair_name][ity] << msg->me() << "_" << ithr << ".bin";
				string infile = sstr.str();
				if(debug){
					cout << "    Node " << me << " had file " << infile
						 << " for " << descs[ity] << " of size "
						 << file_size_string(infile) << endl;
				}
				i.open(infile.c_str(), ios::in | ios::binary);
				copy_buffer(i, o);
				i.close();
				remove(infile.c_str());
			}
			o.close();
		}
    }
	timer.exit("node gather");
    //============================================================
	timer.enter("sync");
	//if(!quiet && me == MASTER)
	//   cout << "  Waiting for all nodes to finish aggregation..." << endl;
	msg->sync();
	timer.exit("sync");
    //============================================================
	// Now gather across multiple nodes.  This could be done in a
	//   tree fashion if it becomes a bottleneck...
	timer.enter("master gather");
	if(me == MASTER) {
		if(!quiet) cout << "  Gathering results from all nodes..." << endl;
		ifstream i;
		for(mapit = dens_pairs.begin(); mapit != dens_pairs.end(); mapit++) {
			string pair_name = (*mapit).first;
			for_each(ity, num_types){
				string fname = prefix + pair_name + descs[ity] + ".bin";
				ofstream o(fname.c_str(), ios::out | ios::binary);
				write_header(o, bs13);
				for_each(inode, msg->n()){
					stringstream sstr;
					sstr << tmp_prefixes[pair_name][ity] << inode << ".bin";
					string infile = sstr.str();
					if(debug){
						cout << "    Node " << inode << " had file " << infile
							 << " for " << descs[ity] << " of size "
							 << file_size_string(infile) << endl;
					}
					i.open(infile.c_str(), ios::in | ios::binary);
					copy_buffer(i, o);
					i.close();
					remove(sstr.str().c_str());
				}
				o.close();
				if(!quiet) cout << "    Collected complete file for " << pair_name << descs[ity] << " of size " << file_size_string(fname.c_str()) << endl;
			}
		}
		// Print the shell pair distribution
		if(!quiet){
			int n_on_line = 0;
			const int n_per_line = 10;
			cout << "Pair distribution:" << setw(5) << endl;
			for_each(iproc,msg->n(), ithr,nthr){
				cout << iproc << " " << ithr << endl;
				if(n_on_line == n_per_line){
					cout << endl;
					n_on_line = 0;
				}
				cout << setw(7) << quartets_processed[iproc][ithr] << " ";
				n_on_line++;
			}
			cout << setw(0) << endl;
		}
	}
	timer.exit("master gather");
	timer.exit("gather results");


    //============================================================
    // Cleanup
    thr->delete_threads();
    for_each(iproc, msg->me()){
    	delete[] quartets_processed[iproc];
    }

}

int
main(int argc, char** argv) {

    //#############################################################################//
    // Setup
    //=========================================================//
    // Read the input file
    
    char *infile = new char[256];
    if (argc == 2) {
        infile = argv[1];
    } else {
        sprintf(infile, "pair_input.dat");
    }

    Ref<KeyVal> pkv(new ParsedKeyVal(infile));
    Ref<KeyVal> keyval(new PrefixKeyVal(pkv, ":test"));

    //=========================================================//
    // Create the message group

    msg = MessageGrp::initial_messagegrp(argc,argv);
    if (msg.null()) {
		#if HAVE_MPI_H
        msg = new MTMPIMessageGrp();
		#else
        msg = new ProcMessageGrp();
        low_comm_mode = true;
		#endif
    }
    MessageGrp::set_default_messagegrp(msg);
    n = msg->n();
    int me = msg->me();

    //=========================================================//
    // Create the timer
    
    Ref<RegionTimer> tim = new ParallelRegionTimer(msg, "AO Ints", 1, 1);
    timer = Timer(tim);

    //=========================================================//
    // Get the basis and the molecule from the input

    timer.enter("setup");
    timer.enter("input");
    Ref<GaussianBasisSet> obs = require_dynamic_cast<GaussianBasisSet*>(
            keyval->describedclassvalue("basis").pointer(), "main\n");
    Ref<GaussianBasisSet> auxbs = require_dynamic_cast<GaussianBasisSet*>(
            keyval->describedclassvalue("aux_basis").pointer(), "main\n");
    Ref<Molecule> mol = obs->molecule();
    print_distribution = keyval->booleanvalue("print_distribution");
    quartets_per_task = keyval->intvalue("quartets_per_task");
    num_masters = keyval->intvalue("num_masters", KeyValValueint(std::max(n/24, 1)));
    debug = keyval->booleanvalue("debug");
    verbose = keyval->booleanvalue("verbose");
    quiet = keyval->booleanvalue("quiet");
    dynamic = keyval->booleanvalue("dynamic", KeyValValueboolean(true));
    bool do_eri = keyval->booleanvalue("do_eri", KeyValValueboolean(true));
    bool do_f12 = keyval->booleanvalue("do_f12", KeyValValueboolean(true));
    bool do_f12g12 = keyval->booleanvalue("do_f12g12", KeyValValueboolean(true));
    bool do_f12sq = keyval->booleanvalue("do_f12sq", KeyValValueboolean(true));
    bool do_dblcomm = keyval->booleanvalue("do_dblcomm", KeyValValueboolean(true));
    bool do_ri = keyval->booleanvalue("do_ri", KeyValValueboolean(true));

    // Figure out where to put scratch files
    char* scratch_dir = getenv("SCRATCH");
    if(scratch_dir == 0){
    	scratch_dir = new char[32];
        sprintf(scratch_dir, "/tmp");
    }
    char* job_id = getenv("PBS_JOBID");
    int pid;
    if(job_id == 0 && me == MASTER){
    	// get unique process id.  (Only if we can't get a JOBID)
    	pid_t mypid = getpid();
    	pid = (int)mypid;
    }
    // broadcast the process id
    msg->bcast(&pid, sizeof(pid_t), 0);
	{
    	stringstream sstr;
		sstr << pid;
		job_id = const_cast<char*>(sstr.str().c_str());
	}
	// Get the temporary directory as a string
    stringstream sstr;
    sstr << scratch_dir << "/" << job_id;
    string tmpdir = sstr.str();
    if(!quiet){
    	if((!debug && me == MASTER) || debug){
    		cout << "Using temporary directory " << tmpdir << endl;
    	}
    }

    // Create the temporary directory
    if(me == MASTER){
    	int result = mkdir(tmpdir.c_str(), 0777);
    	assert(result==0);
    }


    timer.exit("input");
    
    //=========================================================//
    // Create the thread group object
    
    timer.enter("threadgrp");
    thr = ThreadGrp::initial_threadgrp(argc, argv);
    if (thr.null()) {
    	thr = ThreadGrp::get_default_threadgrp();
    }
    ThreadGrp::set_default_threadgrp(thr);
    int nthr = thr->nthread();
    if (!quiet && me == MASTER) {
        cout << "Running on " << n << " node" 
             << (n == 1 ? "" : "s") << " with " << nthr << " thread"
             << (nthr == 1 ? "" : "s") << " per node." << endl;
    }
    timer.exit("threadgrp");

    //=========================================================//
    // Construct an array of the number of basis functions
    // per center, which will be used in writing the binary files


    // Come up with a name for the files
    string basname = obs->name();
    string abasname = auxbs->name();
    string output_dir = keyval->stringvalue("output_dir");
    string molname = keyval->stringvalue("molname");
    string prefix = output_dir + "/" + molname + "_" + basname + "_" + abasname + "_";
    prefix += "allmax_";

    timer.exit("setup");

    //=========================================================//
    // Construct the CLHF object

    Ref<CLHF> hf = require_dynamic_cast<CLHF*>(
            keyval->describedclassvalue("hf").pointer(), "main\n");
    Ref<Integral> integral = hf->integral();
    Ref<PetiteList> pl;

    RefSymmSCMatrix P, Q;
    bool need_P = true;
    bool need_Q = true;
    bool need_O = true;

    // Only get P and Q if we need them.
    if(need_P || need_Q || need_O){
    	// Get P
    	timer.enter("HF");
    	pl = integral->petite_list();
		RefSymmSCMatrix Pso = hf->density();
		timer.exit("HF");
		Pso.scale(0.5);
		P = pl->to_AO_basis(Pso);

		if(need_Q){
			// Get Q
			timer.enter("Q");
			RefSymmSCMatrix S = hf->overlap();
			RefSymmSCMatrix Sinv = S.i();
			RefSymmSCMatrix Qso = Sinv - Pso;
			RefSymmSCMatrix Q = pl->to_AO_basis(Qso);
			timer.exit("Q");
		}
    }

    //=========================================================//
    // make R12WavefunctionWorld
    //
    timer.enter("build cf");
	Ref<WavefunctionWorld> world = require_dynamic_cast<WavefunctionWorld*>(
			keyval->describedclassvalue("world").pointer(), "main\n");
	Ref<RefWavefunction> refwfn = new SD_RefWavefunction(world, hf);
	Ref<R12WavefunctionWorld> r12world = new R12WavefunctionWorld(keyval, refwfn);
	Ref<R12Technology> r12tech = r12world->r12tech();
	Ref<OrbitalSpace> ribs_space = r12world->ribs_space();
	Ref<GaussianBasisSet> ribs = ribs_space->basis();
	Ref<R12Technology::CorrelationFactor> cf = r12tech->corrfactor();
	Ref<Integral> ri_integral = ribs_space->integral();
	Ref<OrbitalSpace> cabs_space = r12world->cabs_space(Alpha);
    timer.exit("build cf");



    //#############################################################################//

    if(do_eri || do_f12){
    	timer.enter("ERI/F12/F12G12");
    	int num_types = do_eri + do_f12;
    	TwoBodyOper::type otypes[num_types];
    	string* descs = new string[num_types];
    	int ity = 0;
    	if(do_eri){
    		otypes[ity] = cf->tbint_type_eri();
    		descs[ity] = "g";
    		++ity;
    	}
    	if(do_f12){
    		otypes[ity] = cf->tbint_type_f12();
    		descs[ity] = "F";
    		++ity;
    	}
    	// TODO Compute the integrals only once
    	Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
    	DensityMap dens_pairs;
    	dens_pairs["PP"] = SymmSCMatrixPair(P, P);
    	dens_pairs["QQ"] = SymmSCMatrixPair(Q, Q);
    	dens_pairs["PQ"] = SymmSCMatrixPair(P, Q);
    	compute_ints_threaded(
    			descr,
    			otypes, descs, num_types,
    			prefix, tmpdir,
    			obs, obs,
    			dens_pairs
    	);
    	if(do_ri){
			// Get O
			timer.enter("O");
			RefSCMatrix Ccabs = cabs_space->coefs();
			RefSymmSCMatrix O(Ccabs->rowdim(), Ccabs->kit());
			O.accumulate_symmetric_product(Ccabs);

			// Compute and tranform the integrals
    		ri_integral->set_basis(obs, ribs, obs, ribs);
    		DensityMap ri_pairs;
			ri_pairs["OO"] = SymmSCMatrixPair(O, O);
			Ref<TwoBodyIntDescr> ri_descr = cf->tbintdescr(ri_integral, 0);
			compute_ints_threaded(
					ri_descr,
					otypes, descs, num_types,
					prefix, tmpdir,
					obs, ribs,
					ri_pairs,
					Ccabs->kit()
			);

    	}
    	delete[] descs;
    	timer.exit("ERI/F12/F12G12");
    }

    //=========================================================//
    // Clean up

    timer.print();

    // Delete the temporary directory
    if(me == MASTER){
    	int error = remove(tmpdir.c_str());
    	if(error != 0){
    		cout << "WARNING: Could not delete temporary directory " << tmpdir << "." << endl;
    	}
    }

    return 0;

}
