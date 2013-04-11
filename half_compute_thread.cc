/*
 * half_compute_thread.cc
 *
 *  Created on: Feb 11, 2013
 *      Author: David Hollman
 */

#include "compute_thread.h"
#include "utils.h"

// macros
#define DBG_MSG(mymsg) \
	if(opts.debug){\
		print_lock->lock(); \
		dbg_out_ << mymsg << endl;\
		dbg_out_.flush(); \
		print_lock->unlock(); \
	}

using namespace sparse_ints;
using namespace sc;
using namespace std;



/////////////////////////////////////////////////////////////////////////////
// HalfTransCommThread class

HalfTransComputeThread::HalfTransComputeThread(
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
		Ref<LocalSCMatrixKit>& kit,
		int* quartets_processed
) : ComputeThread(num, intdescr, lock, bs1, bs2, bs3, bs4, kit)
{
	quartets_processed_ = quartets_processed;
	(*quartets_processed_) = 0;

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

}

void
HalfTransComputeThread::run()
{
	//=========================================================//
	// open the output file(s)
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
			lock_, basis1_, basis3_, opts.dynamic);
	// Setup the print frequency of the DistShellPair object
	shellpairs.set_print_percent(10);
	if(opts.quiet) shellpairs.set_print_percent(1000);
	if(opts.debug) shellpairs.set_debug(1);

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
		//DBG_MSG("Computing shell pair (" << sh1 << ", " << sh3 << ").")
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
				if (opts.debug) {
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
					//tmpval[ipair] = max_abs<value_t>(myhalf);
					tmpval[ipair] = myhalf->maxabs();
				}
				if (opts.debug) {
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




