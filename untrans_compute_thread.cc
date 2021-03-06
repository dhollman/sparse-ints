/*
 * half_compute_thread.cc
 *
 *  Created on: Feb 11, 2013
 *      Author: David Hollman
 */

#include "compute_thread.h"
#include "utils.h"
#include "histogram.h"


using namespace sparse_ints;
using namespace sc;
using namespace std;



/////////////////////////////////////////////////////////////////////////////
// UntransCommThread class

UntransComputeThread::UntransComputeThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const Ref<TwoBodyIntDescr>& intdescr,
		const Ref<ThreadLock>& lock,
		const Ref<GaussianBasisSet>& bs1,
		const Ref<GaussianBasisSet>& bs2,
		const Ref<GaussianBasisSet>& bs3,
		const Ref<GaussianBasisSet>& bs4,
		TwoBodyOper::type* otypes,
		vector<string> prefixes,
		int num_types,
		Ref<LocalSCMatrixKit>& kit,
		int* quartets_processed
) : ComputeThread(msg, thr, num, intdescr, lock, bs1, bs2, bs3, bs4, kit)
{
	quartets_processed_ = quartets_processed;
	(*quartets_processed_) = 0;

	otypes_ = otypes;
	prefixes_ = prefixes;
	num_types_ = num_types;
	kit_ = kit;
	int maxn1 = 0, maxn2 = 0;
}

void
UntransComputeThread::run()
{
	/*=========================================================*/
	/* Open the output file(s)	                          {{{1 */ #if fold_begin

	ofstream* o = new ofstream[num_types_];
	for_each(ity, num_types_){
		stringstream sstr;
		sstr << prefixes_[ity] << msg->me() << "_" << threadnum_ << ".bin";
		o[ity].open(sstr.str().c_str(), ios::binary | ios::out);
	}

	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
	/* Create the DistShellPair object and local vars     {{{1 */ #if fold_begin

	DistShellPair shellpairs(msg, thr->nthread(), threadnum_,
			lock_, basis1_, basis3_, opts.dynamic);
	// Setup the print frequency of the DistShellPair object
	shellpairs.set_print_percent(10);
	if(opts.quiet) shellpairs.set_print_percent(1000);
	if(opts.debug) shellpairs.set_debug(1);

	int sh1 = 0, sh3 = 0;
	int sh2, sh4, nsh2, nsh4;
	idx_t identifier[4];
	const double* buffers[num_types_];
	for_each(ity, num_types_){
		buffers[ity] = inteval_->buffer(otypes_[ity]);
	}
	bool bs1eqbs3 = basis1_ == basis3_;

	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
	/* Loop over tasks and compute untransformed ints     {{{1 */ #if fold_begin

	size_t total_mem_alloc = 0;
	size_t current_mem_alloc = 0;
	int nbftot2 = basis2_->nbasis();
	int nbftot4 = basis4_->nbasis();
	value_t nbftot24 = (value_t)(nbftot2*nbftot4);
	while(shellpairs.get_task(sh1, sh3)) {
		// TODO utilize permutational symmetry when possible
		// For now, just unroll the permutational symmetry (which
		//   only needs to be done if the basis sets are equal; if
		//   not, both permutations will be passed.)
		for(int iperm = 0; iperm < ((sh1==sh3 || !bs1eqbs3) ? 1 : 2); ++iperm){
			if(iperm == 1){
				int tmp = sh3;
				sh3 = sh1;
				sh1 = tmp;
			}
			nsh2 = basis2_->nshell();
			nsh4 = basis4_->nshell();
			(*quartets_processed_) += nsh2*nsh4;
			identifier[0] = (idx_t)sh1;
			identifier[2] = (idx_t)sh3;
			int nbf1 = basis1_->shell(sh1).nfunction();
			int nbf3 = basis3_->shell(sh3).nfunction();
			int nbfpairs = nbf1*nbf3;

			value_t* max_vals[num_types_];
			if(opts.out_type & MaxAbs){
				for_each(ity, num_types_){
					max_vals[ity] = allocate<value_t>(nbfpairs);
					for_each(ipair, nbfpairs)
						max_vals[ity][ipair] = -1.0;
				}
			}
			value_t* averages[num_types_];
			if(opts.out_type & Average){
				for_each(ity, num_types_){
					averages[ity] = allocate<value_t>(nbfpairs);
					for_each(ipair, nbfpairs)
						averages[ity][ipair] = -1.0;
				}
			}
			if(opts.out_type & Median){
				cout << opts.out_type << endl;
				assert(not_implemented);
			}
			if(opts.out_type & StdDev){
				cout << opts.out_type << endl;
				assert(not_implemented);
			}
			LogHistogram* histograms[num_types_];
			if(opts.out_type & Histogram){
				for_each(ity, num_types_){
					histograms[ity] = new LogHistogram[nbfpairs];
				}
			}

			for(sh2 = 0; sh2 < nsh2; ++sh2){
				for(sh4 = 0; sh4 < nsh4; ++sh4){
					identifier[1] = (idx_t)sh2;
					identifier[3] = (idx_t)sh4;

					const int nbf2 = basis2_->shell(sh2).nfunction();
					const int nbf4 = basis4_->shell(sh4).nfunction();

					// Compute the shell
					timer.enter("compute shell", threadnum_);
					inteval_->compute_shell(sh1, sh2, sh3, sh4);
					timer.exit("compute shell", threadnum_);
					(*quartets_processed_)++;

					// Loop over basis functions and store
					int nfunc = nbf1*nbf2*nbf3*nbf4;
					for_each(ity, num_types_){
						const double* buff = buffers[ity];
						if(opts.out_type == AllInts){
							// Cast to a value_t if value_t is not a double
							#if WRITE_AS_FLOAT
							value_t converted[nfunc];
							int bf1234 = 0;
							for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
								value_t val = (value_t)buff[bf1234];
								if(opts.use_fake_ints){
									val = (value_t)fake_int(opts.use_fake_ints);
								}
								converted[bf1234] = val;
								bf1234++;
							}
							assert(bf1234==nfunc);
							o[ity].write((char*)identifier, 4*sizeof(idx_t));
							o[ity].write((char*)converted, nfunc*sizeof(value_t));
							#else
							// Otherwise just write the buffer as is
							o[ity].write((char*)identifier, 4*sizeof(idx_t));
							if(opts.use_fake_ints){
								int bf1234 = 0;
								value_t converted[nfunc];
								for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
									converted[bf1234] = (value_t)fake_int(opts.use_fake_ints);
									bf1234++;
								}
								o[ity].write((char*)converted, nfunc*sizeof(double));
							}
							else{
								o[ity].write((char*)buff, nfunc*sizeof(double));
							}
							#endif
						}
						else{
							if(opts.out_type & MaxAbs){
								// Just keep track of the maximum value for the current
								//   basis function pair in the given shell
								int bf1234 = 0;
								for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
									int pairnum = bf1*nbf3 + bf3;
									value_t aval = fabs(buff[bf1234++]);
									if(aval > max_vals[ity][pairnum])
										max_vals[ity][pairnum] = aval;
								}
							}
							if(opts.out_type & Average){
								int bf1234 = 0;
								for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
									int pairnum = bf1*nbf3 + bf3;
									value_t aval = fabs(buff[bf1234++]);
									averages[ity][pairnum] += aval / nbftot24;
								}
							}
							if(opts.out_type & Histogram){
								int bf1234 = 0;
								for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
									int pairnum = bf1*nbf3 + bf3;
									histograms[ity][pairnum].insert(buff[bf1234++]);
								}
							}

						}
					} // end loop over types
				} // end loop over index 4
			} // end loop over index 2

			bool ident_written = false;
			if(opts.out_type & MaxAbs){
				for_each(ity, num_types_){
					o[ity].write((char*)identifier, 4*sizeof(idx_t));
					ident_written = true;
					o[ity].write((char*)(max_vals[ity]), nbfpairs*sizeof(value_t));
					deallocate(max_vals[ity]);
				}
			}
			if(opts.out_type & Average){
				for_each(ity, num_types_){
					if(!ident_written){
						o[ity].write((char*)identifier, 4*sizeof(idx_t));
						ident_written = true;
					}
					o[ity].write((char*)(averages[ity]), nbfpairs*sizeof(value_t));
					deallocate(averages[ity]);
				}
			}
			if(opts.out_type & Histogram){
				for_each(ity, num_types_){
					if(!ident_written){
						o[ity].write((char*)identifier, 4*sizeof(idx_t));
						ident_written = true;
					}
					for_each(ipair, nbfpairs){
						histograms[ity][ipair].write(o[ity]);
					}
					delete[] histograms[ity];
				}
			}


		} // End loop over permutations
	} // End while get task

	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
	/* Cleanup	                 	                      {{{1 */ #if fold_begin

	// Close the output files
	for_each(ity, num_types_){
		o[ity].close();
	}
	delete[] o;

	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
}




