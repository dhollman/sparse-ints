/*
 * compute_thread.cc
 *
 *  Created on: Jan 7, 2013
 *      Author: David Hollman
 *
 */

#include "compute_thread.h"
#include "utils.h"

using namespace std;
using namespace sc;
using namespace sparse_ints;


/////////////////////////////////////////////////////////////////////////////
// ComputeThread class

ComputeThread::ComputeThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::Ref<sc::LocalSCMatrixKit>& kit
) : SparseIntsThread(msg, thr, num, bs1, bs2, bs3, bs4, kit)
{
	lock_ = lock;

	timer.enter("get inteval");
	inteval_ = intdescr->inteval();
	timer.exit("get inteval");
}



/////////////////////////////////////////////////////////////////////////////
// FullTransComputeThread class

FullTransComputeThread::FullTransComputeThread(
		// Superclass arguments
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		// Arguments for this class
		sc::TwoBodyOper::type otype,
		std::string prefix,
		RefSymmSCMatrix& P1,
		RefSymmSCMatrix& P2,
		RefSymmSCMatrix& P3,
	    RefSymmSCMatrix& P4,
		Ref<LocalSCMatrixKit>& kit,
		SendThread* send_thread,
		ReceiveThread* recv_thread,
		int* quartets_computed
) : ComputeThread(msg, thr, num, intdescr, lock, bs1, bs2, bs3, bs4, kit),
		P1_(P1),
		P2_(P2),
		P3_(P3),
		P4_(P4),
		bsdim1_(bs1->basisdim()),
		bsdim2_(bs2->basisdim()),
		bsdim3_(bs3->basisdim()),
		bsdim4_(bs4->basisdim())
{
	send_thread_ = send_thread;
	recv_thread_ = recv_thread;
	prefix_ = prefix;
	otype_ = otype;
	quartets_computed_ = quartets_computed;
	(*quartets_computed) = 0;
}

void
FullTransComputeThread::run(){

	DBG_MSG("Compute thread starting up.")

	//=========================================================//
	// Create the DistShellPair object
	DistShellPair shellpairs(msg, thr->nthread()-2, threadnum_, lock_, basis1_, basis2_, opts.dynamic);

	// Setup the print frequency of the DistShellPair object
	shellpairs.set_print_percent(10);
	if(opts.quiet) shellpairs.set_print_percent(1000);
	if(opts.debug) shellpairs.set_debug(1);

	//=========================================================//
	// Compute the integrals assigned to this thread and
	//   perform the first two quarter transformations
	int sh3 = 0, sh4 = 0;
	const int nsh1 = basis1_->nshell();
	const int nsh2 = basis2_->nshell();
	const int nsh3 = basis3_->nshell();
	const int nbf1tot = basis1_->nbasis();
	const int nbf2tot = basis2_->nbasis();
	const int nbf3tot = basis3_->nbasis();
	const int nbf4tot = basis4_->nbasis();

	idx_t identifier[4];

	// Compute stuff until there's no more work left
	timer.enter("trans1q2q", threadnum_);
	// NOTE: USING SYMMETRY DOESN'T WORK YET!!!
	bool ignore_symmetry = true;
	bool bs3_eq_bs4 = basis3_ == basis4_;
	bool bs1_eq_bs2 = basis1_ == basis2_;
	while(shellpairs.get_task(sh3, sh4)){
		DBG_MSG("Got task compute pair (" << sh3 << ", " << sh4 << ")");
		for_each(iperm, ((ignore_symmetry && bs3_eq_bs4 && sh3 != sh4) ? 2 : 1)){
			if(iperm == 1){
				int shtmp = sh3;
				sh3 = sh4;
				sh4 = shtmp;
			}
			int nbf3 = basis3_->shell(sh3).nfunction();
			int nbf4 = basis4_->shell(sh4).nfunction();
			int sh3begin = basis3_->shell_to_function(sh3);
			int sh4begin = basis4_->shell_to_function(sh4);
			int nbfpairs = nbf3*nbf4;
			const double* buffer = inteval_->buffer(otype_);

			// Initialize the matrices that will hold the 1Q transformed integrals
			SCMatrix* ints1q[nbfpairs];
			for_each(ipair, nbfpairs){
				ints1q[ipair] = kit_->matrix(bsdim1_, bsdim2_);
				ints1q[ipair]->assign(0.0);
			}
			// Do the transformation of index 1.
			//   Loop over all shells for indexes 1 and 2
			for_each(sh1,nsh1){
				int sh2max = (bs1_eq_bs2 && !ignore_symmetry) ? sh1+1 : nsh2;
				for_each(sh2,sh2max){
					(*quartets_computed_)++;

					// Compute the shell
					timer.enter("compute shell", threadnum_);
					inteval_->compute_shell(sh1, sh2, sh3, sh4);
					timer.exit("compute shell", threadnum_);

					// Get the number of functions in the respective shells
					int nbf1 = basis1_->shell(sh1).nfunction();
					int nbf2 = basis2_->shell(sh2).nfunction();

					// Get the relevant block of the density matrix for the first quarter transform
					int sh1begin = basis1_->shell_to_function(sh1);
					int sh2begin = basis2_->shell_to_function(sh2);

					timer.enter("1q", threadnum_);

					int bf1234 = 0;
					for_each(bf1,nbf1, bf2,nbf2){
						int idx1 = sh1begin+bf1, idx2 = sh2begin+bf2;
						for_each(bf3,nbf3, bf4,nbf4){
							int bfpair = bf3*nbf4 + bf4;
							double val = buffer[bf1234];
							if(opts.use_fake_ints){
								val = fake_int(opts.use_fake_ints);
							}
							for_each(ibf1,nbf1tot){
								ints1q[bfpair]->accumulate_element(ibf1, idx2, P1_.get_element(ibf1, idx1) * val);
							}
							if(bs1_eq_bs2 && sh1 != sh2 && !ignore_symmetry){
								for_each(ibf1,nbf1tot){
									// Utilize the M-N symmetry in (MN|RS)
									ints1q[bfpair]->accumulate_element(ibf1, idx1, P1_(ibf1, idx2) * val);
								}
							}
							bf1234++;
						}
					}
					timer.exit("1q", threadnum_);
				} // end loop over all shell 2
			} // End loop over shell 1

			DBG_MSG("  Starting 2q transform")

			timer.enter("2q", threadnum_);
			for_each(iperm2, ((!ignore_symmetry && bs3_eq_bs4 && sh3 != sh4) ? 2 : 1)){
				if(iperm2 == 1){
					int shtmp = sh3;
					sh3 = sh4;
					sh4 = shtmp;
				}
				nbf3 = basis3_->shell(sh3).nfunction();
				nbf4 = basis4_->shell(sh4).nfunction();
				const int sh3begin = basis3_->shell_to_function(sh3);
				const int sh4begin = basis4_->shell_to_function(sh3);
				int nbfpairs = nbf3*nbf4;
				for_each(ish1,nsh1, jsh3,nsh3){
					const int nbf1 = basis1_->shell(ish1).nfunction();
					const int outernbf3 = basis3_->shell(jsh3).nfunction();
					const int bf1off = basis1_->shell_to_function(ish1);
					const int bf3off = basis3_->shell_to_function(jsh3);
					int nbfpair2q = nbf1*outernbf3;
					SCMatrix* ints2qS[nbfpair2q];
					RefSCDimension dim3(new SCDimension(nbf3));
					RefSCDimension dim4(new SCDimension(nbf4));
					for_each(bf1,nbf1, outerbf3,outernbf3){
						int ibf1 = bf1off+bf1;
						int jbf3 = bf3off+outerbf3;
						// Not a leaked pointer; deleted in send_thread_->distribute_shell_pair()
						SCMatrix* tmp2qS = kit_->matrix(bsdim2_, dim4);
						tmp2qS->assign(0.0);
						for_each(bf2n,nbf2tot, innerbf3,nbf3, bf4,nbf4){
							int idx3 = sh3begin + innerbf3;
							int pair1q = innerbf3*nbf4 + bf4;
							assert(ibf1<nbf1tot);
							assert(bf2n<nbf2tot);
							assert(ibf1>=0);
							assert(bf2n>=0);
							DBG_MSG("Doing I2Q[" << bf2n << ", " << bf4 << "] += P3[" << jbf3 << "," << idx3 << "] * ints1q[" << pair1q << "][" << ibf1 << "," << bf2n << "]");
							tmp2qS->accumulate_element(bf2n, bf4, P3_.get_element(jbf3, idx3) * ints1q[pair1q]->get_element(ibf1, bf2n));
							/*SCMatrix* matnum = ints1q[pair1q];
							DBG_MSG("  Matrix ints1q has nrow = " << matnum->nrow() << " and ncol = " << matnum->ncol());
							double q1part = matnum->get_element(ibf1, bf2n);
							double p3part = P3_.get_element(jbf3, idx3);
							tmp2qS->accumulate_element(bf2n, bf4, p3part * q1part);*/
							DBG_MSG("...success!");
						}
						int ipair2q = bf1*outernbf3 + outerbf3;
						assert(ipair2q < nbfpair2q);
						ints2qS[ipair2q] = tmp2qS;
					}
					// This subroutine transfers the data and deletes the SCMatrix objects in ints2qS,
					//   as well as posting the send request
					send_thread_->distribute_shell_pair(ints2qS, ish1, jsh3, sh4, threadnum_);
				}
			}
			timer.exit("2q", threadnum_);
			for_each(ipair, nbfpairs){
				delete ints1q[ipair];
			}
		} // end loop over permutations
	} // End while get_task
	timer.exit("trans1q2q", threadnum_);


	//=========================================================//

	if(msg->me() == MASTER){
		// There's no reason the master node can't participate in the 3rd and 4th quarter transforms
		//   This would be more trouble than it's worth though.  It doesn't gain much in terms of
		//   aggregate memory except in the case of small computations, in which case it doesn't
		//   matter anyway.
		DBG_MSG("Master compute thread shutting down");
		return;
	}


	timer.enter("half-trans sync", threadnum_);
	// Tell the send thread that we are done computing data
	DBG_MSG("No more pairs, sending ComputeThreadDone message.");
	msg->sendt(msg->me(), NeedsSend, threadnum_);

	// Wait until all of the data is ready and available
	DBG_MSG("No more pairs, sending HaveAllData message.");
	int garbage = 0;
	msg->recvt(msg->me(), HaveAllData, garbage);
	DBG_MSG("Waking up after waiting for comm threads.")
	timer.exit("half-trans sync", threadnum_);



	//=========================================================//
	// Open the output file for writing
	timer.enter("trans3q4q", threadnum_);
	char name[512];
	sprintf(name, "%s%d_%d.bin", prefix_.c_str(), msg->me(), threadnum_);
	ofstream o(name, ios::binary | ios::out);

	//=========================================================//

	// Get our thread's next shell assignment
	int sh1 = 0;
	while(recv_thread_->get_my_next_pair(sh1, sh3)){
		DBG_MSG("Got local pair (" << sh1 << "," << sh3 << ") to work on");
		const int nbf1 = basis1_->shell(sh1).nfunction();
		const int nbf3 = basis3_->shell(sh3).nfunction();
		IntPair sh13(sh1, sh3);
		SCMatrix** ints24 = recv_thread_->my_ints2q[sh13];
		int nbfpairs = nbf1*nbf3;
		idx_t identifier[4];
		value_t maxvals[nbfpairs];
		cr->consume_memory(nbfpairs*sizeof(value_t));
		value_t averages[nbfpairs];
		cr->consume_memory(nbfpairs*sizeof(value_t));
		value_t medians[nbfpairs];
		cr->consume_memory(nbfpairs*sizeof(value_t));
		value_t stddevs[nbfpairs];
		cr->consume_memory(nbfpairs*sizeof(value_t));

		// Write the identifier beforehand, since it applies to both the allints and maxabs schemes
		identifier[0] = (idx_t)sh1;
		identifier[2] = (idx_t)sh3;
		o.write((char*)&identifier, 4*sizeof(idx_t));
		for_each(bf1,nbf1, bf3,nbf3){
			int ipair = bf1*nbf3 + bf3;
			RefSCMatrix ints3q = kit_->matrix(bsdim2_, bsdim4_);
			RefSCMatrix ints4q = kit_->matrix(bsdim2_, bsdim4_);
			assert(P2_.nonnull());
			ints3q = P2_ * ints24[ipair];
			ints4q = ints3q * P4_;
			if(opts.out_type == AllInts){
				// write all ints
				value_t vals[nbf2tot*nbf4tot];
				for_each(idx2,nbf2tot, idx4,nbf4tot){
					vals[idx2*nbf4tot + idx4] = (value_t)ints4q.get_element(idx2, idx4);
				}
				o.write((char*)&vals, nbf2tot*nbf4tot*sizeof(value_t));
			}
			else {
				// TODO average, stddev, median, pair max etc (all at once, and write all files in one go)
				if(opts.out_type & MaxAbs){
					maxvals[ipair] = ints4q->maxabs();
				}
				//--------------------//
				if(opts.out_type & Average){
					averages[ipair] = average<value_t>(ints4q);
					if(opts.out_type & StdDev){
						stddevs[ipair] = stddev<value_t>(ints4q, averages[ipair]);
					}
				}
				//--------------------//
				else if(opts.out_type & StdDev){
					stddevs[ipair] = stddev<value_t>(ints4q);
				}
				//--------------------//
				if(opts.out_type & Median){
					medians[ipair] = median<value_t>(ints4q);
				}
				//--------------------//
				if(opts.out_type & Histogram){
					assert(not_implemented);
				}

			}
		}
		if(opts.out_type & MaxAbs){
			o.write((char*)maxvals, nbfpairs*sizeof(value_t));
		}
		if(opts.out_type & Average){
			o.write((char*)averages, nbfpairs*sizeof(value_t));
		}
		if(opts.out_type & StdDev){
			o.write((char*)stddevs, nbfpairs*sizeof(value_t));
		}
		if(opts.out_type & Median){
			o.write((char*)medians, nbfpairs*sizeof(value_t));
		}
		if(opts.out_type & Histogram){
			assert(not_implemented);
		}
		// Release memory since the following variables are going out of scope:
		// maxvals
		cr->release_memory(nbfpairs*sizeof(value_t));
		// averages
		cr->release_memory(nbfpairs*sizeof(value_t));
		// medians
		cr->release_memory(nbfpairs*sizeof(value_t));
		// stddevs
		cr->release_memory(nbfpairs*sizeof(value_t));
	}
	o.close();

	timer.exit("trans3q4q", threadnum_);

	DBG_MSG("Compute thread shutting down.");

}

