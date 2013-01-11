/*
 * compute_thread.cc
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "compute_thread.h"
#include "utils.h"


using namespace std;
using namespace sc;
using namespace sparse_ints;

// Not actually used yet
//static size_t FullTransCommThread::max_queue_size = 10 * 1024 * 1024 * sizeof(double);



/////////////////////////////////////////////////////////////////////////////
// ComputeThread class

ComputeThread::ComputeThread(
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4
) :
		basis1_(bs1),
		basis2_(bs2),
		basis3_(bs3),
		basis4_(bs4)
{
	lock_ = lock;
	threadnum_ = num;

	timer.enter("get inteval");
	inteval_ = intdescr->inteval();
	timer.exit("get inteval");

}


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
		Ref<SCMatrixKit>& kit,
		int* quartets_processed
) : ComputeThread(num, intdescr, lock, bs1, bs2, bs3, bs4)
{
	quartets_processed_ = quartets_processed;
	(*quartets_processed) = 0;

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
		if (opts.debug) {
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
					// SCMatrix->maxabs doesn't work for some reason
					tmpval[ipair] = max_abs<value_t>(myhalf);

					DBG_MSG(":: Found maxabs for " << pairname)
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


/////////////////////////////////////////////////////////////////////////////
// FullTransCommThread class

FullTransComputeThread::FullTransComputeThread(
		// Superclass arguments
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
		Ref<SCMatrixKit>& kit,
		FullTransCommThread* send_thread,
		FullTransCommThread* recv_thread
) :
		ComputeThread(num, intdescr, lock, bs1, bs2, bs3, bs4),
		P1_(P1),
		P2_(P2),
		P3_(P3),
		P4_(P4),
		bsdim1_(new SCDimension(bs1->nbasis(), 1)),
		bsdim2_(new SCDimension(bs2->nbasis(), 1)),
		bsdim3_(new SCDimension(bs3->nbasis(), 1)),
		bsdim4_(new SCDimension(bs4->nbasis(), 1))
{
	send_thread_ = send_thread;
	recv_thread_ = recv_thread;
	prefix_ = prefix;
	otype_ = otype;
	bsdim1_->blocks()->set_subdim(0, new SCDimension(bs1->nbasis()));
	bsdim2_->blocks()->set_subdim(0, new SCDimension(bs2->nbasis()));
	bsdim3_->blocks()->set_subdim(0, new SCDimension(bs3->nbasis()));
	bsdim4_->blocks()->set_subdim(0, new SCDimension(bs4->nbasis()));
	kit_ = kit;
}


void
FullTransComputeThread::run(){

	//=========================================================//
	// Create the DistShellPair object
	DistShellPair shellpairs(msg, thr->nthread(), threadnum_, lock_, basis1_, basis3_, opts.dynamic);

	// Setup the print frequency of the DistShellPair object
	shellpairs.set_print_percent(10);
	if(opts.quiet) shellpairs.set_print_percent(1000);
	if(opts.debug) shellpairs.set_debug(1);

	//=========================================================//
	// Compute the integrals assigned to this thread and
	//   perform the first two quarter transformations
	int sh1 = 0, sh2 = 0;
	const int nsh3 = basis3_->nshell();
	const int nsh4 = basis4_->nshell();
	const int nbf3tot = basis3_->nbasis();
	const int nbf4tot = basis4_->nbasis();

	idx_t identifier[4];

	// Compute stuff until there's no more work left
	while(shellpairs.get_task(sh1, sh2)){
		const int nbf1 = basis1_->shell(sh1).nfunction();
		const int nbf2 = basis2_->shell(sh2).nfunction();
		int nbfpairs = nbf1*nbf2;
		const double* buffer = inteval_->buffer(otype_);

		// Initialize the matrices that will hold the 1Q transformed integrals
		vector<RefSCMatrix> ints1q;
		for_each(ipair, nbfpairs){
			ints1q[ipair] = kit_->matrix(bsdim3_, bsdim4_);
			ints1q[ipair].assign(0.0);
		}
		// Do the transformation of index 3.
		//   Loop over all shells for indexes 3 and 4
		for_each(sh3,nsh3, sh4,nsh4){
			// Compute the shell
			inteval_->compute_shell(sh1, sh2, sh3, sh4);

			// Get the number of functions in the respective shells
			const int nbf3 = basis3_->shell(sh3).nfunction();
			const int nbf4 = basis4_->shell(sh4).nfunction();

			// Get the relevant block of the density matrix for the first quarter transform
			int sh3begin = basis3_->shell_to_function(sh3);
			int sh4begin = basis4_->shell_to_function(sh4);
			RefSCMatrix P3part = P3_.get_subblock(0, nbf3tot, sh3begin, sh3begin + nbf3);

			// Create the matrix to hold the untransformed integrals temporarily
			RefSCDimension bfdim3(new SCDimension(nbf3, 1));
			bfdim3->blocks()->set_subdim(0, new SCDimension(nbf3));
			RefSCDimension bfdim4(new SCDimension(nbf4, 1));
			bfdim4->blocks()->set_subdim(0, new SCDimension(nbf4));
			RefSCMatrix g12 = kit_->matrix(bfdim3, bfdim4);

			int bf1234 = 0;
			for_each(bf1,nbf1, bf2,nbf2){
				for_each(bf3,nbf3, bf4,nbf4){
					g12.set_element(bf3, bf4, buffer[bf1234]);
					bf1234++;
				}
				// Better:
				// g12->assign(&buffer(bf1*nbf2*nbf3*nbf4 + bf2*nbf3*nbf4))

				// accumulate the first quarter transform
				int pairnum = bf1*nbf2 + bf2;
				ints1q[pairnum].accumulate_subblock(P3part * g12, 0, nbf3tot, sh4begin, sh4begin + nbf4);
			}
		}

		// Initialize the matrices that will hold the 2Q transformed integrals
		vector<RefSCMatrix> ints2q;
		for_each(ipair, nbfpairs){
			ints2q[ipair] = kit_->matrix(bsdim3_, bsdim4_);
			ints2q[ipair] = ints1q[ipair] * P4_;
		}

		// Now tell the comm thread to send the data to the node assigned to handle it.
		send_thread_->distribute_shell_pair(ints2q, sh1, sh2, nbf1, nbf2, nbf3tot, nbf4tot, this);

	} // End while get_task

	//=========================================================//

	// Tell the send thread that we are done computing data
	int thridx = threadnum_ - 2;
	msg->sendt(msg->me(), ComputeThreadDone, thridx);

	// Wait until all of the data is ready and available
	int garbage = 0;
	msg->recvt(msg->me(), HaveAllData, garbage);


	//=========================================================//
	// Open the output file for writing
	ofstream o;
	{
		stringstream sstr;
		sstr << prefix_ << msg->me() << "_" << threadnum_ << ".bin";
		o.open(sstr.str().c_str(), ios::binary | ios::out);
	}

	//=========================================================//

	// Get our thread's next shell assignment
	int sh3 = 0, sh4 = 0;
	while(recv_thread_->get_my_next_pair(sh3, sh4)){
		const int nbf3 = basis3_->shell(sh3).nfunction();
		const int nbf4 = basis4_->shell(sh4).nfunction();
		IntPair sh34(sh3, sh4);
		vector<RefSCMatrix> ints34 = recv_thread_->my_ints2q[sh34];
		// TODO write only unique pairs
		int nbfpairs = nbf3*nbf4;
		idx_t identifier[4];
		value_t maxvals[nbfpairs];
		for_each(bf3,nbf3, bf4,nbf4){
			int ipair = bf3*nbf4 + bf4;
			RefSCMatrix ints3q = P1_ * ints34[ipair];
			RefSCMatrix ints4q = ints3q * P2_;
			if(opts.max_only){
				//value_t maxval = max_abs<value_t, idx_t>(ints4q, idx_t[2], idx_t[3]);
				maxvals[ipair] = max_abs<value_t>(ints4q);
			}
			else{
				assert(not_implemented);
			}
			identifier[0] = bf3;
			identifier[1] = bf4;
			o.write((char*)&identifier, 4*sizeof(idx_t));
			o.write((char*)maxvals, nbfpairs*sizeof(value_t));
		}
	}

}


/////////////////////////////////////////////////////////////////////////////
// FullTransCommThread class

FullTransCommThread::FullTransCommThread(
	const sc::Ref<sc::ThreadLock>& comm_lock,
	const sc::Ref<sc::ThreadLock>& queue_lock,
	const sc::Ref<sc::GaussianBasisSet>& bs1,
	const sc::Ref<sc::GaussianBasisSet>& bs2,
	const sc::Ref<sc::GaussianBasisSet>& bs3,
	const sc::Ref<sc::GaussianBasisSet>& bs4,
	int thr_type,
	Ref<SCMatrixKit>& kit
) :
		bf_per_node_(),
		pair_assignments_(),
		pairs_sorted_(),
		task_queue_(),
		basis1_(bs1),
		basis2_(bs2),
		basis3_(bs3),
		basis4_(bs4),
		kit_(kit)
{
	comm_lock_ = comm_lock;
	queue_lock_ = queue_lock;
	queue_size_ = 0;
	thread_type_ = thr_type;
	my_ints2q_complete_ = false;
	sorted_pairs_position_ = 0;

	if(msg->me() == MASTER){
		for(int inode = 1; inode < msg->n(); ++inode){
			bf_per_node_.push(IntPair(0, inode));
		}
	}

}

void
FullTransCommThread::distribute_shell_pair(
	std::vector<RefSCMatrix> pair_mats,
	int sh1, int sh2, int nbf1, int nbf2, int nbf3tot, int nbf4tot,
	FullTransComputeThread* compute_thread
)
{
	// Note that this is called by the Compute thread, not the Comm thread!

	DataSendTask task;
	task.sh1 = sh1;
	task.sh2 = sh2;
	task.nbf1 = nbf1;
	task.nbf2 = nbf2;
	task.nbftot3 = nbf3tot;
	task.nbftot4 = nbf4tot;
	int stride = nbf3tot*nbf4tot;
	int ndata = nbf1*nbf2*stride;
	task.data = new double[ndata];
	for_each(bf1,nbf1, bf2,nbf2){
		int ipair = bf1*nbf2 + bf2;
		pair_mats[ipair].convert(task.data);
		task.data += stride;
	}

	// TODO Implement a queue size limit using e.g. pthread contitionals
	// Acquire the queue lock.
	queue_lock_->lock();

	// Push the task onto the queue
	task_queue_.push(task);
	queue_size_ += sizeof(DataSendTask) + ndata*sizeof(double);

	int send_code = -1;
	msg->sendt(msg->me(), NeedsSend, send_code);

	// Now release the lock
	queue_lock_->unlock();

}

int
FullTransCommThread::get_pair_assignment(int sh3, int sh4)
{
	if(msg->me() == MASTER){
		IntPair sh34(sh3, sh4);
		const int nbf3 = basis3_->shell(sh3).nfunction();
		const int nbf4 = basis4_->shell(sh4).nfunction();
		std::map<IntPair, int>::iterator item = pair_assignments_.find(sh34);
		if(item == pair_assignments_.end()){
			IntPair node_assignment = bf_per_node_.top();
			pair_assignments_[sh34] = node_assignment.second;
			node_assignment.first += nbf3*nbf4;
			return node_assignment.second;
		}
		else{
			return pair_assignments_[sh34];
		}
	}
	else{
		int sh34_sender[] = {sh3, sh4, msg->me()};
		// Lock to prevent overlapping threads from getting the wrong message.
		//   Could also (probably) be done with different typed messages for each thread...
		comm_lock_->lock();
		msg->sendt(MASTER, NeedPairAssignment, sh34_sender, 3);
		int assignment = 0;
		msg->recvt(MASTER, PairAssignment, assignment);
		comm_lock_->unlock();
		return assignment;
	}
}

bool
FullTransCommThread::get_my_next_pair(int& sh3, int& sh4){
	queue_lock_->lock();
	if(pairs_sorted_.size() == sorted_pairs_position_){
		queue_lock_->unlock();
		return false;
	}
	else{
		IntPair sh34 = pairs_sorted_[sorted_pairs_position_];
		sorted_pairs_position_++;
		queue_lock_->unlock();
		sh3 = sh34.first;
		sh4 = sh34.second;
		return true;
	}
}

void
FullTransCommThread::run()
{
	// This will not work if n = 1
	assert(msg->n() > 1);


	if(opts.dynamic and msg->me() == MASTER){
		// We're the manager for the comm threads; get to work assigning places to put pairs
		if(thread_type_ == FullTransCommThread::SendThread){
			// First, post a non-blocking recieve from each of the worker nodes
			//   of the HaveAllData message
			int done_nodes = 0;

			int n = msg->n();
			int sh34_sender[3];
			bool nodes_still_need_data = true;
			while(done_nodes < n-1){
				msg->recvt(MessageGrp::AnySender, NeedPairAssignment, sh34_sender, 3);
				if(sh34_sender[0] == -1){
					done_nodes++;
				}
				else{
					int node_assignment = get_pair_assignment(sh34_sender[0], sh34_sender[1]);
					msg->sendt(sh34_sender[2], PairAssignment, node_assignment);
				}
			}
		}
	}
	else{
		// We're the comm thread on a non master node, or static task distribution is enabled.
		// static not implemented yet
		assert(opts.dynamic);

		const int nsh3 = basis3_->nshell();
		const int nsh4 = basis4_->nshell();

		if(thread_type_ == FullTransCommThread::SendThread){
			bool tasks_done[thr->nthread()-2];
			for_each(ithr, thr->nthread()-2){
				tasks_done[ithr] = false;
			}
			while(true){

				// TODO this would be done better using e.g. pthread conditionals
				// Receive a message alerting us that there is communication work to do
				int task_code = -1;
				msg->recvt(msg->me(), NeedsSend, task_code);
				if(task_code >= 0){
					// This was a "done" message. Make note of it and continue if needed.
					tasks_done[task_code] = true;
					bool all_tasks_done = true;
					for_each(ithr, thr->nthread()-2){
						if(!tasks_done[ithr]){
							all_tasks_done = false;
							break;
						}
					}
					if(all_tasks_done)
						break;
					else
						continue;
				}

				// Get the task
				queue_lock_->lock();
				DataSendTask t = task_queue_.front();
				task_queue_.pop();
				queue_lock_->unlock();

				// Send the data
				// Put all of the data belonging to a given index 3, 4 pair together
				for_each(sh3,nsh3, sh4,nsh4){
					int dest_node = get_pair_assignment(sh3, sh4);

					int sh3off = basis3_->shell_to_function(sh3);
					int sh4off = basis4_->shell_to_function(sh4);

					const int nbf3 = basis3_->shell(sh3).nfunction();
					const int nbf4 = basis4_->shell(sh4).nfunction();

					int pair_ndata = t.nbf1 * t.nbf2 * nbf3 * nbf4;
					double send_data[pair_ndata];

					int idxs[5];
					idxs[0] = t.sh1;
					idxs[1] = t.sh2;
					idxs[2] = sh3;
					idxs[3] = sh4;
					idxs[4] = msg->me();
					int ibf = 0;
					double *data_ptr;
					// This could probably be done more efficiently (i.e. with a better stride etc.)
					for_each(bf1,t.nbf1, bf2,t.nbf2){
						int pair12 = bf1*t.nbf2 + bf2;
						data_ptr = &t.data[pair12];
						for_each(bf3,nbf3, bf4,nbf4){
							send_data[ibf] = data_ptr[(bf3+sh3off)*t.nbftot4 + sh4off + bf4];
							ibf++;
						}
					}
					msg->sendt(dest_node, IndexData, idxs, 5);
					msg->sendt(dest_node, PairData, send_data, pair_ndata);
				}
				delete[] t.data;
			}

			int me = msg->me();
			// Send a message to the receive thread letting it know that all of the compute threads are
			//   done and just waiting to have all of the data for the third and fourth quarter transformations.
			int done_msg[] = {-1, -1, -1, -1, msg->me()};
			msg->sendt(me, ComputeThreadDone, done_msg, 5);
			// Tell the master node we don't need any more pair assignments
			int pair_assignment_end[] = {-1, -1, me};
			msg->sendt(MASTER, NeedPairAssignment, pair_assignment_end, 3);
		}
		else if(thread_type_ == FullTransCommThread::ReceiveThread){
			bool nodes_done[msg->n()-1];
			for_each(inode, msg->n()-1){
				nodes_done[inode] = false;
			}

			while(true){

				int idxs[5];
				msg->recvt(MessageGrp::AnySender, IndexData, idxs, 5);
				if(idxs[0] == -1){
					nodes_done[idxs[4]-1] = true;
					bool all_nodes_done = true;
					for_each(inode, msg->n()-1){
						if(!nodes_done[inode]){
							all_nodes_done = false;
							break;
						}
					}
					if(all_nodes_done)
						break;
					else
						continue;
				}

				const int sh1 = idxs[0];
				const int sh2 = idxs[1];
				const int sh3 = idxs[2];
				const int sh4 = idxs[3];
				const int nbf1 = basis1_->shell(sh1).nfunction();
				const int nbf2 = basis2_->shell(sh2).nfunction();
				const int nbf3 = basis3_->shell(sh3).nfunction();
				const int nbf4 = basis4_->shell(sh4).nfunction();
				int pair_ndata = nbf1 * nbf2 * nbf3 * nbf4;
				int sh1off = basis1_->shell_to_function(sh1);
				int sh2off = basis2_->shell_to_function(sh2);

				double recv_data[pair_ndata];
				msg->recvt(idxs[4], PairData, recv_data, pair_ndata);

				IntPair sh34(sh3,sh4);
				vector<RefSCMatrix> ints34 = my_ints2q[sh34];

				if(ints34.empty()){
					for_each(bf3,nbf3, bf4,nbf4){
						RefSCMatrix tmp = kit_->matrix(basis1_->basisdim(), basis2_->basisdim());
						tmp.assign(0.0);
						ints34.push_back(tmp);
					}
				}
				int ibf = 0;
				for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
					int ipair = bf3*nbf4 + bf4;
					ints34[ipair].accumulate_element(sh1off+bf1, sh2off+bf2, recv_data[ibf]);
					ibf++;
				}
				assert(ibf==pair_ndata);
			}

			my_ints2q_complete_ = true;
			// Create a sorted list of my pairs decending by number of basis functions
			std::map<IntPair, vector<RefSCMatrix> >::iterator mapit;
			iterate(mapit, my_ints2q){
				pairs_sorted_.push_back((*mapit).first);
			}
			sorted_pairs_position_ = 0;
			bf_compare_ comparator(basis3_, basis4_);
			std::stable_sort(pairs_sorted_.begin(), pairs_sorted_.end(), comparator);


			// Now that all of the nodes are done sending us stuff,
			//   we can tell the compute threads to get back to work.
			int me = msg->me();
			for_each(ithr, thr->nthread()-2){
				// This would be better done with pthreads conditionals
				msg->sendt(msg->me(), HaveAllData, me);
			}
		}

	}

}

