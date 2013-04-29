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

#define QUEUE_MAX_ENABLED 0

size_t SendThread::max_queue_size = size_t(2e8);

#define assert_equal(a, b) \
	if(a != b){\
		print_lock->lock(); \
		cout << "Assertion failed on node " << msg->me() << "." << threadnum_ << ":  " << a << " == " << b << " returns false" << endl; \
		assert(a == b); \
		print_lock->unlock(); \
	}





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
			vector<RefSCMatrix> ints1q;
			for_each(ipair, nbfpairs){
				ints1q.push_back(kit_->matrix(bsdim1_, bsdim2_));
				ints1q[ipair].assign(0.0);
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
								ints1q[bfpair].accumulate_element(ibf1, idx2, P1_.get_element(ibf1, idx1) * val);
							}
							if(bs1_eq_bs2 && sh1 != sh2 && !ignore_symmetry){
								for_each(ibf1,nbf1tot){
									// Utilize the M-N symmetry in (MN|RS)
									ints1q[bfpair].accumulate_element(ibf1, idx1, P1_(ibf1, idx2) * val);
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
					vector<RefSCMatrix> ints2qS, ints2qR;
					RefSCDimension dim3(new SCDimension(nbf3));
					RefSCDimension dim4(new SCDimension(nbf4));
					for_each(bf1,nbf1, outerbf3,outernbf3){
						int ibf1 = bf1off+bf1;
						int jbf3 = bf3off+outerbf3;
						RefSCMatrix tmp2qS = kit_->matrix(bsdim2_, dim4);
						tmp2qS.assign(0.0);
						for_each(bf2n,nbf2tot, innerbf3,nbf3, bf4,nbf4){
							int idx3 = sh3begin + innerbf3;
							int pair1q = innerbf3*nbf4 + bf4;
							tmp2qS.accumulate_element(bf2n, bf4, P3_.get_element(jbf3, idx3) * ints1q[pair1q].get_element(ibf1, bf2n));
						}
						ints2qS.push_back(tmp2qS);
						assert_equal(ints2qS.size()-1, bf1*outernbf3 + outerbf3);
					}
					send_thread_->distribute_shell_pair(ints2qS, ish1, jsh3, sh4, threadnum_);
				}
			}
			timer.exit("2q", threadnum_);
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
		vector<RefSCMatrix> ints24 = recv_thread_->my_ints2q[sh13];
		int nbfpairs = nbf1*nbf3;
		idx_t identifier[4];
		value_t maxvals[nbfpairs];
		// Write the identifier beforehand, since it applies to both the allints and maxabs schemes
		identifier[0] = (idx_t)sh1;
		identifier[2] = (idx_t)sh3;
		o.write((char*)&identifier, 4*sizeof(idx_t));
		for_each(bf1,nbf1, bf3,nbf3){
			int ipair = bf1*nbf3 + bf3;
			RefSCMatrix ints3q = kit_->matrix(bsdim2_, bsdim4_);
			RefSCMatrix ints4q = kit_->matrix(bsdim2_, bsdim4_);
			assert(P2_.nonnull());
			assert(ints24[ipair].nonnull());
			ints3q = P2_ * ints24[ipair];
			ints4q = ints3q * P4_;
			if(opts.out_type == MaxAbs){
				// TODO average, stddev, median, pair max etc (all at once, and write all files in one go)
				maxvals[ipair] = ints4q->maxabs();
				//maxvals[ipair] = mean<value_t>(ints24[ipair]);
			}
			else if(opts.out_type == AllInts){
				// write all ints
				value_t vals[nbf2tot*nbf4tot];
				for_each(idx2,nbf2tot, idx4,nbf4tot){
					vals[idx2*nbf4tot + idx4] = (value_t)ints4q.get_element(idx2, idx4);
				}
				o.write((char*)&vals, nbf2tot*nbf4tot*sizeof(value_t));
			}
			else {
				assert(not_implemented);
			}
		}
		if(opts.max_only){
			o.write((char*)&maxvals, nbfpairs*sizeof(value_t));
		}
	}
	o.close();

	timer.exit("trans3q4q", threadnum_);

	DBG_MSG("Compute thread shutting down.");

}


/////////////////////////////////////////////////////////////////////////////
// FullTransCommThread class

FullTransCommThread::FullTransCommThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const Ref<GaussianBasisSet>& bs1,
		const Ref<GaussianBasisSet>& bs2,
		const Ref<GaussianBasisSet>& bs3,
		const Ref<GaussianBasisSet>& bs4,
		Ref<LocalSCMatrixKit>& kit
) : SparseIntsThread(msg, thr, num, bs1, bs2, bs3, bs4, kit)
{

	bf_per_node_ = new int[msg->n()-1];
	if(msg->me() == MASTER){
		for_each(ind, msg->n()-1){
			bf_per_node_[ind] = 0;
		}
	}
	const int nsh1 = bs1->nshell();
	const int nsh3 = bs3->nshell();

	pair_assignments_ = new int[nsh1*nsh3];
	int npair = basis1_->nshell()*basis3_->nshell();
	if(msg->me() == MASTER){
		for_each(ish1,nsh1, ish3,nsh3){
			int assignment = master_pair_assignment(ish1, ish3);
			pair_assignments_[ish1*nsh3 + ish3] = assignment;
		}
		msg->bcast(pair_assignments_, npair, 0);
	}
	else
		msg->bcast(pair_assignments_, npair, 0);
}

FullTransCommThread::~FullTransCommThread(){
	delete[] bf_per_node_;
	delete[] pair_assignments_;
}

int
FullTransCommThread::master_pair_assignment(int sh3, int sh4)
{
	int sh34 = sh3 * basis4_->nshell() + sh4;
	const int nbf3 = basis3_->shell(sh3).nfunction();
	const int nbf4 = basis4_->shell(sh4).nfunction();
	// This could be done more efficiently by keeping the list sorted by bf per node
	int min_val = -1, min_node = 0;
	for_each(ind, msg->n()-1){
		int nodebf = bf_per_node_[ind];
		if(min_val < 0 || nodebf < min_val){
			min_node = ind;
			min_val = nodebf;
		}
	}
	pair_assignments_[sh34] = min_node + 1;
	bf_per_node_[min_node] += nbf3*nbf4;
	return min_node + 1;
}

void
FullTransCommThread::master_run()
{
	int done_nodes = 0;
	bool halfway_message_done = false;
	int n = msg->n();
	int sh34_sender[3];
	bool nodes_still_need_data = true;
	int done_message = 0;
	while(true){
		DBG_MSG("Send thread waiting for DoneSending message");
		msg->recvt(MessageGrp::AnySender, DoneSending, done_message);
		DBG_MSG("Send thread got DoneSending message");
		done_nodes++;
		if(done_nodes == 1)
			cout << "  First finished node message received.  Still waiting on " << n-2 << " nodes." << endl;
		else if(done_nodes > n/2 && !halfway_message_done){
			halfway_message_done = true;
			cout << "  " << setprecision(1) << fixed << float(done_nodes) / float(n) * 100.0 << "% of nodes done." << endl;
		}
		if(n-1-done_nodes == 10 && n/4 > 10)
			cout << "  Waiting for 10 more nodes..." << endl;
		else if(n-1-done_nodes == 5 && n/4 > 5)
			cout << "  Waiting for 5 more nodes..." << endl;
		else if((n-1-done_nodes < 5 && n > 5) || n-1-done_nodes < 3) {
			if(n-1-done_nodes > 1)
				cout << "  Waiting for " << n-1-done_nodes << " more nodes..." << endl;
			else if(n-1-done_nodes == 1)
				cout << "  Waiting for 1 more node..." << endl;
		}
		if(done_nodes == n-1)
			break;
	}
	DBG_MSG("Master pair assignment thread shutting down");
}


/////////////////////////////////////////////////////////////////////////////
// SendThread class

SendThread::SendThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		Ref<LocalSCMatrixKit>& kit
) : FullTransCommThread(msg, thr, thr->nthread()-2, bs1, bs2, bs3, bs4, kit),
		task_queue_(),
		needs_send_message_(-1)
{
	queue_lock_ = thr->new_lock();
	comm_lock_ = thr->new_lock();
	queue_size_ = 0;
	const int nsh3 = basis3_->nshell();
	const int nsh4 = basis4_->nshell();

	if(msg->me() == MASTER){
		// Print out pair distributions and memory requirements
		if(!opts.quiet) {
			cout << "Pair distribution by nodes:" << endl;
			for_each(ind, msg->n()-1){
				cout << setw(7) << bf_per_node_[ind];
				if((ind+1) % 10 == 0){
					cout << endl;
				}
			}
			if((msg->n()-1)%10 != 0){
				cout << endl;
			}
			cout << "(Loose lower bound of) Memory requirement per node:" << endl;
			int nperpair = basis1_->nbasis() * basis2_->nbasis();
			for_each(ind, msg->n()-1){
				int size = bf_per_node_[ind] * sizeof(double) * nperpair;
				cout << setw(12) << memory_size_string(size);
				if((ind+1) % 6 == 0){
					cout << endl;
				}
			}
			if((msg->n()-1)%10 != 0){
				cout << endl;
			}
		}
	}

}

SendThread::~SendThread()
{
}

void
SendThread::run()
{
	const int me = msg->me();
	int pairs_sent = 0;

	if(me == MASTER) {
		master_run();
	}
	else {

		// We're the comm thread on a non master node, or static task distribution is enabled.
		const int nsh3 = basis3_->nshell();
		const int nsh4 = basis4_->nshell();
		vector<MessageGrp::MessageHandle> handles;
		vector<int> sizes;

		bool tasks_done[thr->nthread()-2];
		for_each(ithr, thr->nthread()-2){
			tasks_done[ithr] = false;
		}
		while(true){

#if QUEUE_MAX_ENABLED
			// Send a message indicating that we have room for more work
			MessageGrp::MessageHandle hndl;
			handles.push_back(hndl);
			int available = int(max_queue_size - queue_size_);
			sizes.push_back(available);
			DBG_MSG("Telling compute nodes that I have " << available << " bytes available.")
			msg->nb_sendt(msg->me(), QueueHasSpace, &available, 1, hndl);
#endif

			// Receive a message indicating that we can pop some work off of the queue
			int task_code = -1;
			DBG_MSG("Waiting for a NeedsSend")
			msg->recvt(msg->me(), NeedsSend, task_code);
			DBG_MSG("Recieved a NeedsSend")
			if(task_code >= 0){
				// TODO move this part to a separate method
				if(task_code > thr->nthread()-2 || tasks_done[task_code]){
					print_lock->lock();
					cout << "Received funny NeedsSend message: " << task_code << "." << endl;
					cout << "Tasks done: ";
					for_each(ithr, thr->nthread()-2){
						cout << " " << (tasks_done[ithr] ? 1 : 0);
					}
					cout << endl;
					print_lock->unlock();
					assert(false);
				}
				DBG_MSG("Received task done message from thread " << task_code);
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
			else {
				DBG_MSG("Locking queue")
				// Lock the queue
				queue_lock_->lock();

				// Get a task
				DataSendTask t = task_queue_.front();
				task_queue_.pop();

				DBG_MSG("Unlocking queue")
				// Unlock the queue
				queue_lock_->unlock();

				//const int nbf2tot = basis2_->nbasis();
				//const int nbf4 = basis4_->shell(t.sh4).nfunction();
				int sh13 = t.ish1 * nsh3 + t.jsh3;
				int dest_node = pair_assignments_[sh13];
				assert(dest_node > 0);
				assert(dest_node < msg->n());

				int idxs[5];
				idxs[0] = t.ish1;
				idxs[1] = t.jsh3;
				idxs[2] = t.sh4;
				idxs[3] = t.ndata;
				idxs[4] = msg->me();

				msg->sendt(dest_node, IndexData, idxs, 5);
				msg->sendt(dest_node, PairData, t.data, t.ndata);

				// Update the queue size
				DBG_MSG("Locking queue")
				queue_lock_->lock();
				queue_size_-= sizeof(DataSendTask) + t.ndata*sizeof(double);
				DBG_MSG("Unlocking queue")
				queue_lock_->unlock();

				++pairs_sent;

				delete[] t.data;


			}

		}

		int me = msg->me();
		// Send a message to all of the receive threads letting them know that all of the compute threads are
		//   done and just waiting to have all of the data for the third and fourth quarter transformations.
		int done_msg[] = {-1, -1, -1, -1, msg->me()};
		for(int ind = 1; ind < msg->n(); ++ind){
			msg->sendt(ind, IndexData, done_msg, 5);
		}

		// Tell the master node we don't need any more pair assignments
		int done_message = -1;
		msg->sendt(MASTER, DoneSending, done_message);
		DBG_MSG("Send thread shutting down...");
	}
	msg->sum(pairs_sent);
	const int nsh1 = basis1_->nshell();
	const int nsh3 = basis3_->nshell();
	const int nsh4 = basis4_->nshell();
	assert_equal(pairs_sent, nsh1 * nsh3 * nsh3 * nsh4);
	assert_equal(queue_size_, 0);

}

void
SendThread::distribute_shell_pair(
	vector<RefSCMatrix> ints2q,
	int ish1, int jsh3, int sh4,
	int threadnum
)
{
	assert(ish1 >= 0);
	assert(jsh3 >= 0);
	assert(sh4 >= 0);
	// Note that this is called by the Compute thread, not the Comm thread!
	timer.enter("distribute prepare", threadnum);
	DataSendTask task;
	task.ish1 = ish1;
	task.jsh3 = jsh3;
	task.sh4 = sh4;
	const int nbf1 = basis1_->shell(ish1).nfunction();
	const int nbf3 = basis3_->shell(jsh3).nfunction();
	const int nbf2tot = basis2_->nbasis();
	const int nbf4 = basis4_->shell(sh4).nfunction();
	task.ndata = nbf1 * nbf3 * nbf2tot * nbf4;
	task.data = new double[task.ndata];
	int idata = 0;
	double maxval = 0.0;
	for_each(bf1,nbf1, bf3,nbf3){
		int ipair = bf1*nbf3 + bf3;
		for_each(bf2n,nbf2tot, bf4, nbf4){
			double val = ints2q[ipair].get_element(bf2n, bf4);
			task.data[idata++] = val;
			if(val > maxval) maxval = val;
		}
	}
	timer.exit("distribute prepare", threadnum);


	timer.enter("queue wait", threadnum);
	int space_required = int(sizeof(DataSendTask) + task.ndata*sizeof(double));
#if QUEUE_MAX_ENABLED
	assert(space_required < max_queue_size);

	// Wait until there is enough space on the queue.
	int space_amt = 0;
	DBG_MSG("Waiting for space in queue.")
	while(space_amt < space_required){
		msg->recvt(msg->me(), QueueHasSpace, space_amt);
		DBG_MSG("  Queue has " << space_amt << " bytes of space, need " << space_required)
	}
	// Pass on the remaining space to the next waiting node
	int space_left = space_required - space_amt;
	MessageGrp::MessageHandle tmphndl;
	//handles_.push_back(tmphndl);
	// TODO Think through the possibility of a race condition caused by this and SendThread both posting conflicting QueueHasSpace messages
	msg->nb_sendt(msg->me(), QueueHasSpace, &space_left, 1, tmphndl);
#endif

	// Acquire the queue lock.
	queue_lock_->lock();
	timer.exit("queue wait", threadnum);

	// Push the task onto the queue
	task_queue_.push(task);
	queue_size_ += space_required;

	// Now release the lock
	queue_lock_->unlock();

	// Post a needs send message

	// The MessageHandle won't be used again unless we need to call
	//   msg->wait() or something, so it's safe to let it be deleted.
	MessageGrp::MessageHandle hndl;
	msg->nb_sendt(msg->me(), NeedsSend, &needs_send_message_, 1, hndl);


}


/////////////////////////////////////////////////////////////////////////////
// ReceiveThread class

ReceiveThread::ReceiveThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		Ref<LocalSCMatrixKit>& kit
) : FullTransCommThread(msg, thr, thr->nthread()-1, bs1, bs2, bs3, bs4, kit),
		assigned_pairs_(),
		nbf_pairs_()
{
	pairs_lock_ = thr->new_lock();
	my_ints2q_complete_ = false;
	sorted_pairs_position_ = 0;

	if(msg->me() == MASTER){
		return;
	}

	const int nsh1 = bs1->nshell();
	const int nsh3 = bs3->nshell();
	int ipair = 0;
	for_each(sh1,nsh1, sh3,nsh3){
		if(pair_assignments_[ipair++] == msg->me()) {
			IntPair mypair(sh1,sh3);
			const int nbf1 = bs1->shell(sh1).nfunction();
			const int nbf3 = bs3->shell(sh3).nfunction();
			const int my_nbf = nbf1*nbf3;
			// Insert ordered by nbf, decending
			bool inserted = false;
			for_each(ipr, assigned_pairs_.size()){
				const int nbf1i = bs1->shell(assigned_pairs_[ipr].first).nfunction();
				const int nbf3i = bs3->shell(assigned_pairs_[ipr].second).nfunction();
				if(my_nbf >= nbf1i*nbf3i){
					assigned_pairs_.insert(assigned_pairs_.begin()+ipr, mypair);
					inserted = true;
					break;
				}
			}
			if(!inserted){
				assigned_pairs_.push_back(mypair);
			}
			for_each(bf1,nbf1, bf3,nbf3){
				RefSCMatrix pmat = kit_->matrix(bs2->basisdim(), bs4->basisdim());
				pmat.assign(0.0);
				my_ints2q[mypair].push_back(pmat);
			}
		}
	}

}

bool
ReceiveThread::get_my_next_pair(int& sh3, int& sh4){
	pairs_lock_->lock();
	if(assigned_pairs_.size() == sorted_pairs_position_){
		pairs_lock_->unlock();
		return false;
	}
	else{
		IntPair sh34 = assigned_pairs_[sorted_pairs_position_];
		sorted_pairs_position_++;
		pairs_lock_->unlock();
		sh3 = sh34.first;
		sh4 = sh34.second;
		return true;
	}
}

void
ReceiveThread::run()
{
	const int me = msg->me();
	int pair_parts_received = 0;
	if(me == MASTER){
		// This only needs to run on nodes participating in the compute process
		// If assignment of pairs becomes a bottleneck, we could have multiple
		//   pair assignment threads ("master_run()"), though I'm not sure that
		//   would necessarily solve the problem.
		DBG_MSG("Master receive thread exiting.")
	}
	else{

		bool nodes_done[msg->n()-1];
		for_each(inode, msg->n()-1){
			nodes_done[inode] = false;
		}

		int memory_used = 0;

		const int nsh1 = basis1_->nshell();
		const int nsh3 = basis3_->nshell();
		const int nsh4 = basis4_->nshell();
		int pairs_recd[nsh1*nsh3];
		for_each(i, nsh1*nsh3) pairs_recd[i] = 0;

		while(true){

			int idxs[5];
			int ish1, ibf1, jsh3, jbf3, sh4, sender, ndata;
			DBG_MSG("Waiting for IndexData")
			msg->recvt(MessageGrp::AnySender, IndexData, idxs, 5);
			DBG_MSG("Received for IndexData")
			if(idxs[0] == -1){
				nodes_done[idxs[4]-1] = true;
				bool all_nodes_done = true;
				for_each(inode, msg->n()-1){
					if(!nodes_done[inode]){
						all_nodes_done = false;
						break;
					}
				}
				DBG_MSG("Receive thread received done message from node " << idxs[4])
				if(all_nodes_done)
					break;
			}
			else {
				ish1 = idxs[0];
				jsh3 = idxs[1];
				sh4 = idxs[2];
				ndata = idxs[3];
				sender = idxs[4];

				// Receive the data
				double data[ndata];
				DBG_MSG("Waiting for PairData from sender " << sender)
				msg->recvt(sender, PairData, data, ndata);
				DBG_MSG("Received for PairData")

				const int nbf2tot = basis2_->nbasis();
				const int nbf1 = basis1_->shell(ish1).nfunction();
				const int nbf3 = basis3_->shell(jsh3).nfunction();
				const int nbf4 = basis4_->shell(sh4).nfunction();
				const int sh4off = basis4_->shell_to_function(sh4);
				double maxval = 0;
				IntPair sh13(ish1, jsh3);
				int idata = 0;
				for_each(ibf1,nbf1, jbf3,nbf3){
					int bf13pair = ibf1 * nbf3 + jbf3;
					RefSCMatrix ij_ints2q = my_ints2q[sh13][bf13pair];
					for_each(bf2n,nbf2tot, bf4,nbf4){
						int idx4 = sh4off + bf4;
						double val = data[idata++];
						ij_ints2q.accumulate_element(bf2n, idx4, val);
						if(fabs(val) > maxval) maxval = fabs(val);
					}
				}
				pairs_recd[ish1*nsh3+jsh3]++;
				if(maxval > 1e-15) DBG_MSG("  Rec'd shell pair " << ish1 << ", " << jsh3 << " from " << sender << ": " << maxval)
				assert(idata==ndata);
				++pair_parts_received;
			}
		}

		for_each(ipair, assigned_pairs_.size()){
			int sh1 = assigned_pairs_[ipair].first;
			int sh3 = assigned_pairs_[ipair].second;
			assert_equal(pairs_recd[sh1*nsh3 + sh3], nsh4*nsh3)
		}

		my_ints2q_complete_ = true;

		// Now that all of the nodes are done sending us stuff,
		//   we can tell the compute threads to get back to work.
		DBG_MSG("Receive thread shutting down, waking up compute threads");
		for_each(ithr, thr->nthread()-2){
			// This would be better done with pthreads conditionals
			msg->sendt(msg->me(), HaveAllData, me);
		}
	}

	const int nsh1 = basis1_->nshell();
	const int nsh3 = basis3_->nshell();
	const int nsh4 = basis4_->nshell();
	int num_assigned_pairs = assigned_pairs_.size();
	int total_pairs_received = nsh3 * nsh4 * num_assigned_pairs;
	assert_equal(pair_parts_received, total_pairs_received);
	//msg->sum(pair_parts_received);
	//msg->sum(num_assigned_pairs);
	//assert_equal(num_assigned_pairs, nsh1 * nsh3);
	//assert_equal(pair_parts_received, nsh1 * nsh3 * nsh4);
}
