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

size_t SendThread::max_queue_size = size_t(5e7);

// macros
#define DBG_MSG(mymsg) \
	if(opts.debug){\
		/*print_lock->lock();*/ \
		dbg_out_ << mymsg << endl;\
		dbg_out_.flush(); \
		/*print_lock->unlock();*/ \
	}


/////////////////////////////////////////////////////////////////////////////
// ComputeThread class

ComputeThread::ComputeThread(
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::Ref<sc::LocalSCMatrixKit>& kit
) : SparseIntsThread(num, bs1, bs2, bs3, bs4, kit)
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
) : ComputeThread(num, intdescr, lock, bs1, bs2, bs3, bs4, kit),
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
	int sh1 = 0, sh2 = 0;
	const int nsh3 = basis3_->nshell();
	const int nsh4 = basis4_->nshell();
	const int nbf3tot = basis3_->nbasis();
	const int nbf4tot = basis4_->nbasis();

	idx_t identifier[4];

	// Compute stuff until there's no more work left
	timer.enter("trans1q2q", threadnum_);
	while(shellpairs.get_task(sh1, sh2)){
		DBG_MSG("Got task compute pair (" << sh1 << ", " << sh2 << ")");
		const int nbf1 = basis1_->shell(sh1).nfunction();
		const int nbf2 = basis2_->shell(sh2).nfunction();
		int nbfpairs = nbf1*nbf2;
		const double* buffer = inteval_->buffer(otype_);

		// Initialize the matrices that will hold the 1Q transformed integrals
		vector<RefSCMatrix> ints1q;
		for_each(ipair, nbfpairs){
			ints1q.push_back(RefSCMatrix(bsdim3_, bsdim4_, kit_));
			ints1q[ipair].assign(0.0);
		}
		// Do the transformation of index 3.
		//   Loop over all shells for indexes 3 and 4
		for_each(sh3,nsh3, sh4,nsh4){
			// Compute the shell
			timer.enter("compute shell", threadnum_);
			(*quartets_computed_)++;
			inteval_->compute_shell(sh1, sh2, sh3, sh4);
			timer.exit("compute shell", threadnum_);

			// Get the number of functions in the respective shells
			const int nbf3 = basis3_->shell(sh3).nfunction();
			const int nbf4 = basis4_->shell(sh4).nfunction();

			// Get the relevant block of the density matrix for the first quarter transform
			int sh3begin = basis3_->shell_to_function(sh3);
			int sh4begin = basis4_->shell_to_function(sh4);

			timer.enter("1q", threadnum_);
			RefSCMatrix P3part = P3_.get_subblock(0, nbf3tot - 1, sh3begin, sh3begin + nbf3 - 1);

			// Create the matrix to hold the untransformed integrals temporarily
			RefSCDimension bfdim3(new SCDimension(nbf3, 1));
			RefSCDimension bfdim4(new SCDimension(nbf4, 1));
			RefSCMatrix g12 = kit_->matrix(bfdim3, bfdim4);

			int bf1234 = 0;
			for_each(bf1,nbf1, bf2,nbf2){
				//for_each(bf3,nbf3, bf4,nbf4){
				//	g12.set_element(bf3, bf4, buffer[bf1234]);
				//	bf1234++;
				//}
				// Better:
				g12->assign(&(buffer[bf1*nbf2*nbf3*nbf4 + bf2*nbf3*nbf4]));

				// accumulate the first quarter transform
				int pairnum = bf1*nbf2 + bf2;
				//DBG_MSG("Accumulating product of P3part (" << P3part.rowdim().n() << "x" << P3part.coldim().n() << ") and g12 ("
				//			<< g12.rowdim().n() << "x" << g12.coldim().n() << ") into subblock [" << 0 << "," <<  nbf3tot - 1 << "] x ["
				//			<< sh4begin << "," << sh4begin + nbf4 - 1 << "]");
				ints1q[pairnum].accumulate_subblock(P3part * g12, 0, nbf3tot - 1, sh4begin, sh4begin + nbf4 - 1);
				//DBG_MSG("Accumulation completed.")
			}
			timer.exit("1q", threadnum_);
		}

		timer.enter("2q", threadnum_);
		// Initialize the matrices that will hold the 2Q transformed integrals
		vector<RefSCMatrix> ints2q;
		for_each(ipair, nbfpairs){
			ints2q.push_back(kit_->matrix(bsdim3_, bsdim4_));
			//ints1q[ipair].coldim().print(cout);
			//P4_.dim().print(cout);
			//DBG_MSG("Accumulating product of ints1q[ipair] (" << ints1q[ipair].rowdim().n() << "x" << ints1q[ipair].coldim().n() << ") and P4_ ("
			//		<< P4_.dim().n() << "x" << P4_.dim().n() << ") into matrix " << ints2q[ipair].rowdim().n() << "x"<< ints2q[ipair].coldim().n());
			ints2q[ipair].assign(ints1q[ipair] * P4_);
		}
		timer.exit("2q", threadnum_);

		// Now tell the comm thread to send the data to the node assigned to handle it.
		DBG_MSG("About to distribute pair " << sh1 << ", " << sh2 << ".  See send thread log.");
		timer.enter("distribute shell pair", threadnum_);
		send_thread_->distribute_shell_pair(ints2q, sh1, sh2, nbf1, nbf2, nbf3tot, nbf4tot, threadnum_);
		timer.exit("distribute shell pair", threadnum_);
		DBG_MSG("Done distributing pair" << sh1 << ", " << sh2 << ".");

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
	//MessageGrp::MessageHandle hndl;
	//msg->nb_sendt(msg->me(), NeedsSend, threadnum_, hndl);
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
		DBG_MSG("Got local pair (" << sh3 << "," << sh4 << ") to work on");
		const int nbf3 = basis3_->shell(sh3).nfunction();
		const int nbf4 = basis4_->shell(sh4).nfunction();
		IntPair sh34(sh3, sh4);
		vector<RefSCMatrix> ints34 = recv_thread_->my_ints2q[sh34];
		// TODO write only unique pairs
		int nbfpairs = nbf3*nbf4;
		idx_t identifier[4];
		value_t maxvals[nbfpairs];
		//DBG_MSG("Working on shell pair " << sh3 << "," << sh4 << ".  Vector size is " << ints34.size() << ", should be " << nbfpairs);
		for_each(bf3,nbf3, bf4,nbf4){
			int ipair = bf3*nbf4 + bf4;
			//DBG_MSG("Transforming pair "<<sh3<<"."<<bf3<<","<<sh4<<"."<<bf4
			//		<<": ints34[ipair] is a "<<ints34[ipair]->rowdim()->n()
			//		<<"x"<<ints34[ipair]->coldim()->n()<<" matrix");
			RefSCMatrix ints3q = P1_ * ints34[ipair];
			RefSCMatrix ints4q = ints3q * P2_;
			if(opts.max_only){
				//value_t maxval = max_abs<value_t, idx_t>(ints4q, idx_t[2], idx_t[3]);
				maxvals[ipair] = ints4q->maxabs();
			}
			else{
				assert(not_implemented);
			}
		}
		identifier[0] = sh3;
		identifier[2] = sh4;
		o.write((char*)&identifier, 4*sizeof(idx_t));
		o.write((char*)&maxvals, nbfpairs*sizeof(value_t));
	}
	o.flush();
	o.close();

	timer.exit("trans3q4q", threadnum_);

	DBG_MSG("Compute thread shutting down.");

}


/////////////////////////////////////////////////////////////////////////////
// FullTransCommThread class

FullTransCommThread::FullTransCommThread(
		int num,
		const Ref<GaussianBasisSet>& bs1,
		const Ref<GaussianBasisSet>& bs2,
		const Ref<GaussianBasisSet>& bs3,
		const Ref<GaussianBasisSet>& bs4,
		Ref<LocalSCMatrixKit>& kit
) : SparseIntsThread(num, bs1, bs2, bs3, bs4, kit),
		pair_assignments_()
{

	bf_per_node_ = new int[msg->n()-1];
	if(msg->me() == MASTER){
		for_each(ind, msg->n()-1){
			bf_per_node_[ind] = 0;
		}
	}

}

FullTransCommThread::~FullTransCommThread(){
	delete[] bf_per_node_;
}

int
FullTransCommThread::master_pair_assignment(int sh3, int sh4)
{
	IntPair sh34(sh3, sh4);
	const int nbf3 = basis3_->shell(sh3).nfunction();
	const int nbf4 = basis4_->shell(sh4).nfunction();
	std::map<IntPair, int>::iterator item = pair_assignments_.find(sh34);
	// This could be done more efficiently by keeping the list sorted by bf per node
	if(item == pair_assignments_.end()){
		int min_val = -1, min_node = 0;
		for_each(ind, msg->n()-1){
			int nodebf = bf_per_node_[ind];
			if(min_val < 0 || nodebf < min_val){
				min_node = ind;
				min_val = nodebf;
			}
		}
		pair_assignments_[sh34] = min_node;
		bf_per_node_[min_node] += nbf3*nbf4;
		return min_node + 1;
	}
	else{
		return pair_assignments_[sh34];
	}
}

void
FullTransCommThread::master_run()
{
	// We're the manager for the comm threads; get to work assigning places to put pairs
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
		cout << "Memory requirement per node:" << endl;
		int nperpair = basis1_->nbasis() * basis2_->nbasis();
		for_each(ind, msg->n()-1){
			int size = bf_per_node_[ind] * sizeof(double) * nperpair;
			cout << setw(9) << memory_size_string(size);
			if((ind+1) % 10 == 0){
				cout << endl;
			}
		}
		if((msg->n()-1)%10 != 0){
			cout << endl;
		}
	}
	int done_nodes = 0;
	bool halfway_message_done = false;
	int n = msg->n();
	int sh34_sender[3];
	bool nodes_still_need_data = true;
	while(true){
		DBG_MSG("Send thread waiting for NeedPairAssignment");
		msg->recvt(MessageGrp::AnySender, NeedPairAssignment, sh34_sender, 3);
		DBG_MSG("Send thread got NeedPairAssignment from sender " << sh34_sender[2]);
		if(sh34_sender[0] == -1){
			DBG_MSG("Send thread got shutdown message");
			done_nodes++;
			if(done_nodes == 1){
				cout << "  First finished node message received.  Still waiting on " << n-2 << " nodes." << endl;
			}
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
			else
				continue;
		}
		else{
			int node_assignment = master_pair_assignment(sh34_sender[0], sh34_sender[1]);
			DBG_MSG("Master sending PairAssignment to " << sh34_sender[2])
			msg->sendt(sh34_sender[2], PairAssignment, node_assignment);
			DBG_MSG("PairAssignment sent")
		}
	}
	if(!opts.quiet) {
		// TODO Shell distribution
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
	}
	DBG_MSG("Master pair assignment thread shutting down");
}


/////////////////////////////////////////////////////////////////////////////
// SendThread class

SendThread::SendThread(
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		Ref<LocalSCMatrixKit>& kit
) : FullTransCommThread(thr->nthread()-2, bs1, bs2, bs3, bs4, kit),
		task_queue_(),
		needs_send_message_(-1)
{
	queue_lock_ = thr->new_lock();
	comm_lock_ = thr->new_lock();
	queue_size_ = 0;
	const int nsh3 = basis3_->nshell();
	const int nsh4 = basis4_->nshell();
	cached_pair_assignments_ = new int[nsh3*nsh4];
	for_each(ish3,nsh3, ish4,nsh4){
		cached_pair_assignments_[ish3*nsh4 + ish4] = master_pair_assignment(ish3, ish4);
	}

}

SendThread::~SendThread()
{
	delete[] cached_pair_assignments_;
}

void
SendThread::run()
{
	const int me = msg->me();

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

			// TODO this would be done better using e.g. pthread conditionals
			// Receive a message alerting us that there is communication work to do

			MessageGrp::MessageHandle hndl;
			handles.push_back(hndl);
			int available = int(max_queue_size - queue_size_);
			sizes.push_back(available);
			msg->nb_sendt(msg->me(), QueueHasSpace, &available, 1, hndl);
			// TODO figure out how to go through and clean up unused entries in the vectors?

			int task_code = -1;
			DBG_MSG("Send thread waiting for NeedsSend");
			msg->recvt(msg->me(), NeedsSend, task_code);
			DBG_MSG("Send thread received NeedsSend");
			if(task_code >= 0){
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

			// Get the task
			DBG_MSG("Trying to lock queue_lock_")
			queue_lock_->lock();
			DBG_MSG("Successfully locked queue_lock_")
			DataSendTask t = task_queue_.front();
			task_queue_.pop();
			DBG_MSG("Unlocking queue_lock_")
			queue_lock_->unlock();

			// Send the data
			// Put all of the data belonging to a given index 3, 4 pair together
			for_each(sh3,nsh3, sh4,nsh4){
				int sh34 = sh3*nsh4 + sh4;
				int dest_node = cached_pair_assignments_[sh34];
				assert(dest_node < msg->n() && dest_node > 0);

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
				DBG_MSG("Send thread about to send shell quartet (" << idxs[0] << "," << idxs[1] << "," << idxs[2] << "," << idxs[3] << ") with " << pair_ndata << " entries.");
				for_each(bf1,t.nbf1, bf2,t.nbf2){
					int pair12 = bf1*t.nbf2 + bf2;
					data_ptr = &t.data[pair12];
					for_each(bf3,nbf3, bf4,nbf4){
						//DBG_MSG("  Putting data_ptr[" << (bf3+sh3off)*t.nbftot4 + sh4off + bf4 << "] in send_data["
						//		<< ibf << "] (pair_ndata = " << pair_ndata << "), val = " << data_ptr[(bf3+sh3off)*t.nbftot4 + sh4off + bf4]);
						send_data[ibf] = data_ptr[(bf3+sh3off)*t.nbftot4 + sh4off + bf4];
						ibf++;
					}
				}
				msg->sendt(dest_node, IndexData, idxs, 5);
				msg->sendt(dest_node, PairData, send_data, pair_ndata);
				DBG_MSG("Send index data successfully.");
				queue_lock_->lock();
				queue_size_-= sizeof(DataSendTask) + pair_ndata*sizeof(double);
				queue_lock_->unlock();
			}
			delete[] t.data;
		}

		int me = msg->me();
		// Send a message to all of the receive threads letting them know that all of the compute threads are
		//   done and just waiting to have all of the data for the third and fourth quarter transformations.
		int done_msg[] = {-1, -1, -1, -1, msg->me()};
		for(int ind = 1; ind < msg->n(); ++ind){
			msg->sendt(ind, IndexData, done_msg, 5);
		}

		// Tell the master node we don't need any more pair assignments
		int pair_assignment_end[] = {-1, -1, me};
		msg->sendt(MASTER, NeedPairAssignment, pair_assignment_end, 3);
		DBG_MSG("Send thread shutting down...");
	}

}

void
SendThread::distribute_shell_pair(
	std::vector<RefSCMatrix> pair_mats,
	int sh1, int sh2, int nbf1, int nbf2, int nbf3tot, int nbf4tot,
	int threadnum
)
{
	// Note that this is called by the Compute thread, not the Comm thread!

	timer.enter("rearrange data", threadnum);
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
	int idata = 0;
	for_each(bf1,nbf1, bf2,nbf2){
		int ipair = bf1*nbf2 + bf2;
		for_each(bf3,nbf3tot, bf4,nbf4tot){
			task.data[idata++] = pair_mats[ipair].get_element(bf3, bf4);
		}
		// This doesn't work!
		//pair_mats[ipair].convert(task.data);
		//task.data += stride;
	}
	timer.exit("rearrange data", threadnum);

	size_t space_required = sizeof(DataSendTask) + ndata*sizeof(double);
	assert(space_required < max_queue_size);
	// TODO Implement a queue size limit using e.g. pthread contitionals
	// Acquire the queue lock.
	int space_amt = 0;
	while(space_amt < space_required)
		msg->recvt(msg->me(), QueueHasSpace, space_amt);
	DBG_MSG("(outside compute thread) " << sh1 << "," << sh2 << ": trying to lock queue_lock_")
	timer.enter("push task", threadnum);
	queue_lock_->lock();
	DBG_MSG("(outside compute thread) " << sh1 << "," << sh2 << ": successfully locked queue_lock_")

	// Push the task onto the queue
	task_queue_.push(task);
	queue_size_ += space_required;

	// Now release the lock
	queue_lock_->unlock();
	timer.exit("push task", threadnum);
	DBG_MSG("(outside compute thread) " << sh1 << "," << sh2 << ": unlocked queue_lock_")

	// Post a needs send message
	DBG_MSG("(outside compute thread) " << sh1 << "," << sh2 << ": posting NeedsSend of message " << needs_send_message_)
	timer.enter("posting needs send", threadnum);
	MessageGrp::MessageHandle hndl;
	handles_.push_back(hndl);
	//messages_.push_back(-1);
	msg->nb_sendt(msg->me(), NeedsSend, &needs_send_message_, 1, hndl);
	//msg->sendt(msg->me(), NeedsSend, needs_send_message_);
	timer.exit("posting needs send", threadnum);
	DBG_MSG("(outside compute thread) " << sh1 << "," << sh2 << ": sent NeedsSend successfully.")

}

/////////////////////////////////////////////////////////////////////////////
// ReceiveThread class

ReceiveThread::ReceiveThread(
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		Ref<LocalSCMatrixKit>& kit
) : FullTransCommThread(thr->nthread()-1, bs1, bs2, bs3, bs4, kit),
		pairs_sorted_()
{
	pairs_lock_ = thr->new_lock();
	my_ints2q_complete_ = false;
	sorted_pairs_position_ = 0;
}

bool
ReceiveThread::get_my_next_pair(int& sh3, int& sh4){
	pairs_lock_->lock();
	if(pairs_sorted_.size() == sorted_pairs_position_){
		pairs_lock_->unlock();
		return false;
	}
	else{
		IntPair sh34 = pairs_sorted_[sorted_pairs_position_];
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
	if(me == MASTER){
		// This only needs to run on nodes participating in the compute process
		// If assignment of pairs becomes a bottleneck, we could have multiple
		//   pair assignment threads ("master_run()"), though I'm not sure that
		//   would necessarily solve the problem.
		DBG_MSG("Master receive thread exiting.")
		return;
	}

	bool nodes_done[msg->n()-1];
	for_each(inode, msg->n()-1){
		nodes_done[inode] = false;
	}

	while(true){

		int idxs[5];
		DBG_MSG("Receive thread waiting for IndexData");
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
			DBG_MSG("Receive thread received done message from node " << idxs[4])
			if(all_nodes_done)
				break;
			else
				continue;
		}
		DBG_MSG("Receive thread received IndexData: (" << idxs[0] << "," << idxs[1] << "," << idxs[2] << "," << idxs[3] << ")");

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
				my_ints2q[sh34].push_back(tmp);
			}
		}
		int ibf = 0;
		for_each(bf1,nbf1, bf2,nbf2, bf3,nbf3, bf4,nbf4){
			int ipair = bf3*nbf4 + bf4;
			my_ints2q[sh34][ipair].accumulate_element(sh1off+bf1, sh2off+bf2, recv_data[ibf]);
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
	/*
	if(opts.debug){
		cout << "Pairs assigned to node " << me << ":" << endl;
		for_each(ipair, pairs_sorted_.size()){
			IntPair pr = pairs_sorted_[ipair];
			cout << "  (" << setw(3) << pr.first << ", " << setw(3) << pr.second << ")";
			if(ipair > 0 && (ipair+1) % 5 == 0)
				cout << endl;
		}
		if(pairs_sorted_.size() % 5 != 0){
			cout << endl;
		}
	}
	*/

	// Now that all of the nodes are done sending us stuff,
	//   we can tell the compute threads to get back to work.
	DBG_MSG("Receive thread shutting down, waking up compute threads");
	for_each(ithr, thr->nthread()-2){
		// This would be better done with pthreads conditionals
		msg->sendt(msg->me(), HaveAllData, me);
	}
}
