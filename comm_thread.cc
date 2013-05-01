/*
 * comm_thread.cc
 *
 *  Created on: May 1, 2013
 *      Author: dhollman
 */

#include "compute_thread.h"
#include "utils.h"

using namespace std;
using namespace sc;
using namespace sparse_ints;

// Doesn't work with 1
#define QUEUE_MAX_ENABLED 0

size_t SendThread::max_queue_size = size_t(1e7);

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

	bf_per_node_ = allocate<int>(msg->n()-1);
	if(msg->me() == MASTER){
		for_each(ind, msg->n()-1){
			bf_per_node_[ind] = 0;
		}
	}
	const int nsh1 = bs1->nshell();
	const int nsh3 = bs3->nshell();

	pair_assignments_ = allocate<int>(nsh1*nsh3);
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
	deallocate(bf_per_node_);
	deallocate(pair_assignments_);
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
		needs_send_message_(-1),
		queue_full_warning_printed_(false)
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
			cout << "Rough memory requirement per node:" << endl;
			int nperpair = basis1_->nbasis() * basis2_->nbasis();
			Ref<ConsumableResources> cr = ConsumableResources::get_default_instance();
			size_t maxmem = cr->max_memory();
			for_each(ind, msg->n()-1){
				size_t size = size_t(bf_per_node_[ind]) * sizeof(double) * size_t(nperpair);
				if(size + max_queue_size > maxmem){
					throw LimitExceeded<size_t>(
							"Not enough memory to do full transformation",
							__FILE__,
							__LINE__,
							maxmem,
							size + max_queue_size
					);
				}
				cout << setw(12) << memory_size_string(int(size + max_queue_size));
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

bool
SendThread::do_send_task(bool* tasks_done, int& pairs_sent){

	const int nsh3 = basis3_->nshell();
	const int nsh4 = basis4_->nshell();

	#if QUEUE_MAX_ENABLED
	// Send a message indicating that we have room for more work
	MessageGrp::MessageHandle hndl;
	handles.push_back(hndl);
	size_t available = size_t(max_queue_size - queue_size_);
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
			return true;
	}
	else {
		timer.enter("send thread: 01 queue lock", threadnum_);
		DBG_MSG("Locking queue")
		// Lock the queue
		queue_lock_->lock();
		timer.exit("send thread: 01 queue lock", threadnum_);

		// Get a task
		timer.enter("send thread: 02 get task", threadnum_);
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
		timer.exit("send thread: 02 get task", threadnum_);

		timer.enter("send thread: 03 send data", threadnum_);
		DBG_MSG("Sending IndexData " << idxs[0] << ", " << idxs[1] << ", " << idxs[2] << ", " << idxs[3] << ", " << idxs[4] << ".");
		msg->sendt(dest_node, IndexData, idxs, 5);
		DBG_MSG("done sending index data.  Sending integral data of size " << t.ndata);
		msg->sendt(dest_node, PairData, t.data, t.ndata);
		DBG_MSG("done sending integral data");
		timer.exit("send thread: 03 send data", threadnum_);


		// Update the queue size
		timer.enter("send thread: 04 queue update", threadnum_);
		deallocate(t.data);
		DBG_MSG("Locking queue")
		queue_lock_->lock();
		queue_size_-= sizeof(DataSendTask) + t.ndata*sizeof(double);
		DBG_MSG("Unlocking queue")
		queue_lock_->unlock();

		++pairs_sent;
		timer.exit("send thread: 04 queue update", threadnum_);
	}

	return false;
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
		#if QUEUE_MAX_ENABLED
		vector<MessageGrp::MessageHandle> handles;
		vector<int> sizes;
		#endif

		bool tasks_done[thr->nthread()-2];
		for_each(ithr, thr->nthread()-2){
			tasks_done[ithr] = false;
		}
		timer.enter("send thread main loop", threadnum_);
		while(true){
			if(do_send_task(tasks_done, pairs_sent)){
				break;
			}
		}
		timer.exit("send thread main loop", threadnum_);

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
	SCMatrix** ints2q,
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
	task.data = allocate<double>(task.ndata);
	int idata = 0;
	double maxval = 0.0;
	for_each(bf1,nbf1, bf3,nbf3){
		int ipair = bf1*nbf3 + bf3;
		for_each(bf2n,nbf2tot, bf4, nbf4){
			double val = ints2q[ipair]->get_element(bf2n, bf4);
			task.data[idata++] = val;
			if(val > maxval) maxval = val;
		}
		delete ints2q[ipair];
		ints2q[ipair] = 0;
	}
	assert(idata==task.ndata);
	timer.exit("distribute prepare", threadnum);


	timer.enter("queue wait", threadnum);
	size_t space_required = sizeof(DataSendTask) + size_t(task.ndata)*sizeof(double);
#if QUEUE_MAX_ENABLED
	assert(space_required < max_queue_size);

	// Wait until there is enough space on the queue.
	size_t space_amt = 0;
	DBG_MSG("Waiting for space in queue.")
	while(space_amt < space_required){
		msg->recvt(msg->me(), QueueHasSpace, &space_amt, 1);
		DBG_MSG("  Queue has " << space_amt << " bytes of space, need " << space_required)
		assert(space_amt >= 0);
	}
	// Pass on the remaining space to the next waiting node
	size_t space_left = space_required - space_amt;
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
	timer.enter("notify send thread", threadnum);
	if(queue_size_ > max_queue_size){
		// Do a blocking send to let the send thread catch up
		DBG_MSG("(thread " << threadnum << ") queue was full, queue_size_ was " << queue_size_ << ".  Using blocking send.");
		if(!queue_full_warning_printed_){
			queue_full_warning_printed_ = true;
			print_lock->lock();
			ExEnv::outn() << "Warning: Send queue full on node " << msg->me() << ".  Blocking sends will now be used to notify send thread."<< endl
					<< "The program may run slowly now; increase memory if possible." << endl;
			print_lock->unlock();
		}
		msg->sendt(msg->me(), NeedsSend, &needs_send_message_, 1);
	}
	else{
		// Post a non-blocking send and get out of here
		// The MessageHandle won't be used again unless we need to call
		//   msg->wait() or something, so it's safe to let it be deleted.
		MessageGrp::MessageHandle hndl;
		DBG_MSG("(thread " << threadnum << ") queue has space, so nb_sendt was done.  queue_size_ was " << queue_size_ << ".");
		msg->nb_sendt(msg->me(), NeedsSend, &needs_send_message_, 1, hndl);
	}
	timer.exit("notify send thread", threadnum);


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
			int nbfpairs = nbf1*nbf3;
			my_ints2q[mypair] = new SCMatrix*[nbfpairs];
			for_each(bf1,nbf1, bf3,nbf3){
				RefSCMatrix pmat = kit_->matrix(bs2->basisdim(), bs4->basisdim());
				pmat->assign(0.0);
				int bfpairnum = bf1*nbf3 + bf3;
				my_ints2q[mypair][bfpairnum] = pmat;
				pmat->reference();
			}
		}
	}
}

ReceiveThread::~ReceiveThread()
{
	// Delete the SCMatrix objects allocated in the contructor
	if(msg->me() != MASTER){
		const int nsh1 = basis1_->nshell();
		const int nsh3 = basis3_->nshell();
		int ipair = 0;
		std::map<IntPair, SCMatrix**>::iterator iter;
		for (iter = my_ints2q.begin(); iter != my_ints2q.end(); ++iter){
			const int nbf1 = basis1_->shell((*iter).first.first).nfunction();
			const int nbf3 = basis3_->shell((*iter).first.second).nfunction();
			for_each(bf1,nbf1, bf3,nbf3){
				int bfpairnum = bf1*nbf3 + bf3;
				SCMatrix* ptr = (*iter).second[bfpairnum];
				ptr->dereference();
				assert(ptr->nreference() == 0);
				delete ptr;
			}
			delete[] (*iter).second;
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
		cr->consume_memory((msg->n()-1)*sizeof(bool));
		for_each(inode, msg->n()-1){
			nodes_done[inode] = false;
		}

		int memory_used = 0;

		const int nsh1 = basis1_->nshell();
		const int nsh3 = basis3_->nshell();
		const int nsh4 = basis4_->nshell();

		int pairs_recd[nsh1*nsh3];
		cr->consume_memory(nsh1*nsh3*sizeof(int));

		for_each(i, nsh1*nsh3) pairs_recd[i] = 0;

		while(true){

			int idxs[5];
			cr->consume_memory(5*sizeof(int));

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
				if(all_nodes_done){
					// for idxs leaving memory:
					cr->release_memory(5*sizeof(int));
					break;
				}
			}
			else {
				ish1 = idxs[0];
				jsh3 = idxs[1];
				sh4 = idxs[2];
				ndata = idxs[3];
				sender = idxs[4];

				// Receive the data
				double data[ndata];
				cr->consume_memory(ndata*sizeof(double));

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
					SCMatrix* ij_ints2q = my_ints2q[sh13][bf13pair];
					for_each(bf2n,nbf2tot, bf4,nbf4){
						int idx4 = sh4off + bf4;
						double val = data[idata++];
						ij_ints2q->accumulate_element(bf2n, idx4, val);
						if(fabs(val) > maxval) maxval = fabs(val);
					}
				}
				pairs_recd[ish1*nsh3+jsh3]++;
				if(maxval > 1e-15) DBG_MSG("  Rec'd shell pair " << ish1 << ", " << jsh3 << " from " << sender << ": " << maxval)
				assert(idata==ndata);
				++pair_parts_received;

				// data goes out of scope and we recover that memory
				cr->release_memory(ndata*sizeof(double));
			}
			// idxs is going out of memory
			cr->release_memory(5*sizeof(int));
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

		// The following variables are going out of scope:
		// nodes_done:
		cr->release_memory((msg->n()-1)*sizeof(bool));
		// pairs_recd:
		cr->release_memory(nsh1*nsh3*sizeof(int));
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



