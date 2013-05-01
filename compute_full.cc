/*
 * compute_ints.cc
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "compute_full.h"
#include "utils.h"
#include "binfiles.h"
#include "compute_thread.h"

using namespace sparse_ints;
using namespace sc;
using namespace std;

#define DBG(...) if(opts.debug) cout << "  " << "NODE" << me << "::" << __VA_ARGS__ << endl
#define DONE_COLLECTING 10

enum { DoneCollecting = 10 };

void
sparse_ints::compute_full_trans_ints(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
        const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
        sc::TwoBodyOper::type otype,
        std::string prefix, std::string matname, std::string tmpdir,
        const sc::Ref<sc::GaussianBasisSet>& bs1,
        const sc::Ref<sc::GaussianBasisSet>& bs2,
        const sc::Ref<sc::GaussianBasisSet>& bs3,
        const sc::Ref<sc::GaussianBasisSet>& bs4,
        sc::RefSymmSCMatrix& P1,
        sc::RefSymmSCMatrix& P2,
        sc::RefSymmSCMatrix& P3,
        sc::RefSymmSCMatrix& P4,
    	const sc::Ref<sc::LocalSCMatrixKit>& inkit
)
{
	/*=========================================================*/
	/* Set up local variables  	                          {{{1 */ #if fold_begin
    int me = msg->me();
    int nthr = thr->nthread();

    int* quartets_processed[msg->n()];
    for_each(iproc, msg->n()){
    	quartets_processed[iproc] = allocate<int>(nthr);
    	for_each(ithr, thr->nthread()){
    		quartets_processed[iproc][ithr] = 0;
    	}
    }

    /***********************************************************/ #endif //1}}}
    /*=========================================================*/
	/* Set up the worker threads                          {{{1 */ #if fold_begin

	// Crete the locks
    Ref<ThreadLock> compute_lock = thr->new_lock();
    Ref<ThreadLock> comm_lock = thr->new_lock();
    Ref<ThreadLock> queue_lock = thr->new_lock();

    // Create the matrix kit if needed
    Ref<LocalSCMatrixKit> kit;
    if(inkit.nonnull())
    	kit = inkit;
    else
    	kit = new LocalSCMatrixKit();

    // Make a prefix for the temporary files used by the compute nodes
    string tmp_prefix;
    tmp_prefix = tmpdir + "/node_";

    // Create the send and receive thread objects
    DBG("Creating send thread");
    SendThread* send_thread = new SendThread(msg, thr, bs1, bs2, bs3, bs4, kit);
    thr->add_thread(nthr-2, send_thread);
    DBG("Creating receive thread");
    ReceiveThread* recv_thread = new ReceiveThread(msg, thr, bs1, bs2, bs3, bs4, kit);
    thr->add_thread(nthr-1, recv_thread);

    // Create the comptue threads
    for_each(ithr, nthr-2){
    	DBG("Creating compute thread " << ithr);
    	FullTransComputeThread* comp_thread = new FullTransComputeThread(msg, thr,
    			ithr, intdescr, compute_lock,
    			bs1, bs2, bs3, bs4,
    			otype, tmp_prefix, P1, P2, P3, P4, kit,
    			send_thread, recv_thread,
    			&(quartets_processed[me][ithr])
    	);
    	thr->add_thread(ithr, comp_thread);
    }

    /***********************************************************/ #endif //1}}}
    /*=========================================================*/
	/* Run the threads     	                              {{{1 */ #if fold_begin

	timer.enter("run threads");
	DBG("Running " << nthr << " total threads, two of which are comm threads.");
    thr->start_threads();
    thr->wait_threads();
	timer.exit("run threads");

	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
	/* Gather results on each node                        {{{1 */ #if fold_begin

    timer.enter("gather results");
    timer.enter("node gather");
	if(!opts.quiet && me == MASTER) cout << "  Gathering results on a per-node basis..." << endl;
	string outfile, infile;
	{
		stringstream sstr;
		sstr << tmp_prefix << msg->me() << ".bin";
		outfile = sstr.str();
	}
	//const char* filename = sstr.str().c_str();
	if(me != MASTER) {
		ofstream o(outfile.c_str(), ios::out | ios::binary);
		// TODO Thread this part?  It seems to be taking a while...
		// Loop over the temporary files created by the threads
		ifstream i;
		for(int ithr = 0; ithr < nthr-2; ++ithr){
			{
				stringstream sstr;
				sstr << tmp_prefix << msg->me() << "_" << ithr << ".bin";
				infile = sstr.str();
			}
			if(opts.debug){
				cout << "        Node " << me << " had file " << infile << " of size "
					 << file_size_string(infile) << endl;
			}
			i.open(infile.c_str(), ios::in | ios::binary);
			copy_buffer(i, o);
			i.close();
			remove(infile.c_str());
		}
		o.close();
		msg->sendt(MASTER, DoneCollecting, me);
	}
	else{
		// Master node listens for nodes to be done.
		int ndone = 0;
		int node_done = 0;
		bool done_nodes[msg->n()-1];
		for_each(ind, msg->n()-1)
			done_nodes[ind] = false;
		while(ndone < msg->n() - 1){
			msg->recvt(MessageGrp::AnySender, DoneCollecting, node_done);
			ndone++;
			done_nodes[node_done-1] = true;
			if(msg->n() - 1 - ndone < msg->n()/10){
				cout << "    Node " << node_done << " reported finished collecting." << endl;
				if(ndone < msg->n() - 1){
					cout << "    Still waiting on " << (msg->n()-1-ndone) << " nodes.";
					int tmpline = 0;
					for_each(ind, msg->n() - 1) {
						if(!done_nodes[ind]){
							if(tmpline++ % 15 == 0){
								cout << endl << "    ";
							}
							cout << setw(4) << ind+1;
						}
					}
					cout << endl;
				}
				else
					cout << endl;
			} else if(msg->n() - 1 - ndone == 10){
				cout << ndone << " nodes reported finished collecting. Still waiting on 10 nodes." << endl;
			}
		}
	}
	// Report memory status:
	//ExEnv::outn() << indent << "Node " << msg->me() << " consumable resources summary after half transform and all sends:" << endl;
	//ExEnv::outn() << incindent;
	//ConsumableResources::get_default_instance()->print_summary(ExEnv::outn(), true, true);
	//ExEnv::outn() << decindent;

	timer.exit("node gather");
	DBG("Waiting for other nodes to finish");
	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
	/* Sync nodes, waiting for all nodes to finish        {{{1 */ #if fold_begin
	timer.enter("sync");
	msg->sync();
	timer.exit("sync");
	for_each(inode, msg->n()){
		msg->sum(quartets_processed[inode], thr->nthread());
    }
	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
	/* Gather results from all nodes into one file        {{{1 */ #if fold_begin
	// This could be done in a tree fashion if it becomes a bottleneck...
	// TODO Tar instead of gather
	timer.enter("master gather");
	if(me == MASTER) {
		if(!opts.quiet) cout << "  Gathering results from all nodes..." << endl;
		ifstream i;
		string fname = prefix + matname + ".bin";
		ofstream o(fname.c_str(), ios::out | ios::binary);
		write_header(o, bs1, bs2, bs3, bs4, opts.out_type);
		for(int inode = 1; inode < msg->n(); ++inode){
			{
				stringstream sstr;
				sstr << tmp_prefix << inode << ".bin";
				infile = sstr.str();
			}
			if(!opts.quiet)
				cout << "    Node " << inode << " had file " << infile << " of size " << file_size_string(infile) << endl;
			i.open(infile.c_str(), ios::in | ios::binary);
			copy_buffer(i, o);
			i.close();
			remove(infile.c_str());
		}
		o.close();

		// Now print out the quartet distribution
		cout << "  Quartets processed per thread:" << endl << "    ";
		int n_on_line = 0;
		for(int inode = 1; inode < msg->n(); ++inode){
			for(int ithr = 0; ithr < thr->nthread()-2; ++ithr){
				if(n_on_line % 10 == 0){
					n_on_line = 0;
					cout << endl << "    ";
				}
				cout << setw(7) << quartets_processed[inode][ithr];
				n_on_line++;
			}
		}
		cout << endl;
	}
	timer.exit("master gather");
	timer.exit("gather results");

	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
	/* Cleanup             		                          {{{1 */ #if fold_begin
    thr->delete_threads();
    for_each(iproc, msg->n()){
    	deallocate(quartets_processed[iproc]);
    }
	/***********************************************************/ #endif //1}}}
	/*=========================================================*/
}

