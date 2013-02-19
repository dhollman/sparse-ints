/*
 * compute_ints.cc
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "compute_ints.h"
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
sparse_ints::compute_half_trans_ints(
        const Ref<TwoBodyIntDescr>& intdescr,
        TwoBodyOper::type* otypes, string* descs, int num_types,
        string prefix, string tmpdir,
        const Ref<GaussianBasisSet>& bs13,
        const Ref<GaussianBasisSet>& bs24,
    	DensityMap dens_pairs,
    	const Ref<LocalSCMatrixKit>& inkit
)
{
    //============================================================
    // Setup

	//SparseIntOptions opts = Globals::opts;

	timer.enter("compute threads setup");
    int me = msg->me();
    int nthr = thr->nthread();
    int* quartets_processed[msg->n()];
    for_each(iproc, msg->n()){
    	quartets_processed[iproc] = new int[nthr];
    }

    if(!opts.quiet && me == MASTER){
    	cout << "Computing integrals ";
    	for_each(ity, num_types-1)
    		cout << descs[ity] << ", ";
    	cout << descs[num_types-1] << "..." << endl;
    }

    // Setup the worker threads
    Ref<ThreadLock> lock = thr->new_lock();
    std::map<string, vector<string> > tmp_prefixes;
    DensityMapIterator mapit = dens_pairs.begin();
    Ref<LocalSCMatrixKit> kit;
    if(inkit.nonnull())
    	kit = inkit;
    else
    	kit = new LocalSCMatrixKit();
    for(; mapit != dens_pairs.end(); mapit++) {
    	string pair_name = (*mapit).first;
		for_each(ity, num_types) {
			tmp_prefixes[pair_name].push_back(tmpdir + "/" + pair_name + descs[ity] + "_");
		}
    }
    for_each(ithr,nthr){
    	timer.enter("create compute threads");
        HalfTransComputeThread* thread = new HalfTransComputeThread(
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
	if(!opts.quiet && me == MASTER)
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
				if(opts.debug){
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
	msg->sync();
	timer.exit("sync");
    //============================================================
	// Now gather across multiple nodes.  This could be done in a
	//   tree fashion if it becomes a bottleneck...
	timer.enter("master gather");
	if(me == MASTER) {
		if(!opts.quiet) cout << "  Gathering results from all nodes..." << endl;
		ifstream i;
		for(mapit = dens_pairs.begin(); mapit != dens_pairs.end(); mapit++) {
			string pair_name = (*mapit).first;
			for_each(ity, num_types){
				string fname = prefix + pair_name + descs[ity] + ".bin";
				ofstream o(fname.c_str(), ios::out | ios::binary);
				write_header(o, bs13, opts);
				for_each(inode, msg->n()){
					stringstream sstr;
					sstr << tmp_prefixes[pair_name][ity] << inode << ".bin";
					string infile = sstr.str();
					if(opts.debug){
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
				if(!opts.quiet) cout << "    Collected complete file for " << pair_name << descs[ity] << " of size " << file_size_string(fname.c_str()) << endl;
			}
		}
		// Print the shell pair distribution
		if(!opts.quiet){
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

void
sparse_ints::compute_full_trans_ints(
        const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
        sc::TwoBodyOper::type otype,
        std::string prefix, std::string tmpdir,
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
    int me = msg->me();
    int nthr = thr->nthread();

    int* quartets_processed[msg->n()];
    for_each(iproc, msg->n()){
    	quartets_processed[iproc] = new int[nthr];
    	for_each(ithr, thr->nthread()){
    		quartets_processed[iproc][ithr] = 0;
    	}
    }

    // Setup the worker threads
    Ref<ThreadLock> compute_lock = thr->new_lock();
    Ref<ThreadLock> comm_lock = thr->new_lock();
    Ref<ThreadLock> queue_lock = thr->new_lock();
    Ref<LocalSCMatrixKit> kit;
    if(inkit.nonnull())
    	kit = inkit;
    else
    	kit = new LocalSCMatrixKit();
    DBG("Creating send thread");
    SendThread* send_thread = new SendThread(bs1, bs2, bs3, bs4, kit);
    thr->add_thread(nthr-2, send_thread);
    DBG("Creating receive thread");
    ReceiveThread* recv_thread = new ReceiveThread(bs1, bs2, bs3, bs4, kit);
    thr->add_thread(nthr-1, recv_thread);
    string tmp_prefix;
    {
    	stringstream sstr;
    	sstr << tmpdir << "/tmp_";
    	tmp_prefix = sstr.str();
    }
    for_each(ithr, nthr-2){
    	DBG("Creating compute thread " << ithr);
    	FullTransComputeThread* comp_thread = new FullTransComputeThread(
    			ithr, intdescr, compute_lock,
    			bs1, bs2, bs3, bs4,
    			otype, tmp_prefix, P1, P2, P3, P4, kit,
    			send_thread, recv_thread,
    			&(quartets_processed[me][ithr])
    	);
    	thr->add_thread(ithr, comp_thread);
    }

    //============================================================
    // Run the threads
	timer.enter("run threads");
	DBG("Running " << nthr << " total threads, two of which are comm threads.");
    thr->start_threads();
    thr->wait_threads();
	timer.exit("run threads");

    //============================================================
    // Gather the results
    // First gather on each node...
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
	timer.exit("node gather");
	DBG("Waiting for other nodes to finish");
    //============================================================
	timer.enter("sync");
	msg->sync();
	timer.exit("sync");
	for_each(inode, msg->n()){
		msg->sum(quartets_processed[inode], thr->nthread());
    }
    //============================================================
	// Now gather across multiple nodes.  This could be done in a
	//   tree fashion if it becomes a bottleneck...
	timer.enter("master gather");
	if(me == MASTER) {
		if(!opts.quiet) cout << "  Gathering results from all nodes..." << endl;
		ifstream i;
		string fname = prefix + ".bin";
		ofstream o(fname.c_str(), ios::out | ios::binary);
		//write_header(o, bs1, bs2, bs3, bs4, MaxAbs);
		write_header(o, bs1, opts);
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


    //============================================================
    // Cleanup
    thr->delete_threads();

}

