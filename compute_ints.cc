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

void
sparse_ints::compute_half_trans_ints(
        const Ref<TwoBodyIntDescr>& intdescr,
        TwoBodyOper::type* otypes, string* descs, int num_types,
        string prefix, string tmpdir,
        const Ref<GaussianBasisSet>& bs13,
        const Ref<GaussianBasisSet>& bs24,
    	DensityMap dens_pairs,
    	const Ref<SCMatrixKit>& inkit
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
    	const sc::Ref<sc::SCMatrixKit>& inkit
)
{
    int me = msg->me();
    int nthr = thr->nthread();

    // Setup the worker threads
    Ref<ThreadLock> compute_lock = thr->new_lock();
    Ref<ThreadLock> comm_lock = thr->new_lock();
    Ref<ThreadLock> queue_lock = thr->new_lock();
    Ref<SCMatrixKit> kit;
    if(inkit.nonnull())
    	kit = inkit;
    else
    	kit = P1->kit();
    FullTransCommThread* send_thread = new FullTransCommThread(
    		comm_lock, queue_lock,
    		bs1, bs2, bs3, bs4, FullTransCommThread::SendThread, kit
    );
    thr->add_thread(0, send_thread);
    FullTransCommThread* recv_thread = new FullTransCommThread(
    		comm_lock, queue_lock,
    		bs1, bs2, bs3, bs4, FullTransCommThread::ReceiveThread, kit
    );
    thr->add_thread(1, recv_thread);
    for(int ithr = 2; ithr < nthr; ++ithr){
    	FullTransComputeThread* comp_thread = new FullTransComputeThread(
    			ithr, intdescr, compute_lock,
    			bs1, bs2, bs3, bs4,
    			otype, prefix, P1, P2, P3, P4, kit,
    			send_thread, recv_thread
    	);
    	thr->add_thread(1, comp_thread);
    }

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
	if(!opts.quiet && me == MASTER) cout << "  Gathering results on a per-node basis..." << endl;
	stringstream sstr;
	sstr << prefix << msg->me() << ".bin";
	string outfile = sstr.str();
	//const char* filename = sstr.str().c_str();
	ofstream o(outfile.c_str(), ios::out | ios::binary);
	// Loop over the temporary files created by the threads
	ifstream i;
	for(int ithr = 2; ithr < nthr; ++ithr){
		stringstream sstr;
		sstr << prefix << msg->me() << "_" << ithr << ".bin";
		string infile = sstr.str();
		if(opts.debug){
			cout << "    Node " << me << " had file " << infile
				 << " for " << prefix << " of size "
				 << file_size_string(infile) << endl;
		}
		i.open(infile.c_str(), ios::in | ios::binary);
		copy_buffer(i, o);
		i.close();
		remove(infile.c_str());
	}
	o.close();
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
		string fname = prefix + ".bin";
		ofstream o(fname.c_str(), ios::out | ios::binary);
		// TODO !!! change this when using multiple basis sets!!!
		write_header(o, bs1, opts);
		for_each(inode, msg->n()){
			stringstream sstr;
			sstr << prefix << inode << ".bin";
			string infile = sstr.str();
			if(opts.debug)
				cout << "    Node " << inode << " had file " << infile << " for " << prefix << " of size " << file_size_string(infile) << endl;
			i.open(infile.c_str(), ios::in | ios::binary);
			copy_buffer(i, o);
			i.close();
			remove(sstr.str().c_str());
		}
		o.close();
	}
	timer.exit("master gather");
	timer.exit("gather results");


    //============================================================
    // Cleanup
    thr->delete_threads();

}

