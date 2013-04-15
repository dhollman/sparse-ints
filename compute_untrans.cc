/*
 * compute_hti.cc
 *
 *  Created on: Apr 11, 2013
 *      Author: dhollman
 */

#include "compute_untrans.h"
#include "utils.h"
#include "binfiles.h"
#include "compute_thread.h"

using namespace sparse_ints;
using namespace sc;
using namespace std;

void
sparse_ints::compute_untrans_threaded(
        const Ref<TwoBodyIntDescr>& intdescr,
        TwoBodyOper::type* otypes, string* descs, int num_types,
        string prefix, string tmpdir,
        const Ref<GaussianBasisSet>& bs1,
        const Ref<GaussianBasisSet>& bs2,
        const Ref<GaussianBasisSet>& bs3,
        const Ref<GaussianBasisSet>& bs4,
    	Ref<LocalSCMatrixKit>& kit
)
{
    //============================================================
    // Setup

	timer.enter("compute threads setup");
    int me = msg->me();
    int nthr = thr->nthread();

    // Create the quartets processed array
    int quartets_processed[msg->n()*nthr];
    for_each(ii,msg->n()*nthr){
    	quartets_processed[ii] = 0;
    }

    if(!opts.quiet && me == MASTER){
    	cout << "Computing integrals ";
    	for_each(ity, num_types-1)
    		cout << descs[ity] << ", ";
    	cout << descs[num_types-1] << "..." << endl;
    }

    // Setup the worker threads
    Ref<ThreadLock> lock = thr->new_lock();
    vector<string> tmp_prefixes;
	for_each(ity, num_types) {
		tmp_prefixes.push_back(tmpdir + "/" + descs[ity] + "_");
	}
    for_each(ithr,nthr){
    	timer.enter("create compute threads");
        UntransComputeThread* thread = new UntransComputeThread(
        		ithr, intdescr, lock,
        		bs1, bs2, bs3, bs4,
        		otypes, tmp_prefixes, num_types,
        		kit, &(quartets_processed[me*nthr + ithr])
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
	if(me != MASTER){
		for_each(ity, num_types){
			stringstream sstr;
			sstr << tmp_prefixes[ity] << msg->me() << ".bin";
			string outfile = sstr.str();
			//const char* filename = sstr.str().c_str();
			ofstream o(outfile.c_str(), ios::out | ios::binary);
			// Loop over the temporary files created by the threads
			ifstream i;
			for_each(ithr, nthr){
				stringstream sstr;
				sstr << tmp_prefixes[ity] << msg->me() << "_" << ithr << ".bin";
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
	msg->sum(quartets_processed, msg->n()*nthr);
	if(me == MASTER) {
		if(!opts.quiet) cout << "  Gathering results from all nodes..." << endl;
		ifstream i;
		for_each(ity, num_types){
			string fname = prefix + descs[ity] + ".bin";
			ofstream o(fname.c_str(), ios::out | ios::binary);
			if(opts.out_type == AllInts)
				write_header(o, bs1, bs2, bs3, bs4, AllIntsUntransformed);
			else
				write_header(o, bs1, bs2, bs3, bs4, opts.out_type);
			for_each(inode, msg->n()){
				stringstream sstr;
				sstr << tmp_prefixes[ity] << inode << ".bin";
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
			if(!opts.quiet) cout << "    Collected complete file for " << descs[ity] << " of size " << file_size_string(fname.c_str()) << endl;
		}
		// Print the number of shell quartets computed per node
		if(!opts.quiet){
			cout << "Quartet distribution:" << setw(5) << endl;
			PRINT_LIST(quartets_processed, msg->n()*nthr, 7, 10);
			cout << endl;
		}
	}
	timer.exit("master gather");
	timer.exit("gather results");


    //============================================================
    // Cleanup
    thr->delete_threads();

}

