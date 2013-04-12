/*
 * compute_hti.cc
 *
 *  Created on: Apr 11, 2013
 *      Author: dhollman
 */

#include "compute_hti.h"
#include "utils.h"
#include "binfiles.h"
#include "compute_thread.h"

using namespace sparse_ints;
using namespace sc;
using namespace std;

void
sparse_ints::compute_hti_threaded(
        const Ref<TwoBodyIntDescr>& intdescr,
        TwoBodyOper::type* otypes, string* descs, int num_types,
        string prefix, string tmpdir,
        const Ref<GaussianBasisSet>& bs13,
        const Ref<GaussianBasisSet>& bs24,
    	DensityMap dens_pairs,
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
    std::map<string, vector<string> > tmp_prefixes;
    DensityMapIterator mapit = dens_pairs.begin();
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
        		dens_pairs, kit, &(quartets_processed[me*nthr + ithr])
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
	//if(!opts.quiet && me == MASTER)
	//   cout << "  Waiting for all nodes to finish aggregation..." << endl;
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
		for(mapit = dens_pairs.begin(); mapit != dens_pairs.end(); mapit++) {
			string pair_name = (*mapit).first;
			for_each(ity, num_types){
				string fname = prefix + pair_name + descs[ity] + ".bin";
				ofstream o(fname.c_str(), ios::out | ios::binary);
				// TODO update this to modern write
				write_header(o, bs13, bs24, bs13, bs24, opts.out_type);
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
			cout << "Quartet distribution:" << setw(5) << endl;
			for_each(iproc,msg->n(), ithr,nthr){
				if(n_on_line == n_per_line){
					cout << endl;
					n_on_line = 0;
				}
				cout << setw(7) << quartets_processed[iproc*nthr + ithr] << " ";
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

}

