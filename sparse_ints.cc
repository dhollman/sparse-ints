//
// sparse_ints.cc
//
// Copyright (C) 2013 David Hollman
//
// Author: David S. Hollman <dhollman@uga.edu>
// Maintainer: DSH
//
// License
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
////////////////////////////////////////////////////////////////////////////////////

#include "sparse_ints.h"
#include "utils.h"
#include "compute_ints.h"


using namespace sparse_ints;
using namespace sc;
using namespace std;

// Global variables
Ref<MessageGrp> sparse_ints::msg;
Ref<ThreadGrp> sparse_ints::thr;
MultiTimer sparse_ints::timer;
SparseIntOptions sparse_ints::opts;
Ref<ThreadLock> sparse_ints::print_lock;


int
main(int argc, char** argv) {

    //#############################################################################//
    // Setup
    //=========================================================//
    // Read the input file

    
    char *infile = new char[256];
    if (argc == 2) {
        infile = argv[1];
    } else {
        sprintf(infile, "input.dat");
    }

    Ref<KeyVal> pkv(new ParsedKeyVal(infile));
    Ref<KeyVal> keyval(new PrefixKeyVal(pkv, ":all"));

    //=========================================================//
    // Create the message group

    Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc,argv);
    if (msg.null()) {
		#if HAVE_MPI_H
        msg = new MTMPIMessageGrp();
		#else
        msg = new ProcMessageGrp();
		#endif
    }
    MessageGrp::set_default_messagegrp(msg);
    sparse_ints::msg = msg;
    int n = msg->n();
    int me = msg->me();

    //=========================================================//
    // Create the thread group object

    sparse_ints::thr = ThreadGrp::initial_threadgrp(argc, argv);
    if (thr.null()) {
    	thr = ThreadGrp::get_default_threadgrp();
    }
    ThreadGrp::set_default_threadgrp(thr);
    int nthr = thr->nthread();
    if (!opts.quiet && me == MASTER) {
        cout << "Running on " << n << " node"
             << (n == 1 ? "" : "s") << " with " << nthr << " thread"
             << (nthr == 1 ? "" : "s") << " per node." << endl;
    }
    print_lock = thr->new_lock();

    //=========================================================//
    // Create the timer
    
    MultiTimer timer("Sparse Ints", msg, thr);
    sparse_ints::timer = timer;

    //=========================================================//
    // Get the basis and the molecule from the input

    timer.enter("setup");
    timer.enter("input");
    Ref<GaussianBasisSet> obs = require_dynamic_cast<GaussianBasisSet*>(
            keyval->describedclassvalue("basis").pointer(), "main\n");
    Ref<GaussianBasisSet> auxbs = require_dynamic_cast<GaussianBasisSet*>(
            keyval->describedclassvalue("aux_basis").pointer(), "main\n");
    Ref<Molecule> mol = obs->molecule();
    SparseIntOptions opts;
    opts.debug = keyval->booleanvalue("debug", KeyValValueboolean(false));
    opts.verbose = keyval->booleanvalue("verbose", KeyValValueboolean(false));
    opts.quiet = keyval->booleanvalue("quiet", KeyValValueboolean(false));
    opts.dynamic = keyval->booleanvalue("dynamic", KeyValValueboolean(true));
    opts.max_only = keyval->booleanvalue("max_only", KeyValValueboolean(true));
    sparse_ints::opts = opts;

    bool do_eri = keyval->booleanvalue("do_eri", KeyValValueboolean(true));
    bool do_f12 = keyval->booleanvalue("do_f12", KeyValValueboolean(true));
    bool do_f12g12 = keyval->booleanvalue("do_f12g12", KeyValValueboolean(true));
    bool do_f12sq = keyval->booleanvalue("do_f12sq", KeyValValueboolean(true));
    bool do_dblcomm = keyval->booleanvalue("do_dblcomm", KeyValValueboolean(true));
    bool do_ri = keyval->booleanvalue("do_ri", KeyValValueboolean(true));

    // Figure out where to put scratch files
    char* scratch_dir = getenv("SCRATCH");
    if(scratch_dir == 0){
    	scratch_dir = new char[32];
        sprintf(scratch_dir, "/tmp");
    }
    char* job_id = getenv("PBS_JOBID");
    int pid;
    if(job_id == 0 && me == MASTER){
    	// get unique process id.  (Only if we can't get a JOBID)
    	pid_t mypid = getpid();
    	pid = (int)mypid;
    }
    // broadcast the process id
    msg->bcast(&pid, sizeof(pid_t), 0);
	{
    	stringstream sstr;
		sstr << pid;
		job_id = const_cast<char*>(sstr.str().c_str());
	}
	// Get the temporary directory as a string
    stringstream sstr;
    sstr << scratch_dir << "/" << job_id;
    string tmpdir = sstr.str();
    if(!opts.quiet){
    	if((!opts.debug && me == MASTER) || opts.debug){
    		cout << "Using temporary directory " << tmpdir << endl;
    		if(opts.dynamic) cout << "Using dynamic load balancing" << endl;
    		else cout << "Using static load balancing" << endl;
    	}
    }

    // Create the temporary directory
    if(me == MASTER){
    	int result = mkdir(tmpdir.c_str(), 0777);
    	if(result!=0) {
    		// Try deleting the directory and then making it
    		cout << "Deleting existing temporary directory '"<< tmpdir << "' to overwrite..." << endl;
    		system(("rm -rf " + tmpdir).c_str());
    		int result = mkdir(tmpdir.c_str(), 0777);
    	}
    	assert(result==0);
    }


    timer.exit("input");
    

    //=========================================================//
    // Come up with a name for the files

    string basname = obs->name();
    string abasname = auxbs->name();
    string output_dir = keyval->stringvalue("output_dir");
    string molname = keyval->stringvalue("molname");
    string prefix = output_dir + "/" + molname + "_" + basname + "_" + abasname + "_";
    prefix += "allmax_";

    timer.exit("setup");

    //=========================================================//
    // Construct the CLHF object

    Ref<CLHF> hf = require_dynamic_cast<CLHF*>(
            keyval->describedclassvalue("hf").pointer(), "main\n");
    assert(hf.nonnull());
    Ref<Integral> integral = hf->integral();
    assert(integral.nonnull());
    Ref<PetiteList> pl = integral->petite_list();

    RefSymmSCMatrix P, Q, O;
    Ref<LocalSCMatrixKit> kit = new LocalSCMatrixKit();
    RefSCMatrix Ccabs;
    bool need_P = true;
    bool need_Q = true;
    bool need_O = do_ri;

    // Only get P and Q if we need them.
    if(need_P || need_Q || need_O){
    	// Get P
    	timer.enter("HF");
		RefSymmSCMatrix Pso = hf->density();
		timer.exit("HF");
		Pso.scale(0.5);
		RefSymmSCMatrix Ptmp = pl->to_AO_basis(Pso);
		SCMatrixKit::set_default_matrixkit(kit);
		assert(Ptmp.nblock() == 1);
		P = kit->symmmatrix(obs->basisdim());
		P.assign(require_dynamic_cast<ReplSymmSCMatrix*>(Ptmp.block(0), "main\n")->get_data());

		if(need_Q){
			// Get Q
			timer.enter("Q");
			RefSymmSCMatrix S = hf->overlap();
			RefSymmSCMatrix Sinv = S.i();
			RefSymmSCMatrix Qso = Sinv - Pso;
			RefSymmSCMatrix Qtmp = pl->to_AO_basis(Qso);
			assert(Qtmp.nblock() == 1);
			Q = kit->symmmatrix(obs->basisdim());
			Q.assign(require_dynamic_cast<ReplSymmSCMatrix*>(Qtmp.block(0), "main\n")->get_data());
			timer.exit("Q");
		}

    }

    //=========================================================//
    // make R12WavefunctionWorld
    timer.enter("build cf");
	Ref<WavefunctionWorld> world = require_dynamic_cast<WavefunctionWorld*>(
			keyval->describedclassvalue("world").pointer(), "main\n");
	Ref<RefWavefunction> refwfn = new SD_RefWavefunction(world, hf);
	Ref<R12WavefunctionWorld> r12world = new R12WavefunctionWorld(keyval, refwfn);
	Ref<R12Technology> r12tech = r12world->r12tech();
	Ref<OrbitalSpace> ribs_space = r12world->ribs_space();
	Ref<GaussianBasisSet> ribs = ribs_space->basis();
	Ref<R12Technology::CorrelationFactor> cf = r12tech->corrfactor();
	Ref<Integral> ri_integral = ribs_space->integral();
	Ref<OrbitalSpace> cabs_space = r12world->cabs_space(Alpha);
    timer.exit("build cf");

    // Get O if needed
    if(need_O){
    	// Get O
    	RefSCMatrix Ccabs = cabs_space->coefs();
    	RefSymmSCMatrix O(Ccabs->rowdim(), Ccabs->kit());
    	O.accumulate_symmetric_product(Ccabs);
    }

    //#############################################################################//

    if(do_f12){
    	timer.enter("F12");
    	int num_types = do_eri + do_f12;
    	TwoBodyOper::type otype = cf->tbint_type_f12();
    	string* descs = new string[num_types];
    	// Compute and tranform the integrals
    	Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
    	compute_full_trans_ints(
    			descr, otype,
    			prefix + "PQPQF", tmpdir,
    			obs, obs, obs, obs,
    			P, Q, P, Q
    	);
//    	if(do_ri){
//			// Compute and tranform the integrals
//    		ri_integral->set_basis(obs, ribs, obs, ribs);
//    		DensityMap ri_pairs;
//			ri_pairs["OO"] = SymmSCMatrixPair(O, O);
//			Ref<TwoBodyIntDescr> ri_descr = cf->tbintdescr(ri_integral, 0);
//			compute_half_trans_ints(
//					ri_descr,
//					otypes, descs, num_types,
//					prefix, tmpdir,
//					obs, ribs,
//					ri_pairs,
//					Ccabs->kit()
//			);
//
//    	}
    	delete[] descs;
    	timer.exit("F12");
    }

    //=========================================================//
    // Clean up

    timer.print();

    // Delete the temporary directory
    if(me == MASTER){
    	int error = remove(tmpdir.c_str());
    	if(error != 0){
    		cout << "WARNING: Could not delete temporary directory " << tmpdir << "." << endl;
    	}
    }

    return 0;

}
