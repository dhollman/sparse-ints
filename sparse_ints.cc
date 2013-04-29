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
#include "compute_full.h"
#include "compute_hti.h"
#include "compute_untrans.h"
#include "binfiles.h"


using namespace sparse_ints;
using namespace sc;
using namespace std;

// Global variables
MultiTimer sparse_ints::timer;
SparseIntOptions sparse_ints::opts;
Ref<ThreadLock> sparse_ints::print_lock;


int
main(int argc, char** argv) {

    /*=========================================================*/
    /* Create msg, thr, and timer                              */ #if fold_begin

    //---------------------------------------------------------//
    // Create the message group

    Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc,argv);
    if (msg.null()) {
		#if HAVE_MPI_H
        msg = new MPIMessageGrp();
		#else
        msg = new ProcMessageGrp();
		#endif
    }
    MessageGrp::set_default_messagegrp(msg);
    int n = msg->n();
    int me = msg->me();

    //---------------------------------------------------------//
    // Create the thread group object

    Ref<ThreadGrp> thr = ThreadGrp::initial_threadgrp(argc, argv);
    if (thr.null()) {
    	thr = ThreadGrp::get_default_threadgrp();
    }
    ThreadGrp::set_default_threadgrp(thr);
    int nthr = thr->nthread();
    print_lock = thr->new_lock();


    if (!opts.quiet && me == MASTER) {
        cout << "Running on " << n << " node"
             << (n == 1 ? "" : "s") << " with " << nthr << " thread"
             << (nthr == 1 ? "" : "s") << " per node." << endl;
    }

    //---------------------------------------------------------//
    // Create the timer
    
    MultiTimer timer("Sparse Ints");
    sparse_ints::timer = timer;

    /***********************************************************/ #endif
    /*=========================================================*/
    /* Get the stuff from the input                            */ #if fold_begin

    timer.enter("setup");
    timer.enter("input");
    // Read the input file
    char *infile = new char[256];
    if (argc == 2) {
        infile = argv[1];
    } else {
        sprintf(infile, "input.dat");
    }
    Ref<KeyVal> pkv(new ParsedKeyVal(infile));

    // Create keyvalue objects for each of the sections
    bool do_untrans, do_halftrans, do_fulltrans;
    Ref<KeyVal> keyval(new PrefixKeyVal(pkv, ":all"));

    // Get the basis set
    Ref<GaussianBasisSet> obs = require_dynamic_cast<GaussianBasisSet*>(
            keyval->describedclassvalue("basis").pointer(), "main\n");
    // Get rid of general constractions for simplicity
    if (obs->max_ncontraction() > 1) {
        Ref<GaussianBasisSet> split_basis = new SplitBasisSet(obs, obs->name());
        obs = split_basis;
    }

    // Get the auxiliary basis set
    Ref<GaussianBasisSet> auxbs = require_dynamic_cast<GaussianBasisSet*>(
            keyval->describedclassvalue("aux_basis").pointer(), "main\n");
    // Get rid of general constractions for simplicity
    if (auxbs->max_ncontraction() > 1) {
        Ref<GaussianBasisSet> split_basis = new SplitBasisSet(auxbs, auxbs->name());
        auxbs = split_basis;
    }

    // Get the molecule object
    Ref<Molecule> mol = obs->molecule();

    // Get some boolean options
    SparseIntOptions opts;
    opts.debug = keyval->booleanvalue("debug", KeyValValueboolean(false));
    if(opts.debug && me == MASTER) cout << "Debugging enabled." << endl;
    opts.verbose = keyval->booleanvalue("verbose", KeyValValueboolean(false));
    opts.quiet = keyval->booleanvalue("quiet", KeyValValueboolean(false));
    opts.dynamic = keyval->booleanvalue("dynamic", KeyValValueboolean(true));
    opts.max_only = keyval->booleanvalue("max_only", KeyValValueboolean(true));

    // For debugging purposes, force the computation
    //   of various types of "fake" integrals to try
    //   and isolate problems with the transformation
    opts.use_fake_ints = keyval->intvalue("use_fake_ints", KeyValValueint(0));


    //Set the output type of the integrals
    if(opts.max_only)
    	opts.out_type = MaxAbs;
    else
    	opts.out_type = AllInts;
    sparse_ints::opts = opts;

    // Get which integrals to do
    bool do_eri = keyval->booleanvalue("do_eri", KeyValValueboolean(true));
    bool do_f12 = keyval->booleanvalue("do_f12", KeyValValueboolean(true));
    bool do_f12g12 = keyval->booleanvalue("do_f12g12", KeyValValueboolean(false));
    bool do_f12sq = keyval->booleanvalue("do_f12sq", KeyValValueboolean(false));
    bool do_dblcomm = keyval->booleanvalue("do_dblcomm", KeyValValueboolean(false));

    // Get which (extra) density matrices to do, even if we don't need them
    bool do_P = keyval->booleanvalue("do_P", KeyValValueboolean(false));
    bool do_Q = keyval->booleanvalue("do_Q", KeyValValueboolean(false));
    bool do_O = keyval->booleanvalue("do_O", KeyValValueboolean(false));

    /*-----------------------------------------------------*/
	/*  Options specific to untrans int code          {{{2 */ #if fold_begin

    Ref<KeyVal> untranskv(new PrefixKeyVal(pkv, ":untrans"));
    do_untrans = pkv->exists(":untrans");
    bool use_ribs[4];
    if(do_untrans){
    	for_each(nriq, 4){
			use_ribs[nriq] = untranskv->booleanvalue("do_ri", nriq, KeyValValueboolean(false));
    	}
    }

    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/
	/*  Options specific to half trans code           {{{2 */ #if fold_begin

    Ref<KeyVal> halftranskv(new PrefixKeyVal(pkv, ":halftrans"));
    do_halftrans = pkv->exists(":halftrans");

    // Get the list of density matrix transformations to do.
    int num_densmats_half = 0;
    vector<string> densmats_half;
    if(do_halftrans){
		num_densmats_half = halftranskv->count("densmats");
		densmats_half.resize(num_densmats_half);
		for_each(imat, num_densmats_half){
			densmats_half[imat] = halftranskv->stringvalue("densmats", imat, KeyValValuestring("PQ"));
		}
    }

    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/
	/* Options specific to full trans code            {{{2 */ #if fold_begin

    Ref<KeyVal> fulltranskv(new PrefixKeyVal(pkv, ":fulltrans"));
    do_fulltrans = pkv->exists(":fulltrans");

    // Get the list of density matrix transformations to do.
    int num_densmats_full = 0;
    vector<string> densmats_full;
    if(do_fulltrans){
		num_densmats_full = fulltranskv->count("densmats");
		densmats_full.resize(num_densmats_full);
		for_each(imat, num_densmats_full){
			densmats_full[imat] = fulltranskv->stringvalue("densmats", imat, KeyValValuestring("PQPQ"));
		}
    }
    /*******************************************************/ #endif //2}}}
    /*-----------------------------------------------------*/


    /***********************************************************/ #endif
    /*=========================================================*/
    /* Make the scratch directory							   */ #if fold_begin

    char* scratch_dir = getenv("SCRATCH");
    bool alloc_scratch_dir = false, alloc_job_id = false;
    if(scratch_dir == 0){
    	scratch_dir = new char[32];
        alloc_scratch_dir = true;
        sprintf(scratch_dir, "/tmp");
    }
    char* job_id = getenv("PBS_JOBID");
    int pid;
    if(job_id == 0){
    	if(me == MASTER){
			// get unique process id.  (Only if we can't get a JOBID)
			pid_t mypid = getpid();
			pid = (int)mypid;
    	}
    	// broadcast the process id
    	msg->bcast(pid, 0);
    	job_id = new char[32];
        alloc_job_id = true;
    	sprintf(job_id, "%d", pid);
    }
	// Get the temporary directory as a string
    char* ctmpdir = new char[256];
    sprintf(ctmpdir, "%s/%s", scratch_dir, job_id);
    // Pedantically avoid leaking memory
    if(alloc_scratch_dir) delete[] scratch_dir;
    if(alloc_job_id) delete[] job_id;
    string tmpdir(const_cast<const char*>(ctmpdir));
    if(!opts.quiet){
    	if((!opts.debug && me == MASTER) || opts.debug){
    		cout << "Using temporary directory " << tmpdir << endl;
    		if(opts.dynamic) cout << "Using dynamic load balancing" << endl;
    		else cout << "Using static load balancing" << endl;
    	}
    }
    delete[] ctmpdir;

    // Create the temporary directory
    if(me == MASTER){
    	int result = mkdir(tmpdir.c_str(), 0777);
    	int num_retries = 0;
    	while(num_retries < 15 && result != 0){
    		num_retries++;
    		// Try deleting the directory and then making it
    		cout << "Deleting existing temporary directory '"<< tmpdir << "' to overwrite..." << endl;
    		system(("rm -rf " + tmpdir).c_str());
    		sleep(1);
    		int result = mkdir(tmpdir.c_str(), 0777);
    	}
    	assert(result==0);
    }

    timer.exit("input");

    /***********************************************************/ #endif
    /*=========================================================*/
    /* Come up with an name for the output files			   */ #if fold_begin

    string basname = obs->name();
    string abasname = auxbs->name();
    string output_dir = keyval->stringvalue("output_dir");
    string molname = keyval->stringvalue("molname");
    string prefix = output_dir + "/" + molname + "_" + basname + "_" + abasname + "_";
    string densprefix = prefix + "";
    if(opts.out_type == MaxAbs)
		prefix += "allmax_";
    else if(opts.out_type == AllInts)
    	prefix += "allints_";
    else if(opts.out_type == Median)
    	prefix += "median_";
    else if(opts.out_type == Average)
    	prefix += "average_";
    else if(opts.out_type == StdDev)
    	prefix += "stddev_";
    else
    	assert(not_implemented);

    if(opts.use_fake_ints){
    	stringstream sstr;
    	sstr << prefix << "fake" << opts.use_fake_ints << "_";
    	prefix = sstr.str();
    }

    string untrans_prefix = prefix + "";
    if(do_untrans){
    	int bsval = 0;
    	string bs_indicator = "bs";
    	bool need_write = false;
    	for_each(iri, 4){
			if(use_ribs[iri]){
				need_write = true;
				bs_indicator += "r";
			}
			else{
				bs_indicator += "o";
			}
    	}
    	if(need_write){
    		untrans_prefix += bs_indicator + "_";
    	}
    }

    timer.exit("setup");

    /***********************************************************/ #endif
    /*=========================================================*/
    /* Construct the CLHF object 							   */ #if fold_begin

    //Ref<Integral> integral = Integral::initial_integral(argc, argv);
    //if(integral.null()){
    	//integral = Integral::get_default_integral();
    //}
    //assert(integral.nonnull());

    Ref<CLHF> hf = require_dynamic_cast<CLHF*>(
            keyval->describedclassvalue("hf").pointer(), "main\n");
    assert(hf.nonnull());
    Ref<Integral> integral = hf->integral();
    assert(integral.nonnull());
    Ref<PetiteList> pl = integral->petite_list();

    /***********************************************************/ #endif
    /*=========================================================*/
    /* Make R12WavefunctionWorld                               */ #if fold_begin

    timer.enter("build cf");
	Ref<WavefunctionWorld> world = require_dynamic_cast<WavefunctionWorld*>(
			keyval->describedclassvalue("world").pointer(), "main\n");
	Ref<RefWavefunction> refwfn = new SD_RefWavefunction(world, hf);
	Ref<R12WavefunctionWorld> r12world = new R12WavefunctionWorld(keyval, refwfn);
	r12world->initialize();
	Ref<R12Technology> r12tech = r12world->r12tech();
	Ref<OrbitalSpace> ribs_space = r12world->ribs_space();
	Ref<GaussianBasisSet> ribs = ribs_space->basis();
	Ref<R12Technology::CorrelationFactor> cf = r12tech->corrfactor();
	Ref<Integral> ri_integral = ribs_space->integral();
	Ref<OrbitalSpace> cabs_space = r12world->cabs_space(Alpha);
    timer.exit("build cf");

    /***********************************************************/ #endif
    /*=========================================================*/
    /* Determine which density matrices we need and get them   */ #if fold_begin

    // Determine which density matrices we need
    RefSymmSCMatrix P, Q, O, I, R, zero;
    Ref<LocalSCMatrixKit> kit = new LocalSCMatrixKit();
    RefSCMatrix Ccabs;
    bool need_P, need_Q, need_O, need_I, need_R, need_0;
    need_P = need_Q = need_O = need_I = need_R = need_0 = false;
    if(do_halftrans){
		for_each(iset, num_densmats_half){
			for_each(imat, 2){
				switch(densmats_half[iset][imat]) {
				case 'P':
					need_P = true;
					break;
				case 'Q':
					need_Q = true;
					break;
				case 'O':
					need_O = true;
					break;
				case 'I':
					need_I = true;
					break;
				case 'R':
					need_R = true;
					break;
				case '0':
					need_0 = true;
					break;
				}
			}
		}
    } // end if do_halftrans
    if(do_fulltrans){
		for_each(iset, num_densmats_full){
			for_each(imat, 4){
				switch(densmats_full[iset][imat]) {
				case 'P':
					need_P = true;
					break;
				case 'Q':
					need_Q = true;
					break;
				case 'O':
					need_O = true;
					break;
				case 'I':
					need_I = true;
					break;
				case 'R':
					need_R = true;
					break;
				case '0':
					need_0 = true;
					break;
				}
			}
		}
    }

    need_P = need_P || do_P;
    need_Q = need_Q || do_Q;
    need_O = need_O || do_O;

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
    else {
		SCMatrixKit::set_default_matrixkit(kit);
    }


    // Get O if needed
    if(need_O){
    	// Get O
    	RefSCMatrix Ccabs = cabs_space->coefs();
		RefSymmSCMatrix Otmp = Ccabs->kit()->symmmatrix(Ccabs->rowdim());
		Otmp.assign(0.0);
		Otmp.accumulate_symmetric_product(Ccabs);
    	assert(Otmp.nblock() == 1);
    	O = kit->symmmatrix(ribs->basisdim());
    	O.assign(require_dynamic_cast<ReplSymmSCMatrix*>(Otmp.block(0), "main\n")->get_data());
    }

    // Make an identity transformation matrices and the zero transformation matrix.  Mostly for testing...
    if(need_I){
		I = kit->symmmatrix(obs->basisdim());
		I.assign(0.0);
		I->shift_diagonal(1.0);
    }
    if(need_R){
		R = kit->symmmatrix(ribs->basisdim());
		R.assign(0.0);
		R->shift_diagonal(1.0);
    }
    if(need_0){
    	zero = kit->symmmatrix(obs->basisdim());
    	zero.assign(0.0);
    }

    // Print the density matrices to binary files if we have them
    if(need_P && me == MASTER){
    	write_density_binfile(P, densprefix + "P.bin", obs);
    }
    if(need_Q && me == MASTER){
    	write_density_binfile(Q, densprefix + "Q.bin", obs);
    }
    if(need_O && me == MASTER){
    	write_density_binfile(O, densprefix + "O.bin", ribs);
    }

    /***********************************************************/ #endif
    /*=========================================================*/
    /*#########################################################*/
    /*=========================================================*/
    /* Do the untransformed integrals if we need to            */ #if fold_begin

    if(do_untrans){
    	if(me == MASTER){
			cout << "========================================================================" << endl;
			cout << "=                  Computing Untransformed Integrals                   =" << endl;
			cout << "========================================================================" << endl;
			cout << endl;
    	}
		Ref<GaussianBasisSet> basis_sets[4];
		for_each(iri, 4){
			if(use_ribs[iri])
				basis_sets[iri] = ribs;
			else
				basis_sets[iri] = obs;
		}
		integral->set_basis(basis_sets[0], basis_sets[1], basis_sets[2], basis_sets[3]);

		if(do_eri || do_f12 || do_f12g12){
			timer.enter("untrans ERI/F12/F12G12");
			int num_types = do_eri + do_f12 + do_f12g12;
			TwoBodyOper::type otypes[num_types];
			string* descs = new string[num_types];
			int ity = 0;
			if(do_eri){
				otypes[ity] = cf->tbint_type_eri();
				descs[ity] = "g";
				++ity;
			}
			if(do_f12){
				otypes[ity] = cf->tbint_type_f12();
				descs[ity] = "F";
				++ity;
			}
			if(do_f12g12){
				otypes[ity] = cf->tbint_type_f12eri();
				descs[ity] = "Fg";
				++ity;
			}
			Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);

			compute_untrans_threaded(msg, thr,
					descr,
					otypes, descs, num_types,
					untrans_prefix, tmpdir,
					basis_sets[0], basis_sets[1], basis_sets[2], basis_sets[3],
					kit
			);
			timer.exit("untrans ERI/F12/F12G12");
		}

		if(do_f12sq || do_dblcomm){
			timer.enter("untrans F12sq/DC");
			int num_types = do_f12sq + do_dblcomm;
			TwoBodyOper::type otypes[num_types];
			string* descs = new string[num_types];
			int ity = 0;
			if(do_f12sq){
				otypes[ity] = cf->tbint_type_f12f12();
				descs[ity] = "Fsq";
				++ity;
			}
			if(do_dblcomm){
				otypes[ity] = cf->tbint_type_f12t1f12();
				descs[ity] = "DC";
				++ity;
			}
			Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0, 0);

			compute_untrans_threaded(msg, thr,
					descr,
					otypes, descs, num_types,
					untrans_prefix, tmpdir,
					basis_sets[0], basis_sets[1], basis_sets[2], basis_sets[3],
					kit
			);
			timer.exit("untrans F12sq/DC");
		}
    }

    /***********************************************************/ #endif
    /*=========================================================*/
    /* Do the half transformed integrals if we need to         */ #if fold_begin
    if(do_halftrans){
    	if(me == MASTER){
			cout << "========================================================================" << endl;
			cout << "=                Computing Half Transformed Integrals                  =" << endl;
			cout << "========================================================================" << endl;
			cout << endl;
    	}
		if(do_eri || do_f12 || do_f12g12){
			timer.enter("half trans ERI/F12");
			int num_types = do_eri + do_f12 + do_f12g12;
			TwoBodyOper::type otypes[num_types];
			string* descs = new string[num_types];
			int ity = 0;
			if(do_eri){
				otypes[ity] = cf->tbint_type_eri();
				descs[ity] = "g";
				++ity;
			}
			if(do_f12){
				otypes[ity] = cf->tbint_type_f12();
				descs[ity] = "F";
				++ity;
			}
			if(do_f12g12){
				otypes[ity] = cf->tbint_type_f12eri();
				descs[ity] = "Fg";
				++ity;
			}

			DensityMap dens_pairs_oo, dens_pairs_rr, dens_pairs_or;
			bool do_oo, do_rr, do_or;
			do_oo = do_rr = do_or = false;

			// Break the computed pairs into sets: obs-obs, obs-ribs, ribs-ribs
			for_each(iset, num_densmats_half){
				Ref<GaussianBasisSet> basis_sets[2];
				RefSymmSCMatrix dmats[2];
				for_each(imat, 2){
					switch(densmats_half[iset][imat]) {
					case 'P':
						dmats[imat] = P;
						basis_sets[imat] = obs;
						break;
					case 'Q':
						dmats[imat] = Q;
						basis_sets[imat] = obs;
						break;
					case 'O':
						dmats[imat] = O;
						basis_sets[imat] = ribs;
						break;
					case 'I':
						dmats[imat] = I;
						basis_sets[imat] = obs;
						break;
					case 'R':
						dmats[imat] = R;
						basis_sets[imat] = ribs;
						break;
					case '0':
						dmats[imat] = zero;
						basis_sets[imat] = obs;
						break;
					default:
						cout << "Invalid density matrices: " << densmats_half[iset] << endl;
						sleep(1);
						assert(false);
						break;
					}
				}
				if(basis_sets[0] == obs && basis_sets[1] == obs){
					dens_pairs_oo[densmats_half[iset]] = SymmSCMatrixPair(dmats[0], dmats[1]);
					do_oo = true;
				}
				else if(basis_sets[0] == ribs && basis_sets[1] == ribs){
					dens_pairs_rr[densmats_half[iset]] = SymmSCMatrixPair(dmats[0], dmats[1]);
					do_rr = true;
				}
				else if(basis_sets[0] == obs && basis_sets[1] == ribs){
					dens_pairs_or[densmats_half[iset]] = SymmSCMatrixPair(dmats[0], dmats[1]);
					do_or = true;
					assert(not_implemented);
				}
				else{
					assert(not_implemented);
				}

			}
			if(do_oo){
				integral->set_basis(obs, obs, obs, obs);
				Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
				compute_hti_threaded(msg, thr,
						descr,
						otypes, descs, num_types,
						prefix, tmpdir,
						obs, obs,
						dens_pairs_oo,
						kit
				);
			}
			if(do_rr){
				integral->set_basis(ribs, ribs, ribs, ribs);
				Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
				compute_hti_threaded(msg, thr,
						descr,
						otypes, descs, num_types,
						prefix, tmpdir,
						ribs, ribs,
						dens_pairs_rr,
						kit
				);
			}
			if(do_or){
				assert(not_implemented);
				Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
				compute_hti_threaded(msg, thr,
						descr,
						otypes, descs, num_types,
						prefix, tmpdir,
						obs, ribs,
						dens_pairs_or,
						kit
				);
			}
			timer.exit("halftrans ERI/F12");
		}
    }
    /***********************************************************/ #endif
    /*=========================================================*/
    /* Do fully transformed integrals if we need to            */ #if fold_begin

    // Note:  This doesn't reuse the half-transformed integrals,
    //  so it's better to do this in a separate computation if
    //  you need both

    RefSymmSCMatrix dmats[4];
    Ref<GaussianBasisSet> basis_sets[4];
    if(do_fulltrans){
    	if(me == MASTER){
			cout << "========================================================================" << endl;
			cout << "=               Computing Fully Transformed Integrals                  =" << endl;
			cout << "========================================================================" << endl;
			cout << endl;
    	}
		for_each(iset, num_densmats_full){
			if(me == MASTER){
				cout << "========================================" << endl;
				cout << "  Working on transformation " << densmats_full[iset] << endl;
				cout << "========================================" << endl;
			}
			for_each(imat, 4){
				switch(densmats_full[iset][imat]) {
				case 'P':
					dmats[imat] = P;
					basis_sets[imat] = obs;
					break;
				case 'Q':
					dmats[imat] = Q;
					basis_sets[imat] = obs;
					break;
				case 'O':
					dmats[imat] = O;
					basis_sets[imat] = ribs;
					break;
				case 'I':
					dmats[imat] = I;
					basis_sets[imat] = obs;
					break;
				case 'R':
					dmats[imat] = R;
					basis_sets[imat] = ribs;
					break;
				case '0':
					dmats[imat] = zero;
					basis_sets[imat] = obs;
					break;
				default:
					cout << densmats_full[iset] << endl;
					sleep(1);
					assert(false);
					break;
				}
			}

			integral->set_basis(basis_sets[0], basis_sets[1], basis_sets[2], basis_sets[3]);

			if(do_eri){
				timer.enter("ERI");
				TwoBodyOper::type otype = cf->tbint_type_eri();
				// Compute and tranform the integrals
				Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
				compute_full_trans_ints(msg, thr,
						descr, otype,
						prefix, densmats_full[iset] + "g", tmpdir,
						basis_sets[0], basis_sets[1], basis_sets[2], basis_sets[3],
						dmats[0], dmats[1], dmats[2], dmats[3],
						kit
				);
				timer.exit("ERI");
			}
			if(do_f12){
				timer.enter("F12");
				TwoBodyOper::type otype = cf->tbint_type_f12();
				// Compute and tranform the integrals
				Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
				compute_full_trans_ints(msg, thr,
						descr, otype,
						prefix, densmats_full[iset] + "F", tmpdir,
						basis_sets[0], basis_sets[1], basis_sets[2], basis_sets[3],
						dmats[0], dmats[1], dmats[2], dmats[3],
						kit
				);
				timer.exit("F12");
			}
			if(do_f12g12){
				timer.enter("F12G12");
				TwoBodyOper::type otype = cf->tbint_type_f12eri();
				// Compute and tranform the integrals
				Ref<TwoBodyIntDescr> descr = cf->tbintdescr(integral, 0);
				compute_full_trans_ints(msg, thr,
						descr, otype,
						prefix, densmats_full[iset] + "Fg", tmpdir,
						basis_sets[0], basis_sets[1], basis_sets[2], basis_sets[3],
						dmats[0], dmats[1], dmats[2], dmats[3],
						kit
				);
				timer.exit("F12G12");
			}
		}
    } // end if do_fulltrans

    /***********************************************************/ #endif // end fulltrans fold
    /*=========================================================*/
    /*#########################################################*/
    /*=========================================================*/
    /* Clean up                                                */ #if fold_begin

    timer.print();

    // Delete the temporary directory
    if(me == MASTER){
    	int error = remove(tmpdir.c_str());
    	if(error != 0){
    		cout << "WARNING: Could not delete temporary directory " << tmpdir << "." << endl;
    	}
    }

    return 0;

    /***********************************************************/ #endif
    /*=========================================================*/

}
