/*
 * compute_thread.h
 *
 *  Created on: Jan 7, 2013
 *      Author: dhollman
 */

#include "sparse_ints.h"
#include <pthread.h>
#include <queue>
#include <utility>
#include <map>
#include <math/scmat/local.h>

#ifndef COMPUTE_THREAD_H_
#define COMPUTE_THREAD_H_

// macros
#define DBG_MSG(mymsg) \
	if(opts.debug){ \
		print_lock->lock(); \
		dbg_out_ << mymsg << std::endl; \
		dbg_out_.flush(); \
		/*std::cout << msg->me() << "." << threadnum_ << ":" << mymsg << std::endl;*/ \
		/*std::cout.flush();*/ \
		print_lock->unlock(); \
	}

namespace sparse_ints{

class FullTransComputeThread;
class FullTransCommThread;
class SendThread;
class ReceiveThread;

typedef std::pair<int,int> IntPair;

// MessageTypes
enum {
	IndexData=991,
	PairData=992,
	NeedsSend=993,
	HaveAllData=994,
	ComputeThreadDone=995,
	DoneSending=996,
	PairAssignment=997,
	QueueHasSpace=998
};

class SparseIntsThread : public sc::Thread {

protected:

	sc::Ref<MessageGrp> msg;
	sc::Ref<ThreadGrp> thr;

    int threadnum_;
    const sc::Ref<sc::GaussianBasisSet> basis1_;
    const sc::Ref<sc::GaussianBasisSet> basis2_;
    const sc::Ref<sc::GaussianBasisSet> basis3_;
    const sc::Ref<sc::GaussianBasisSet> basis4_;
    sc::Ref<sc::LocalSCMatrixKit> kit_;
    std::ofstream dbg_out_;

	SparseIntsThread(
		sc::Ref<MessageGrp> msg,
		sc::Ref<ThreadGrp> thr,
		int num,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::Ref<sc::LocalSCMatrixKit> kit
	) :
		threadnum_(num),
		basis1_(bs1),
		basis2_(bs2),
		basis3_(bs3),
		basis4_(bs4),
		kit_(kit),
		msg(msg),
		thr(thr)
    {
		if(opts.debug) {
			std::stringstream sstr;
			sstr << "_debug." << msg->me() << "." << threadnum_ << ".dat";
			dbg_out_.open(sstr.str().c_str(), std::ios::out);
			DBG_MSG("Debugging enabled.");
		}

    }

	~SparseIntsThread(){
		if(opts.debug)
			dbg_out_.close();
		msg = 0;
		thr = 0;
	}

};

class ComputeThread : public SparseIntsThread {

protected:

    sc::Ref<sc::ThreadLock> lock_;
    sc::Ref<sc::TwoBodyInt> inteval_;

public:

	ComputeThread(
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
	);

	virtual void run() = 0;

};

class UntransComputeThread : public ComputeThread {

    sc::TwoBodyOper::type* otypes_;
    std::vector<std::string> prefixes_;
    int num_types_;
    int* quartets_processed_;

public:


    UntransComputeThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::TwoBodyOper::type* otypes,
		std::vector<std::string> prefixes,
		int num_types,
		sc::Ref<sc::LocalSCMatrixKit>& kit,
		int* quartets_processed
    );

    ~UntransComputeThread() { inteval_ = 0; }

    void run();

};

class HalfTransComputeThread : public ComputeThread {

    sc::TwoBodyOper::type* otypes_;
    std::map<std::string, std::vector<std::string> > prefixes_;
    int num_types_;
    DensityMap dens_pairs_;
    DensityMapIterator mapit_;
    sc::RefSCDimension dim1_;
    sc::RefSCDimension dim2_;
    int* quartets_processed_;

public:


    HalfTransComputeThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::TwoBodyOper::type* otypes,
		std::map<std::string, std::vector<std::string> > prefixes,
		int num_types,
		DensityMap& dens_pairs,
		sc::Ref<sc::LocalSCMatrixKit>& kit,
		int* quartets_processed
    );

    ~HalfTransComputeThread() { inteval_ = 0; }

    void run();

};


class FullTransComputeThread : public ComputeThread {

	sc::TwoBodyOper::type otype_;
	sc::RefSymmSCMatrix& P1_;
	sc::RefSymmSCMatrix& P2_;
	sc::RefSymmSCMatrix& P3_;
	sc::RefSymmSCMatrix& P4_;
	std::string prefix_;
	SendThread* send_thread_;
	ReceiveThread* recv_thread_;
	RefSCDimension bsdim1_;
	RefSCDimension bsdim2_;
	RefSCDimension bsdim3_;
	RefSCDimension bsdim4_;

	int *quartets_computed_;

public:

    FullTransComputeThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::TwoBodyOper::type otype,
		std::string prefix,
		sc::RefSymmSCMatrix& P1,
		sc::RefSymmSCMatrix& P2,
		sc::RefSymmSCMatrix& P3,
		sc::RefSymmSCMatrix& P4,
		sc::Ref<LocalSCMatrixKit>& kit,
		SendThread* send_thread,
		ReceiveThread* recv_thread,
		int* quartets_computed
    );

    int get_pair_assignment(int sh3, int sh4);

    void run();

};


class FullTransCommThread : public SparseIntsThread {

protected:

	int* pair_assignments_;
	// Doesn't work Priority queue of (num_bf_on_node, node_number) pairs
	//std::priority_queue<IntPair, std::vector<IntPair>, std::greater<IntPair> > bf_per_node_;
	int* bf_per_node_;

	FullTransCommThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		int num,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::Ref<sc::LocalSCMatrixKit>& kit
	);

	~FullTransCommThread();

    void master_run();

    int master_pair_assignment(int sh3, int sh4);

public:

	virtual void run() =0;

};


class SendThread : public FullTransCommThread {

	typedef struct {
		int ish1, jsh3;
		int sh4;
		int ndata;
		double* data;
	} DataSendTask;

	typedef struct {
		int ish1, ibf1;
		int jsh3, jbf3;
		int sh4;
		int ndata;
		double* data;
	} BFDataSendTask;

	std::queue<DataSendTask> task_queue_;
	std::queue<BFDataSendTask> bftask_queue_;
	sc::Ref<sc::ThreadLock> queue_lock_;

	sc::Ref<sc::ThreadLock> comm_lock_;

	long queue_size_;
	//size_t queue_size_;
	static size_t max_queue_size;

	const int needs_send_message_;

	std::vector<MessageGrp::MessageHandle> handles_;
	std::vector<int> messages_;

public:
	SendThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::Ref<sc::LocalSCMatrixKit>& kit
	);

	~SendThread();

    void run();

	void distribute_shell_pair(
		std::vector<RefSCMatrix> pair_mats,
		int ish1, int jsh3, int sh4,
		int threadnum
	);

	void distribute_bf_pair(
		sc::RefSCMatrix ints2q,
		int ibf1, int jbf3, int sh4,
		int threadnum
	);


    int get_pair_assignment(int sh3, int sh4);
};

class ReceiveThread : public FullTransCommThread {

	std::vector<IntPair> assigned_pairs_;
	std::vector<int> nbf_pairs_;
	int sorted_pairs_position_;
	sc::Ref<sc::ThreadLock> pairs_lock_;


	class bf_compare_{
		const sc::Ref<sc::GaussianBasisSet> basis3_;
		const sc::Ref<sc::GaussianBasisSet> basis4_;
	public:
		bf_compare_(
			const sc::Ref<sc::GaussianBasisSet>& basis3,
			const sc::Ref<sc::GaussianBasisSet>& basis4
		) :
			basis3_(basis3),
			basis4_(basis4)
		{
		}

		bool operator()(IntPair sh34a, IntPair sh34b) {
			int nbf3a = basis3_->shell(sh34a.first).nfunction();
			int nbf4a = basis4_->shell(sh34a.second).nfunction();
			int nbf3b = basis3_->shell(sh34b.first).nfunction();
			int nbf4b = basis4_->shell(sh34b.second).nfunction();
			return nbf3a*nbf4a > nbf3b*nbf4b;
		}
	};

public:

	bool my_ints2q_complete_;
	std::map<IntPair, std::vector<RefSCMatrix> > my_ints2q;

	ReceiveThread(
		Ref<MessageGrp> msg,
		Ref<ThreadGrp> thr,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		sc::Ref<sc::LocalSCMatrixKit>& kit
	);

    bool get_my_next_pair(int&, int&);

	void run();

};


} // end namespace sparse_ints


#endif /* COMPUTE_THREAD_H_ */
