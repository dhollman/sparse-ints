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

#ifndef COMPUTE_THREAD_H_
#define COMPUTE_THREAD_H_


namespace sparse_ints{

class FullTransComputeThread;
class FullTransCommThread;

typedef std::pair<int,int> IntPair;

// MessageTypes
enum {
	IndexData,
	PairData,
	NeedsSend,
	HaveAllData,
	ComputeThreadDone,
	NeedPairAssignment,
	PairAssignment
};

class ComputeThread :public sc::Thread {

protected:

    int threadnum_;
    sc::Ref<sc::ThreadLock> lock_;
    const sc::Ref<sc::GaussianBasisSet> basis1_;
    const sc::Ref<sc::GaussianBasisSet> basis2_;
    const sc::Ref<sc::GaussianBasisSet> basis3_;
    const sc::Ref<sc::GaussianBasisSet> basis4_;
    sc::Ref<sc::TwoBodyInt> inteval_;
    sc::Ref<sc::SCMatrixKit> kit_;

public:

	ComputeThread(
		int num,
		const sc::Ref<sc::TwoBodyIntDescr>& intdescr,
		const sc::Ref<sc::ThreadLock>& lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4
	);

	virtual void run() = 0;

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
		sc::Ref<sc::SCMatrixKit>& kit,
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
	FullTransCommThread* send_thread_;
	FullTransCommThread* recv_thread_;
	RefSCDimension bsdim1_;
	RefSCDimension bsdim2_;
	RefSCDimension bsdim3_;
	RefSCDimension bsdim4_;


public:

    FullTransComputeThread(
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
		sc::Ref<SCMatrixKit>& kit,
		FullTransCommThread* send_thread,
		FullTransCommThread* recv_thread
    );

    void run();

};


class FullTransCommThread : public sc::Thread {

	typedef struct {
		int sh1, sh2;
		int nbf1, nbf2, nbftot3, nbftot4;
		double* data;
	} DataSendTask;

	sc::Ref<sc::ThreadLock> comm_lock_;
	sc::Ref<sc::ThreadLock> queue_lock_;

	std::queue<DataSendTask> task_queue_;


	size_t queue_size_;

	const sc::Ref<sc::GaussianBasisSet> basis1_;
	const sc::Ref<sc::GaussianBasisSet> basis2_;
	const sc::Ref<sc::GaussianBasisSet> basis3_;
	const sc::Ref<sc::GaussianBasisSet> basis4_;

	static size_t max_queue_size;

	int thread_type_;
	Ref<SCMatrixKit>& kit_;


	// For use by MASTER
	std::map<IntPair, int> pair_assignments_;
	// Priority queue of (num_bf_on_node, node_number) pairs
	std::priority_queue<IntPair, std::vector<IntPair>, std::greater<IntPair> > bf_per_node_;

	bool my_ints2q_complete_;

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

	// This should be a queue, but sorting was a pain so I just used a vector
	std::vector<IntPair> pairs_sorted_;
	int sorted_pairs_position_;

public:
	enum {
		ReceiveThread,
		SendThread
	};

	std::map<IntPair, std::vector<RefSCMatrix> > my_ints2q;

	FullTransCommThread(
		const sc::Ref<sc::ThreadLock>& comm_lock,
		const sc::Ref<sc::ThreadLock>& queue_lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		int thread_type,
		Ref<SCMatrixKit>& kit_
	);

	void distribute_shell_pair(
		std::vector<RefSCMatrix> pair_mats,
		int sh1, int sh2, int nbf1, int nbf2, int nbf3tot, int nbf4tot,
		FullTransComputeThread* compute_thread
	);

    void run();

    int get_pair_assignment(int sh3, int sh4);

    bool get_my_next_pair(int&, int&);

};


class SendThread : public FullTransCommThread {

public:
	SendThread(
		const sc::Ref<sc::ThreadLock>& comm_lock,
		const sc::Ref<sc::ThreadLock>& queue_lock,
		const sc::Ref<sc::GaussianBasisSet>& bs1,
		const sc::Ref<sc::GaussianBasisSet>& bs2,
		const sc::Ref<sc::GaussianBasisSet>& bs3,
		const sc::Ref<sc::GaussianBasisSet>& bs4,
		Ref<SCMatrixKit>& kit_
	);

};



} // end namespace sparse_ints


#endif /* COMPUTE_THREAD_H_ */
