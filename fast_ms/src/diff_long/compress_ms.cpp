/*
 compress a ms vector
*/

#include <cstdint>
#include <vector>
#include <iostream>
#include <fstream> 
#include <cstdio>
//#include <stdlib>

#include "sdsl/bit_vectors.hpp"
#include "sdsl/util.hpp"

#include "../fd_ms/help.hpp"
#include "../fd_ms/input_spec.hpp"
#include "../fd_ms/opt_parser.hpp"
#include "../fd_ms/counter.hpp"
#include "../fd_ms/p_ms_vector.hpp"
#include "../fd_ms/stree_sct3.hpp"

#include "rlcsa/bits/bitvector.h"
#include "rlcsa/bits/rlevector.h"
#include "rlcsa/bits/deltavector.h"
#include "rlcsa/bits/succinctvector.h"
#include "rlcsa/bits/nibblevector.h"

#include "../../malloc_count/malloc_count.h"

using namespace fdms;
using namespace std;

typedef StreeOhleb<> cst_t;
typedef typename cst_t::size_type size_type;
typedef typename ms_compression::compression_types Compression;
typedef std::map<size_type, size_type> histo_t;


class InputFlags {
public:
    size_type threshold, n_zeros, n_ones;
    bool negative, greedy, verbose;

    InputFlags() { }

    InputFlags(const InputFlags& f) :
        threshold{f.threshold}, n_zeros{f.n_zeros}, n_ones{f.n_ones},
        negative{f.negative}, greedy{f.greedy}, verbose{f.verbose}
    { }

    InputFlags(OptParser input) :
        threshold{static_cast<size_type> (std::stoll(input.getCmdOption("-threshold")))},
        n_zeros{static_cast<size_type> (std::stoll(input.getCmdOption("-nZeros")))},
        n_ones{static_cast<size_type> (std::stoll(input.getCmdOption("-nOnes")))},
        negative{input.getCmdOption("-negative") == "1"},
        greedy{input.getCmdOption("-greedy") == "1"},
        verbose{input.getCmdOption("-verbose") == "1"}
    {}
};


/* read current memory usage */
size_type abs_point() {
    malloc_count_reset_peak();
    return (size_type) malloc_count_peak();
}

/* difference between given and current memory usage */
size_type diff_from(const size_type from){
    size_type to = abs_point();
    if (from > to)
        throw string{"from peak (" + to_string(from) + ") < to peak (" + to_string(to) + ")"};
    return (size_type) (to - from);
}

/**
 * 
 */
template<typename enc_type>
size_type fill_encoder(sdsl::bit_vector ms, enc_type& encoder, histo_t& freq){
    size_type no = 0, i = 0, n_runs = 0;
    while(i < ms.size()){
        if(ms[i] == 1)
            no += 1;
        else{
            assert (i >= no);
            if (no > 0){
                n_runs += 1;
                encoder.addRun(i - no, no);
                freq[no] += 1;
                no = 0;
            }
        }
        i += 1;
    }
    if(no > 0){
        encoder.addRun(i - no, no);
        n_runs += 1;
        freq[no] += 1;
    }
    encoder.flush();
    return n_runs;
}








// ----------------------------------- PERMUTATION ---------------------------------------

int64_t **neighbor;
size_type maxZ, maxO;
size_type out[5];  // Temporary array


/**
 * @param nZeros,nOnes >=1.
 */
size_type coordinates2id(const size_type nZeros, const size_type nOnes, const size_type totalNZeros, const size_type totalNOnes) {
	return (nZeros-1)*totalNOnes+(nOnes-1);
}


void id2coordinates(const size_type nodeID, const size_type totalNZeros, const size_type totalNOnes, size_type out[]) {
	out[0]=1+nodeID/totalNOnes;
	out[1]=1+nodeID%totalNOnes;
}


size_type length(size_type n) {
  size_type b = 0;
  while(n > 0) { b++; n >>= 1; }
  return b;
}


size_type deltaCodeLength(const size_type value) {
	size_type len = length(value);
	size_type llen = length(len);
	return (len + llen + llen - 2);
}


/**
 * @return number of bits for encoding the run pair $(nZeros,nOnes)$.
 */
size_type encodingSize(const size_type nZeros, const size_type nOnes) {
	return deltaCodeLength(nZeros)+deltaCodeLength(nOnes);
}


/**  
 * Let $[i..j]$ be such that $MS[i] \geq \tau$, $MS[j] \geq \tau$, and $MS[k]<\tau$ for 
 * all $k \in [i+1..j-1]$. Permutation should be performed on a range $[p..q]$ such that 
 * $p$ is the first zero-bit to the right of $i$, and $q$ is the first one-bit to the left
 * of $j$. This choice of $p,q$ is done in order not to alter the length of the pairs 
 * $(nZeros_L,nOnes_L),(nZeros_R,nOnes_R)$ that intersect $[i..j]$ on the left and on the 
 * right, which might be similar when MS is large. 
 *
 * @param windowStart,windowEnd values $i,j$ above;
 * @param out output array that describes the permuted range: 0=start; 1=end; 2=nZeros; 
 * 3=nOnes; 4=MS value of the last one-bit before $start$;
 * @return 0 if no range should be permuted, 1 otherwise.
 */
uint8_t getPermutationRange_window(sdsl::bit_vector &ms, const size_type windowStart, const size_type windowEnd, const size_type windowOnes, const size_type threshold, size_type out[]) {
	size_type i, p, q;
	size_type delta, remainingOnes, initialMSValue;
	
	// Not altering the prefix of ones of $[windowStart+1..windowEnd-1]$.
	p=windowStart+1;
	while (ms[p]==1) p++;
	out[0]=p;
	delta=p-windowStart-1;
	initialMSValue=threshold-delta;
	remainingOnes=windowOnes-delta;
	if (remainingOnes==0) return 0;

	// Not altering the suffix of zeros of $[windowStart+1..windowEnd-1]$.
	q=windowEnd-1;
	while (ms[q]==0) q--;
	out[1]=q; out[2]=q-p+1-remainingOnes; out[3]=remainingOnes; out[4]=initialMSValue;
	return 1;
}


/** 
 * Like $getPermutationRange_window()$, but for the case in which $i=-1$.
 */
uint8_t getPermutationRange_firstWindow(sdsl::bit_vector &ms, const size_type windowEnd, const size_type windowOnes, size_type out[]) {
	const size_type initialMSValue = 1;
	size_type i, q;
	size_type remainingOnes;
	
	if (windowOnes==0) return 0;

	// Not altering the suffix of zeros of $[0..windowEnd-1]$.
	q=windowEnd-1;
	while (ms[q]==0) q--;
	remainingOnes=windowOnes;
	out[0]=0; out[1]=q; out[2]=q+1-remainingOnes; out[3]=remainingOnes; out[4]=initialMSValue;
	return 1;
}


/** 
 * Like $getPermutationRange_window()$, but for the case in which $j$ is the length of
 * $ms$.
 */
uint8_t getPermutationRange_lastWindow(sdsl::bit_vector &ms, const size_type windowStart, const size_type windowEnd, const size_type windowZeros, const size_type threshold, size_type out[]) {
	const size_type windowOnes = (windowEnd-windowStart-1)-windowZeros;
	size_type i, p;
	size_type delta, remainingOnes, initialMSValue;
	
	if (windowZeros==0) return 0;
	
	// Not altering the prefix of ones of $[windowStart+1..windowEnd-1]$.
	p=windowStart+1;
	while (ms[p]==1) p++;
	delta=p-windowStart-1;
	initialMSValue=threshold-delta;
	remainingOnes=windowOnes-delta;
	out[0]=p; out[1]=windowEnd-1; out[2]=windowZeros; out[3]=remainingOnes; out[4]=initialMSValue;
	return 1;
}


void loadTables(const string directory, const uint32_t threshold, const size_type nZeros, const size_type nOnes, const uint8_t allowNegativeMS) {
	uint32_t initialMSValue;
	size_t i, result, nNodes;
	int64_t *intBuffer;
	char *charBuffer = (char *)malloc(100*sizeof(char));
	FILE *file;
	
	printf("Loading tables... ");
	maxZ=nZeros; maxO=nOnes;
	neighbor=(int64_t **)malloc(threshold*sizeof(int64_t *));
	for (initialMSValue=0; initialMSValue<threshold; initialMSValue++) {
		neighbor[initialMSValue]=0;
		sprintf(charBuffer,"%s/table-diffLong-%d-%d-%lu-%lu-%d",directory.c_str(),threshold,initialMSValue,nZeros,nOnes,allowNegativeMS);
		file=fopen(charBuffer,"r");
		fread(&nNodes,sizeof(uint64_t),1,file);
		intBuffer=(int64_t *)malloc(nNodes*sizeof(int64_t));
		result=fread(intBuffer,sizeof(int64_t),nNodes,file);
		fclose(file);
		if (result!=nNodes) {
			cerr << "loadTables> Could not read nNodes=" << nNodes << " elements from file, just " << result << ".\n";
			exit(1);
		}
		neighbor[initialMSValue]=intBuffer;
	}
	printf("done \n");
	free(charBuffer);
}


/**
 * Stores in $out[0..1]$ a valid in-neighbor of point $(nZeros,nOnes)$, that maximizes
 * the ratio between uncompressed and compressed size of the added bits (if $strategy=0$),
 * or the total number of added bits (if $strategy=1$).
 * $out[0..1]$ is set to $(nZeros,nOnes)$ if no valid in-neighbor exists.
 */
void greedyInNeighbor(size_type nZeros, size_type nOnes, const size_type maxMSValue, size_type initialMSValue, const uint8_t allowNegativeMS, const uint8_t strategy) {
	size_type maxI, maxJ;
	int64_t i, j;
	int64_t fromI, toI, fromJ, msValue;
	double ratio, maxRatio;
	
	fromJ=initialMSValue+nZeros-maxMSValue-1;  // This is because we want:
	//
	// MS(i,j)+(nZeros-i)-1 <= maxMSValue              // first one of pair (nZeros,nOnes)
	//     
	// initialMSValue-j+i+nZeros-i-1 <= maxMSValue
	// j >= initialMSValue+nZeros-maxMSValue-1
	//
	if (fromJ<1) fromJ=1;
	maxRatio=0.0; maxI=nZeros; maxJ=nOnes;
	for (j=fromJ; j<=nOnes-1; j++) {
		msValue=initialMSValue-j;
		if (allowNegativeMS) fromI=1;
		else {
			fromI=j-initialMSValue;
			if (fromI<1) fromI=1;
		}
		if (fromI>=nZeros) break;
		toI=maxMSValue-msValue;  // This is because we want:
		//
		// initialMSValue-j+i <= maxMSValue            // last one of pair (i,j)
		//
		if (toI>=nZeros) toI=nZeros-1;
		if (strategy==0) {
			for (i=fromI; i<=toI; i++) {
				ratio=nZeros-i+nOnes-j;
				ratio/=encodingSize(nZeros-i,nOnes-j);
				if (ratio>maxRatio) {
					maxRatio=ratio;
					maxI=i; maxJ=j;
				}
			}
		}
		else {
			ratio=nZeros-fromI+nOnes-j;
			if (ratio>maxRatio) {
				maxRatio=ratio;
				maxI=fromI; maxJ=j;
			}
		}
	}
	out[0]=maxI; out[1]=maxJ;
}


/**
 * Transforms a window into a sequence of run-pairs $(nZeros_i,nOnes_i)$, using the 
 * solution stored in the $initialMSValue$-th DAG.
 *
 * @param start,end first and last position of the window to be permuted;
 * @param N_ZEROS,N_ONES total number of zeros and ones in the window;
 * @param permutedMS not assumed to be initialized in any way;
 * @param initialMSValue MS value of the last one-bit before $start$;
 * @param maxMSValue max MS value in a permuted range.
 */
void permuteWindow_impl(const size_type start, const size_type end, const size_type N_ZEROS, const size_type N_ONES, sdsl::bit_vector &permutedMS, const int64_t initialMSValue, size_type maxMSValue, uint8_t allowNegativeMS, uint8_t greedyStrategy) {
	int64_t i, j;
	size_type nZeros, nOnes;
	int64_t node, previousNode;
	
	if (N_ZEROS<=maxZ && N_ONES<=maxO) {
		// Only lookups
		nZeros=N_ZEROS; nOnes=N_ONES; i=end;
		node=coordinates2id(nZeros,nOnes,maxZ,maxO);
		while (node!=-1) {
			previousNode=neighbor[initialMSValue][node];
			if (previousNode!=-1) {
				id2coordinates(previousNode,maxZ,maxO,out);
				for (j=0; j<nOnes-out[1]; j++) permutedMS[i--]=1;
				for (j=0; j<nZeros-out[0]; j++) permutedMS[i--]=0;
				nZeros=out[0]; nOnes=out[1];
			}
			else {
				for (j=0; j<nOnes; j++) permutedMS[i--]=1;
				for (j=0; j<nZeros; j++) permutedMS[i--]=0;
			}
			node=previousNode;
		}
	}
	else {
		// Greedy + lookups
		nZeros=N_ZEROS; nOnes=N_ONES; i=end;
		while (true) {
			if (nZeros<=maxZ && nOnes<=maxO) {
				node=coordinates2id(nZeros,nOnes,maxZ,maxO);
				previousNode=neighbor[initialMSValue][node];
				if (previousNode!=-1) {
					id2coordinates(previousNode,maxZ,maxO,out);
					for (j=0; j<nOnes-out[1]; j++) permutedMS[i--]=1;
					for (j=0; j<nZeros-out[0]; j++) permutedMS[i--]=0;
					nZeros=out[0]; nOnes=out[1];
				}
				else {
					for (j=0; j<nOnes; j++) permutedMS[i--]=1;
					for (j=0; j<nZeros; j++) permutedMS[i--]=0;
					break;
				}
			}
			else {
				greedyInNeighbor(nZeros,nOnes,maxMSValue,initialMSValue,allowNegativeMS,greedyStrategy);
				if (out[0]!=nZeros || out[1]!=nOnes) {
					for (j=0; j<nOnes-out[1]; j++) permutedMS[i--]=1;
					for (j=0; j<nZeros-out[0]; j++) permutedMS[i--]=0;
					nZeros=out[0]; nOnes=out[1];
				}
				else {
					for (j=0; j<nOnes; j++) permutedMS[i--]=1;
					for (j=0; j<nZeros; j++) permutedMS[i--]=0;
					break;
				}
			}
		}
	}
}


/**
 * Permutes the bits in $ms[windowStart+1..windowEnd-1]$, where the $window*$-th one-bit 
 * of $ms$ has an MS value at least $threshold$, and all intermediate one-bits have an MS 
 * value smaller than $threshold$. The permutation is stored in 
 * $permutedMS[windowStart+1..windowEnd-1]$.
 *
 * Remark: the MS value of the $windowStart$-th bit of $ms$ must be equal to $threshold$. 
 * 
 * @param windowOnes number of ones in $[windowStart+1..windowEnd-1]$.
 */
void permuteWindow(sdsl::bit_vector &ms, sdsl::bit_vector &permutedMS, const size_type windowStart, const size_type windowEnd, const size_type windowOnes, const size_type threshold, uint8_t allowNegativeMS, uint8_t greedyStrategy) {
	uint8_t action;
	size_type i;
	
	action=getPermutationRange_window(ms,windowStart,windowEnd,windowOnes,threshold,out);
	for (i=windowStart+1; i<out[0]; i++) permutedMS[i]=1;
	if (action==0) {
		for (i=out[0]; i<windowEnd; i++) permutedMS[i]=0;
		return;
	}
	for (i=out[1]+1; i<windowEnd; i++) permutedMS[i]=0;
	permuteWindow_impl(out[0],out[1],out[2],out[3],permutedMS,out[4],threshold-1,allowNegativeMS,greedyStrategy);
}


/**
 * Remark: the first window of $ms$ contains fewer ones than zeros.
 */
void permuteFirstWindow(sdsl::bit_vector &ms, sdsl::bit_vector &permutedMS, const size_type windowEnd, const size_type windowOnes, const size_type threshold, uint8_t allowNegativeMS, uint8_t greedyStrategy) {
	const size_type initialMSValue = 1;
	uint8_t action;
	size_type i;
	
	action=getPermutationRange_firstWindow(ms,windowEnd,windowOnes,out);
	if (action==0) {
		for (i=0; i<windowEnd; i++) permutedMS[i]=0;
		return;
	}
	for (i=out[1]+1; i<windowEnd; i++) permutedMS[i]=0;
	permuteWindow_impl(0,out[1],out[2],out[3],permutedMS,1,threshold-1,allowNegativeMS,greedyStrategy);
}


/**
 * Remark: the last window of $ms$ contains fewer zeros than ones.
 */
void permuteLastWindow(sdsl::bit_vector &ms, sdsl::bit_vector &permutedMS, const size_type windowStart, const size_type windowEnd, const size_type windowZeros, const size_type threshold, uint8_t allowNegativeMS, uint8_t greedyStrategy) {
	uint8_t action;
	size_type i;
	
	action=getPermutationRange_lastWindow(ms,windowStart,windowEnd,windowZeros,threshold,out);
	if (action==0) {
		for (i=windowStart+1; i<windowEnd; i++) permutedMS[i]=1;
		return;
	}
	for (i=windowStart+1; i<out[0]; i++) permutedMS[i]=1;
	permuteWindow_impl(out[0],windowEnd-1,out[2],out[3],permutedMS,out[4],threshold-1,allowNegativeMS,greedyStrategy);
}


/**
 * @return 0 iff the number of ones is the same before and after $permuteWindow()$.
 */
uint8_t checkPermutation(sdsl::bit_vector &ms, sdsl::bit_vector &permutedMS, const int64_t windowStart, const size_type windowEnd) {
	size_type i, nOnesMS, nOnesPermutedMS;
	
	nOnesMS=0;
	for (i=windowStart+1; i<=windowEnd-1; i++) {
		if (ms[i]==1) nOnesMS++;
	}
	nOnesPermutedMS=0;
	for (i=windowStart+1; i<=windowEnd-1; i++) {
		if (permutedMS[i]==1) nOnesPermutedMS++;
	}
	if (nOnesPermutedMS!=nOnesMS) {
		cerr << "permuteWindow> ERROR: after the permutation, the number of one-bits is different: " << nOnesMS << " != " << nOnesPermutedMS << endl;
		return 1;
	}
	return 0;
}


/**
 * @return the number of bits in the window $[x..y]$ that differ between $ms$ and 
 * $permutedMS$.
 */
size_type windowXOR(sdsl::bit_vector &ms, sdsl::bit_vector &permutedMS, const size_type x, const size_type y) {
	size_type i, out;
	
	out=0;
	for (i=x; i<=y; i++) {
		if (ms[i]!=permutedMS[i]) out++;
	}
	return out;
}


/**
 * @param ms the procedure assumes that at least one MS value is at least $threshold$;
 * @param permutedMS initialized to all zeros; at the end of the procedure, it contains
 * a permutation of $ms$ where all the ones of $ms$ with MS value at least $threshold$ 
 * have the same MS value, and all other ones are assigned a (possibly different) MS value
 * smaller than $threshold$.
 */
void permuteBitvector(sdsl::bit_vector &ms, sdsl::bit_vector &permutedMS, const size_type threshold, const uint8_t allowNegativeMS, const uint8_t greedyStrategy, const uint8_t verbose) {
	const size_type MS_SIZE = ms.size();
	uint8_t error;
	size_type i, j;
	size_type msValue, msValueStart, nZeros, windowZeros, windowOnes, windowSize;
	size_type nOnes, nChangedBits, lastInitialMSValue;
	int64_t windowStart;
	
	cerr << "Permuting... ";
	msValue=1;
	nZeros=0;  // Number of zeros between two consecutive ones
	windowStart=-1;  // First bit whose MS is $>=threshold$.
	msValueStart=0;
	windowOnes=0;  // Number of ones in $[windowStart+1..i-1]$.
	nOnes=0; nChangedBits=0;
	for (i=0; i<MS_SIZE; i++) {
		if (ms[i]==0) {
			nZeros++;
			continue;
		}
		nOnes++;
		msValue+=nZeros-1;
		if (msValue>=threshold) {
			permutedMS[i]=1;
			// Current window
			windowSize=i-windowStart-1;
			if (windowSize>1 && windowOnes>0) {
				if (windowStart!=-1 && msValueStart!=threshold) {
					cerr << "ERROR: msValueStart=" << msValueStart << " != threshold=" << threshold << endl;
					exit(1);
				}
				if (windowOnes==windowSize) {
					cerr << "ERROR: the window contains only ones?! [" << (windowStart+1) << ".." << (i-1) << "]" << endl;
					exit(1);
				}
				if (verbose) {
					cerr << "Permuting interval [" << (windowStart+1) << ".." << (i-1) << "]  length=" << windowSize << " windowOnes=" << windowOnes << " (" << (100*(((double)i)/MS_SIZE)) << "%)  MS_SIZE=" << MS_SIZE << endl;
					cerr << "MS window: ";
					for (size_type x=windowStart+1; x<i; x++) cerr << ms[x] << ",";
					cerr << endl;
				}
				if (windowStart==-1) permuteFirstWindow(ms,permutedMS,i,windowOnes,threshold,allowNegativeMS,greedyStrategy);
				else permuteWindow(ms,permutedMS,windowStart,i,windowOnes,threshold,allowNegativeMS,greedyStrategy);
				if (verbose) {
					cerr << "Permuted window: ";
					for (size_type x=windowStart+1; x<i; x++) cerr << permutedMS[x] << ",";
					cerr << endl;
				}
				error=checkPermutation(ms,permutedMS,windowStart,i);
				if (error) {
					cerr << "Permuting interval [" << (windowStart+1) << ".." << (i-1) << "]  length=" << windowSize << " windowOnes=" << windowOnes << " (" << (100*(((double)i)/MS_SIZE)) << "%)  MS_SIZE=" << MS_SIZE << endl;
					cerr << "MS window: ";
					for (size_type x=windowStart+1; x<i; x++) cerr << ms[x] << ",";
					cerr << endl;
					cerr << "Permuted window: ";
					for (size_type x=windowStart+1; x<i; x++) cerr << permutedMS[x] << ",";
					cerr << endl;
					exit(1);
				}
				nChangedBits+=windowXOR(ms,permutedMS,windowStart+1,i-1);
			}
			// Next window
			windowStart=i; msValueStart=msValue; windowOnes=0;
		}
		else windowOnes++;
		nZeros=0;
	}

	if (windowStart==-1) {
		// No large-enough MS value: permuting the whole bitvector.
		if (verbose) cerr << "Permuting the whole bitvector... length=" << MS_SIZE << " windowOnes=" << windowOnes << endl;
		permuteWindow_impl(0,MS_SIZE-1,MS_SIZE-windowOnes,windowOnes,permutedMS,1,threshold-1,allowNegativeMS,greedyStrategy);
		error=checkPermutation(ms,permutedMS,-1,MS_SIZE);
		if (error) {
			cerr << "MS window: ";
			for (size_type x=0; x<MS_SIZE; x++) cerr << ms[x] << ",";
			cerr << endl;
			cerr << "Permuted window: ";
			for (size_type x=0; x<MS_SIZE; x++) cerr << permutedMS[x] << ",";
			cerr << endl;
			exit(1);
		}
		nChangedBits+=windowXOR(ms,permutedMS,0,MS_SIZE-1);	
	}
	else {
		// Permuting the last window
		windowSize=MS_SIZE-windowStart-1;
		if (windowSize>1 && windowOnes>0) {
			if (windowStart!=-1 && msValueStart!=threshold) {
				cerr << "ERROR: msValueStart=" << msValueStart << " != threshold=" << threshold << endl;
				exit(1);
			}
			if (windowOnes==windowSize) {  // The last window can contain only ones
				for (i=0; i<windowOnes; i++) permutedMS[windowStart+1+i]=1;
			}
			else {		
				if (verbose) cerr << "Permuting interval [" << (windowStart+1) << ".." << (MS_SIZE-1) << "]  length=" << windowSize << " windowOnes=" << windowOnes << " (100%)" << endl;
				permuteLastWindow(ms,permutedMS,windowStart,MS_SIZE,windowSize-windowOnes,threshold,allowNegativeMS,greedyStrategy);
				error=checkPermutation(ms,permutedMS,windowStart,MS_SIZE);
				if (error) {
					cerr << "Permuting interval [" << (windowStart+1) << ".." << (MS_SIZE-1) << "]  length=" << windowSize << " windowOnes=" << windowOnes << endl;
					cerr << "MS window: ";
					for (size_type x=windowStart+1; x<MS_SIZE; x++) cerr << ms[x] << ",";
					cerr << endl;
					cerr << "Permuted window: ";
					for (size_type x=windowStart+1; x<MS_SIZE; x++) cerr << permutedMS[x] << ",";
					cerr << endl;
					exit(1);
				}
				nChangedBits+=windowXOR(ms,permutedMS,windowStart+1,MS_SIZE-1);
			}
		}
	}
	cerr << "done \n";
	if (verbose) cerr << "Changed " << nChangedBits << " bits in total (" << ((100.0*nChangedBits)/MS_SIZE) << "%)" << endl;
}


/**
 * Checks that every one-bit in $ms$ with an MS value at least $threshold$, has a 
 * corresponding one-bit in $permutedMS$ with identical MS value.
 */
void equalMSValues(sdsl::bit_vector &ms, sdsl::bit_vector &permutedMS, const size_type threshold, const uint8_t verbose) {
	const size_type MS_SIZE = ms.size();
	const int64_t THRESHOLD = (int64_t)threshold;
	size_type i, j;
	size_type nZeros, permutedNZeros;
	int64_t msValue, permutedMSValue;
	
	if (verbose) {
		cerr << "Checking that the two bitvectors have the same number of bits..." << endl;
	}
	if (permutedMS.size()!=MS_SIZE) {
		cerr << "ERROR: different size: " << MS_SIZE << " != " << permutedMS.size() << endl;
		exit(1);
	}
	nZeros=0;
	for (i=0; i<MS_SIZE; i++) {
		if (ms[i]==0) nZeros++; 
	}
	permutedNZeros=0;
	for (i=0; i<MS_SIZE; i++) {
		if (permutedMS[i]==0) permutedNZeros++; 
	}
	if (permutedNZeros!=nZeros) {
		cerr << "ERROR: different number of zeros: " << nZeros << " != " << permutedNZeros << endl;
		exit(1);
	}
	if (verbose) cerr << "OK" << endl;
	if (verbose) cerr << "Comparing the MS values of corresponding one-bits..." << endl;
	msValue=1; permutedMSValue=1;
	nZeros=0; permutedNZeros=0;
	i=0; j=0; 
	while (i<MS_SIZE && j<MS_SIZE) {
		if (ms[i]==0) {
			nZeros++;
			i++;
			continue;
		}
		if (permutedMS[j]==0) {
			permutedNZeros++;
			j++;
			continue;
		}
		msValue+=nZeros-1;
		permutedMSValue+=permutedNZeros-1;
		if (msValue>=THRESHOLD) {
			if (permutedMSValue!=msValue) {
				cerr << "ERROR at coordinates (" << i << "," << j << "): msValue=" << msValue << " != permutedMSValue=" << permutedMSValue << endl;
				exit(1);
			}
		}
		else {
			if (permutedMSValue>=THRESHOLD) {
				cerr << "ERROR at coordinates (" << i << "," << j << "): msValue=" << msValue << " is below threshold, but permutedMSValue=" << permutedMSValue << " is above threshold" << endl;
				exit(1);
			}
		}
		nZeros=0; i++;
		permutedNZeros=0; j++;
	}
	if (verbose) cerr << "OK" << endl;
}


/**
 * @return the max MS value encoded in $ms$.
 */
int64_t maxMSValue(sdsl::bit_vector &ms) {
	const size_type MS_SIZE = ms.size();
	size_type i;
	size_type nZeros;
	int64_t msValue, maxValue;
	
	maxValue=0; msValue=1; nZeros=0;
	for (i=0; i<MS_SIZE; i++) {
		if (ms[i]==0) {
			nZeros++;
			continue;
		}
		msValue+=nZeros-1;
		if (msValue>maxValue) maxValue=msValue;
		nZeros=0;
	}
	return maxValue;
}


/**
 * Prints all pairs $(nZeros,nOnes)$.
 */
void runPairsStats(sdsl::bit_vector &permutedMS) {
	const size_type MS_SIZE = permutedMS.size();
	size_type i;
	size_type nZeros, nOnes;
	
	nZeros=1; nOnes=0;
	for (i=1; i<MS_SIZE; i++) {
		if (permutedMS[i]==0 && permutedMS[i-1]==1) {
			cout << nZeros << "," << nOnes << endl;
			nZeros=1; nOnes=0;
			continue;
		}
		if (permutedMS[i]==0) nZeros++;
		else nOnes++;
	}
	cout << nZeros << "," << nOnes << endl;
}


// ------------------------------- END OF PERMUTATION ------------------------------------



template<typename vec_type = CSA::RLEVector, typename enc_type = CSA::RLEEncoder>
size_type comp1(const string ms_path, const InputFlags& flags, const string outputDirectory, const string tablesDirectory) {
	uint32_t initialMSValue;
	sdsl::bit_vector ms;
    sdsl::load_from_file(ms,ms_path);
	sdsl::bit_vector permutedMS(ms.size(),0);
	uint32_t threshold = flags.threshold;
	size_type nZeros = flags.n_zeros;
	size_type nOnes = flags.n_ones;
	uint8_t allowNegativeMS = (uint8_t) flags.negative;
	uint8_t greedyStrategy = flags.greedy;
	uint8_t verbose = flags.verbose;

	
	
	// Permuting
	loadTables(tablesDirectory,threshold,nZeros,nOnes,allowNegativeMS);
	permuteBitvector(ms,permutedMS,threshold,allowNegativeMS,greedyStrategy,verbose);
	cerr << "Permuting completed" << endl;
	for (initialMSValue=0; initialMSValue<threshold; initialMSValue++) free(neighbor[initialMSValue]);
	free(neighbor);
	equalMSValues(ms,permutedMS,threshold,verbose);

    size_type from = abs_point();
    enc_type encoder(32);
    histo_t counter;
    size_type n_runs = fill_encoder<enc_type>(permutedMS, encoder, counter);
    vec_type c_ms(encoder, permutedMS.size());	
	char buffer[100]; 
	sprintf(buffer,"%s/compressedMS-diffLong-%d-%d-%d",outputDirectory.c_str(),threshold,allowNegativeMS,greedyStrategy);
    std::ofstream out{  buffer   /*ms_path + ms_compression::to_str(flags.compression)*/, std::ios::binary};
    c_ms.writeTo(out);

    (cerr << n_runs << " runs over "
          << c_ms.getSize() << " elements ("
          << c_ms.getSize() / static_cast<float>(n_runs) << " elements / run)"
          << endl);
    //for(auto item : counter)
    //    cout << item.first << "," << item.second << endl;
    return diff_from(from);
}


int main(int argc, char **argv){
    OptParser input(argc, argv);
    string ms_path;

    if(argc == 1){
        (cerr << "Compress a ms vector. Creates files <ms_path>.xxx\n"
              << "Args:\n"
              << help__ms_path
              << help__ms_path
              << help__threshold
              << help__nzeros
              << help__nones
              << help__negative
              << help__greddy
              << help__verbose
              << endl);
        exit(0);
    }
    InputFlags flags;
    try{
        flags = InputFlags(input);
    }
    catch (string s) {
        cerr << s << endl;
        return 1;
    }
    ms_path = input.getCmdOption("-ms_path");
    cerr << "check: " << diff_from(abs_point()) << endl;//size_type mem_mark = abs_point();
    cout << comp1<>(ms_path,flags, input.getCmdOption("-outputDir"), input.getCmdOption("-tablesDir"))
         << endl;
    return 0;
}
