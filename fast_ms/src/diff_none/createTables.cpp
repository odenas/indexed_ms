/**
 * Let $T$ be a threshold and let $[i..j]$ be an interval of the MS array such that  
 * $MS[i]>=T, MS[j]>=T$, and $MS[k]<T$ for all $k \in [i+1..j-1]$. The program builds a 
 * DAG that encodes how to permute $[i+1..j]$ into a sequence of run-pairs $(nZeros_0,
 * nOnes_0),...,(nZeros_x,nOnes_x),...$, where $nZeros_0$ might be zero and every other  
 * $nZeros_x$ and $nOnes_x$ is positive, that minimizes the size of the encoding, where 
 * each $nZeros_x$ and $nOnes_x$ is delta-coded separately. The minimal encoding is 
 * equivalent to a shortest path in the DAG: the program computes all shortest paths in 
 * order to answer multiple queries later.
 */
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>


/** 
 * Representation of the DAG
 */
uint8_t ALLOW_NEGATIVE_MS;
uint32_t THRESHOLD;
uint64_t N_ZEROS, N_ONES, N_NODES, N_ARCS;
uint64_t **inNeighbors;
uint64_t *rowLengths;
int64_t *lastInNeighbor;
uint64_t *cost;
int64_t *neighbor;


void allocateMemory() {
	const uint64_t CAPACITY = 10;  // Arbitrary
	uint64_t i;
	
	inNeighbors=(uint64_t **)malloc(N_NODES*sizeof(uint64_t *));
	for (i=0; i<N_NODES; i++) inNeighbors[i]=(uint64_t *)malloc(CAPACITY*sizeof(uint64_t *));
	rowLengths=(uint64_t *)malloc(N_NODES*sizeof(uint64_t));
	for (i=0; i<N_NODES; i++) rowLengths[i]=CAPACITY;
	lastInNeighbor=(int64_t *)malloc(N_NODES*sizeof(uint64_t));
	cost=(uint64_t *)malloc(N_NODES*sizeof(uint64_t));
	neighbor=(int64_t *)malloc(N_NODES*sizeof(int64_t));
	for (i=0; i<N_NODES; i++) lastInNeighbor[i]=-1;
}


void deallocateMemory() {
	uint64_t i;
	
	for (i=0; i<N_NODES; i++) free(inNeighbors[i]);
	free(rowLengths);
	free(lastInNeighbor);
	free(cost);
	free(neighbor);
}


void clean() {
	uint64_t i;
	
	for (i=0; i<N_NODES; i++) lastInNeighbor[i]=-1;
	for (i=0; i<N_NODES; i++) cost[i]=ULLONG_MAX;
	for (i=0; i<N_NODES; i++) neighbor[i]=-1;
}


/**
 * @param nZeros,nOnes >=1.
 */
uint64_t coordinates2id(const uint64_t nZeros, const uint64_t nOnes, const uint64_t totalNZeros, const uint64_t totalNOnes) {
	return nZeros*totalNOnes+(nOnes-1);
}


void id2coordinates(const uint64_t nodeID, const uint64_t totalNZeros, const uint64_t totalNOnes, uint64_t out[]) {
	out[0]=nodeID/totalNOnes;
	out[1]=1+nodeID%totalNOnes;
}


void addArc(const uint64_t nZerosFrom, const uint64_t nOnesFrom, const uint64_t nZerosTo, const uint64_t nOnesTo, const uint64_t totalNZeros, const uint64_t totalNOnes) {
	const uint64_t fromNode = coordinates2id(nZerosFrom,nOnesFrom,totalNZeros,totalNOnes);
	const uint64_t toNode = coordinates2id(nZerosTo,nOnesTo,totalNZeros,totalNOnes);
	
	lastInNeighbor[toNode]++;
	if (lastInNeighbor[toNode]==rowLengths[toNode]) {
		inNeighbors[toNode]=(uint64_t *)realloc(inNeighbors[toNode],(rowLengths[toNode]<<1)*sizeof(uint64_t));
		rowLengths[toNode]<<=1;
	}
	inNeighbors[toNode][lastInNeighbor[toNode]]=fromNode;
	N_ARCS++;
}


uint64_t length(uint64_t n) {
  uint64_t b = 0;
  while (n!=0) { b++; n>>=1; }
  return b;
}


uint64_t deltaCodeLength(const uint64_t value) {
	uint64_t len = length(value);
	uint64_t llen = length(len);
	return len+llen+llen-2;
}


/**
 * @return number of bits for encoding the run pair $(nZeros,nOnes)$.
 *
 * @param $nZeros$ might be zero.
 */
uint64_t encodingSize(const uint64_t nZeros, const uint64_t nOnes) {
	return (nZeros==0?0:deltaCodeLength(nZeros))+deltaCodeLength(nOnes);
}


void buildDAG() {
	const int32_t initialMSValue = THRESHOLD;
	const int32_t maxMSValue = THRESHOLD-1;  // Max MS value in a permuted range
	uint8_t canBeFirst;
	int64_t i, j, p, q;
	uint64_t nZeros, nOnes;
	uint64_t fromNode, toNode, minCost, newCost;
	int64_t ms, msValue, minNeighbor, node, previousNode, fromI, toI, toJ, fromP, toP;
	uint64_t out[5];
	
	// Building the DAG
	printf("Building DAG for initialMSValue=%d, N_NODES=%lu, N_ZEROS=%lu, N_ONES=%lu... ",initialMSValue,N_NODES,N_ZEROS,N_ONES);
	// First column ($nZeros=0$).
	if (ALLOW_NEGATIVE_MS) toJ=N_ONES-1;
	else {
		toJ=initialMSValue;
		if (toJ>N_ONES-1) toJ=N_ONES-1;
	}
	for (j=1; j<=toJ; j++) {
		toP=maxMSValue-initialMSValue+j+1;  // This is because we want:
		//
		// (1) initialMSValue-q+p <= maxMSValue       // last one-bit of pair (p,q)
		// (2) initialMSValue-j+p-1 <= maxMSValue     // first one-bit of pair (p,q)
		//     
		// (1) p <= maxMSValue-initialMSValue+q
		// (2) p <= maxMSValue-initialMSValue+j+1
		//
		if (toP>N_ZEROS) toP=N_ZEROS;
		for (q=j+1; q<=N_ONES; q++) {
			if (ALLOW_NEGATIVE_MS) fromP=1;
			else {
				fromP=q-initialMSValue;
				if (fromP<1) fromP=1;
				else if (fromP>N_ZEROS) break;
			}
			for (p=fromP; p<=toP; p++) addArc(0,j,p,q,N_ZEROS,N_ONES);
		}
	}
	// Other columns
	for (j=1; j<=N_ONES-1; j++) {
		ms=initialMSValue-j;
		if (ALLOW_NEGATIVE_MS) fromI=1;
		else {
			fromI=j-initialMSValue;
			if (fromI<1) fromI=1;
			else if (fromI>N_ZEROS) break;
		}
		toI=maxMSValue-ms;
		if (toI>N_ZEROS) toI=N_ZEROS;
		toP=toI+1;  // This is because we want:
		//
		// (1) initialMSValue-q+p <= maxMSValue     // last one-bit of pair (p,q)
		// (2) MS(i,j)+(p-i)-1 <= maxMSValue        // first one-bit of pair (p,q)
		//     
		// (1) p <= maxMSValue-initialMSValue+q
		// (2) p <= maxMSValue-(initialMSValue-j+i)+i+1
		//     p <= maxMSValue-initialMSValue+j+1
		//
		if (toP>N_ZEROS) toP=N_ZEROS;
		for (i=fromI; i<=toI; i++) {
			for (q=j+1; q<=N_ONES; q++) {
				if (ALLOW_NEGATIVE_MS) fromP=i+1;
				else {
					fromP=q-initialMSValue;
					if (fromP<i+1) fromP=i+1;
					else if (fromP>N_ZEROS) break;
				}
				for (p=fromP; p<=toP; p++) {
					addArc(i,j,p,q,N_ZEROS,N_ONES);
				}
			}
		}
	}
	printf("done   N_ARCS=%lu \n",N_ARCS);
	
	// Computing all single-source shortest paths in the DAG
	printf("Computing shortest paths... ");
	for (i=0; i<N_NODES; i++) {
		id2coordinates(i,N_ZEROS,N_ONES,out);
		nZeros=out[0]; nOnes=out[1];
		canBeFirst=1;
		if (nZeros==0) {
			if (!ALLOW_NEGATIVE_MS && initialMSValue<nOnes) canBeFirst=0;
		}
		else {
			if (initialMSValue+nZeros>maxMSValue+1) canBeFirst=0;
			else {
				if (!ALLOW_NEGATIVE_MS && initialMSValue+nZeros<nOnes) canBeFirst=0;
			}
		}
		if (canBeFirst) minCost=encodingSize(nZeros,nOnes);
		else minCost=ULLONG_MAX;
		minNeighbor=-1;
		for (j=0; j<=lastInNeighbor[i]; j++) {
			fromNode=inNeighbors[i][j];
			if (cost[fromNode]==ULLONG_MAX) continue;
			id2coordinates(fromNode,N_ZEROS,N_ONES,out);
			newCost=cost[fromNode]+encodingSize(nZeros-out[0],nOnes-out[1]);
			if (newCost<minCost) {
				minCost=newCost;
				minNeighbor=fromNode;
			}					
		}
		cost[i]=minCost; neighbor[i]=minNeighbor;
	}
	printf("done \n");
	
	// The graph will be queried for nodes $(nZeros,nOnes)$ such that the MS value of the
	// last one-bit is $>=threshold$, and such nodes have no in-neighbor up to now.
	// Connecting every such node $(nZeros,nOnes)$ to all valid nodes with $nOnes-1$
	// one-bits, disregarding in every connection the constraint on the first one-bit.
	printf("Computing cost of query nodes... ");
	for (i=N_NODES-1; i>=0; i--) {
		id2coordinates(i,N_ZEROS,N_ONES,out);
		nZeros=out[0]; nOnes=out[1];
		if (nZeros==0) break;
		if (nOnes<=2 || initialMSValue+nZeros<=nOnes+maxMSValue) continue;
		if (lastInNeighbor[i]>=0) {
			printf("buildDAG_impl> ERROR: a window with large MS value has already %ld in-neighbors?!  (%lu,%lu) \n",(lastInNeighbor[i]+1),nZeros,nOnes);
			exit(1);
		}
		if (cost[i]!=ULLONG_MAX) {
			printf("buildDAG_impl> ERROR: a window with large MS value has already a finite cost?! %lu (%lu,%lu) \n",cost[i],nZeros,nOnes);
			exit(1);
		}
		ms=nOnes-initialMSValue-1;
		//
		// We want the in-neighbor $(p,nOnes-1)$ to satisfy:
		//
		// (1) p <= nZeros
		// (2) initialMSValue+p-(nOnes-1) <= maxMSValue
		// (3) initialMSValue+p-(nOnes-1) >=0                   // if ALLOW_NEGATIVE_MS==0
		//
		// (1,2) p <= min{ nOnes-initialMSValue-1+maxMSValue, nZeros }
		//   (3) p >= nOnes-initialMSValue-1                    // if ALLOW_NEGATIVE_MS==0
		//
		fromP=ALLOW_NEGATIVE_MS?0:(ms>=0?ms:0);
		toP=ms+maxMSValue;
		if (toP>nZeros) toP=nZeros;
		minCost=ULLONG_MAX; minNeighbor=-1;
		for (p=fromP; p<=toP; p++) {
			fromNode=coordinates2id(p,nOnes-1,N_ZEROS,N_ONES);
			if (cost[fromNode]==ULLONG_MAX) {
				printf("buildDAG_impl> ERROR: a previous valid node has infinite cost?! (%lu,%lu) \n",p,nOnes-1);
				exit(1);
			}
			newCost=cost[fromNode]+encodingSize(nZeros-p,1);
			if (newCost<minCost) {
				minCost=newCost;
				minNeighbor=fromNode;
			}
		}
		cost[i]=minCost; neighbor[i]=minNeighbor;
	}
	printf("done \n");
}


/**
 * @param argv in the following order: 
 * threshold_first
 * threshold_last
 * maxNZeros
 * maxNOnes
 * allowNegativeMS 
 * saveAll
 * destinationDir
 */
int main(int argc, char **argv){
	const uint32_t THRESHOLD_FIRST = (uint32_t)atoi(argv[1]);
	const uint32_t THRESHOLD_LAST = (uint32_t)atoi(argv[2]);
	N_ZEROS=(uint64_t)atol(argv[3]);
	N_ONES=(uint64_t)atol(argv[4]);
	ALLOW_NEGATIVE_MS=(uint8_t)atoi(argv[5]);
	const uint8_t SAVE_ALL = (uint8_t)atoi(argv[6]);
	const char *DIRECTORY = argv[7];
	N_NODES=(N_ZEROS+1)*N_ONES;
	
	uint64_t j;
	FILE *file;
	char *buffer = (char *)malloc(100*sizeof(char));
	
	allocateMemory();
	for (THRESHOLD=THRESHOLD_FIRST; THRESHOLD<=THRESHOLD_LAST; THRESHOLD++) {
		N_ARCS=0;
		clean();
		buildDAG();
		sprintf(buffer,"%s/table-nodiff-%d-%lu-%lu-%d",DIRECTORY,THRESHOLD,N_ZEROS,N_ONES,ALLOW_NEGATIVE_MS);
		file=fopen(buffer,"w");
		fwrite(&N_NODES,sizeof(uint64_t),1,file);
		fwrite(neighbor,sizeof(int64_t),N_NODES,file);
		if (SAVE_ALL) {
			fwrite(cost,sizeof(uint64_t),N_NODES,file);
			fwrite(lastInNeighbor,sizeof(int64_t),N_NODES,file);
			for (j=0; j<N_NODES; j++) fwrite(inNeighbors[j],sizeof(uint64_t),lastInNeighbor[j]+1,file);
		}
		fclose(file);
	}
	free(buffer);
	deallocateMemory();
}
