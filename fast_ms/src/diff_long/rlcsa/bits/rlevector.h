#ifndef RLEVECTOR_H
#define RLEVECTOR_H

#include <fstream>

#include "bitvector.h"
#include <stdlib.h>

namespace CSA
{


/*
  This class is used to construct a RLEVector.
*/

class RLEEncoder : public VectorEncoder
{
  public:
    RLEEncoder(usint block_bytes, usint superblock_size = VectorEncoder::SUPERBLOCK_SIZE);
    ~RLEEncoder();

    void setBit(usint value);
    void setRun(usint start, usint len);

    // These versions try to combine the runs if possible.
    void addBit(usint value);
    void addRun(usint start, usint len);
    void flush(); // Call this when finished.
	
	
	
	
	
	
    inline void RLEncode(usint diff, usint len)  // FC> (distance from previous 1-bit, length of 1-run).
    {
      this->size += diff + len - 1;
      this->items += len;
      this->buffer->writeDeltaCode(diff);
//printf("%lu,%lu -> %lu \n",diff,len,project(diff,len,1,THRESHOLD));
      this->buffer->writeDeltaCode(project(diff,len,1,(usint)atoi(getenv("THRESHOLD"))));  // FC> Remark: $writeDeltaCode()$ assumes positive input.
    }


	/**
 	 * Encodes $nOnes$ as a difference with respect to $nZeros$, if $nZeros>threshold$.
 	 * If $nZeros<=threshold$, the procedure just returns $nOnes$.
     	 *
     	 * @author FC
 	 * @param nOnes >0;
 	 * @param base first integer >0 to be used as the destination of the projection.
  	 */
	inline usint project(usint nZeros, usint nOnes, usint base, usint threshold) {
		usint delta;
		
		if (nZeros<=threshold) return nOnes;
		if (nOnes==nZeros) return base;
		if (nOnes<nZeros) {
			delta=nZeros-nOnes;
			return base-1+(delta<<1);
		}
		delta=nOnes-nZeros;
		if (delta<=nZeros-1) return base+(delta<<1);
		return base+((nZeros-1)<<1)+(delta-nZeros+1);
	}
	
	
	/**
	 * Inverse of $project()$.
     	 *
     	 * @author FC
	 * @param projection output of $project()$.
 	 */
	inline usint antiproject(usint nZeros, usint projection, usint base, usint threshold) {
		usint max, delta;
		
		if (nZeros<=threshold) return projection;
		if (projection==base) return nZeros;
		max=base+((nZeros-1)<<1);
		if (projection>max) return nZeros+(nZeros-1+projection-max);
		delta=projection-base;
		if (delta%2==0) return nZeros+(delta>>1);
		return nZeros-((delta+1)>>1);
	}
	
	
	
	
	


  protected:
    pair_type run;  // FC> ????????????????????????

    // These are not allowed.
    RLEEncoder();
    RLEEncoder(const RLEEncoder&);
    RLEEncoder& operator = (const RLEEncoder&);
};



















/*
  This is a run-length encoded bitvector using delta coding.
*/

class RLEVector : public BitVector
{
  public:
    typedef RLEEncoder Encoder;

    explicit RLEVector(std::ifstream& file);
    explicit RLEVector(FILE* file);
    RLEVector(Encoder& encoder, usint universe_size);
    ~RLEVector();

//--------------------------------------------------------------------------

    usint reportSize() const;

//--------------------------------------------------------------------------

    class Iterator : public BitVector::Iterator
    {
      public:
        explicit Iterator(const RLEVector& par);
        ~Iterator();

        usint rank(usint value, bool at_least = false);  // FC> Remark: rank,select are properties of the iterator rather than of the bitvector.

        usint select(usint index);
        usint selectNext();

        pair_type valueBefore(usint value);
        pair_type valueAfter(usint value);
        pair_type nextValue();

        pair_type selectRun(usint index, usint max_length);
        pair_type selectNextRun(usint max_length);

        bool isSet(usint value);

        usint countRuns();

      protected:

        void valueLoop(usint value);

        // These are not allowed.
        Iterator();
        Iterator(const Iterator&);
        Iterator& operator = (const Iterator&);
    };

//--------------------------------------------------------------------------

  protected:

    // These are not allowed.
    RLEVector();
    RLEVector(const RLEVector&);
    RLEVector& operator = (const RLEVector&);
};


} // namespace CSA


#endif // RLEVECTOR_H
