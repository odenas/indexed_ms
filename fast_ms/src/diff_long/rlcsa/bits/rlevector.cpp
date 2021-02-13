#include <cstdlib>

#include "rlevector.h"
#include "../misc/utils.h"


namespace CSA
{


RLEVector::RLEVector(std::ifstream& file) :
  BitVector(file)
{
}

RLEVector::RLEVector(FILE* file) :
  BitVector(file)
{
}

RLEVector::RLEVector(Encoder& encoder, usint universe_size) :
  BitVector(encoder, universe_size)
{
}

RLEVector::~RLEVector()
{
}

//--------------------------------------------------------------------------

usint
RLEVector::reportSize() const
{
  usint bytes = sizeof(*this);
  bytes += BitVector::reportSize();
  return bytes;
}

//--------------------------------------------------------------------------

RLEVector::Iterator::Iterator(const RLEVector& par) :
  BitVector::Iterator(par)
{
}

RLEVector::Iterator::~Iterator()
{
}



















/** 
 * FC> Remark: $value$ is a position in the bitvector.
 * FC> This is used for rank.
 * FC> The procedure iterates just over one block.
 */
void
RLEVector::Iterator::valueLoop(usint value)  
{
  this->getSample(this->sampleForValue(value));
  this->run = 0;

  if(this->val >= value) { return; }
  usint delta, runLength;
  while(this->cur < this->block_items)  // FC> $cur$ iterates over the ones in the block
  {
    this->val += this->buffer.readDeltaCode();   // FC> distance from previous 1-bit
    this->cur++;
	
	// <FC>
	delta=this->buffer.readDeltaCode()-1;
	if (delta==0) runLength=this->val;
	else if (delta%2==0) runLength=this->val+(delta>>1);
	else runLength=this->val-(delta>>1)-1;
	this->run=runLength-1;  // FC> n. of other bits in the 1-run
	// </FC>
    
    if(this->val >= value) { break; }

    this->cur += this->run;
    this->val += this->run;
    if(this->val >= value)
    {
      this->run = this->val - value;  // FC> $run$ contains the excess ones to the right of position $value$.
      this->val = value;
      this->cur -= this->run;
      break;
    }
    this->run = 0;
  }
}

usint
RLEVector::Iterator::rank(usint value, bool at_least)
{
  const RLEVector& par = (const RLEVector&)(this->parent);

  if(value >= par.size) { return par.items; }

  this->valueLoop(value);

  usint idx = this->sample.first + this->cur + 1;
  if(!at_least && this->val > value)
  {
    idx--;
  }
  if(at_least && this->val < value)
  {
    this->getSample(this->block + 1);
    this->run = 0;
    idx = this->sample.first + this->cur + 1;
  }
  return idx;
}

/**
 * FC> $index$ is of course the rank of a 1-bit.
 */
usint
RLEVector::Iterator::select(usint index)
{
  const RLEVector& par = (const RLEVector&)(this->parent);

  if(index >= par.items) { return par.size; }
  this->getSample(this->sampleForIndex(index));
  this->run = 0;

  usint lim = index - this->sample.first;
  usint delta, runLength, temp;
  while(this->cur < lim)  // FC> $cur$ iterates over the ones in the block
  {
    this->val += this->buffer.readDeltaCode();
	
	// <FC>
	delta=this->buffer.readDeltaCode()-1;
	if (delta==0) runLength=this->val;
	else if (delta%2==0) runLength=this->val+(delta>>1);
	else runLength=this->val-(delta>>1)-1;
	temp=runLength;
	// </FC>
	
    this->val += temp - 1;
    this->cur += temp;
  }
  if(this->cur > lim)
  {
    this->run = this->cur - lim;
    this->cur -= this->run;
    this->val -= this->run;
  }

  return this->val;
}































// FC> -------------------------- BEGIN NOT STUDIED YET ----------------------------------
usint
RLEVector::Iterator::selectNext()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    this->run = 0;
    return this->val;
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer.readDeltaCode();
    this->run = this->buffer.readDeltaCode() - 1;
  }

  return this->val;
}

pair_type
RLEVector::Iterator::valueBefore(usint value)
{
  const RLEVector& par = (const RLEVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }

  this->getSample(this->sampleForValue(value));
  if(this->val > value) { return pair_type(par.size, par.items); }
  this->run = 0;

  while(this->cur < this->block_items && this->val < value)
  {
    usint temp = this->buffer.readDeltaCode(value - this->val);
    if(temp == 0) { break; }
    this->val += temp;

    temp = this->buffer.readDeltaCode();
    this->cur += temp;
    this->val += temp - 1;
  }
  if(this->val > value)
  {
    this->run = this->val - value;
    this->val = value;
    this->cur -= this->run;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
RLEVector::Iterator::valueAfter(usint value)
{
  const RLEVector& par = (const RLEVector&)(this->parent);

  if(value >= par.size) { return pair_type(par.size, par.items); }

  this->valueLoop(value);

  if(this->val < value)
  {
    this->getSample(this->block + 1);
    this->run = 0;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
RLEVector::Iterator::nextValue()
{
  if(this->cur >= this->block_items)
  {
    this->getSample(this->block + 1);
    this->run = 0;
    return pair_type(this->val, this->sample.first);
  }

  this->cur++;
  if(this->run > 0)
  {
    this->val++;
    this->run--;
  }
  else
  {
    this->val += this->buffer.readDeltaCode();
    this->run = this->buffer.readDeltaCode() - 1;
  }

  return pair_type(this->val, this->sample.first + this->cur);
}

pair_type
RLEVector::Iterator::selectRun(usint index, usint max_length)
{
  usint value = this->select(index);

  usint len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

pair_type
RLEVector::Iterator::selectNextRun(usint max_length)
{
  usint value = this->selectNext();

  usint len = std::min(max_length, this->run);
  this->run -= len; this->cur += len; this->val += len;

  return pair_type(value, len);
}

bool
RLEVector::Iterator::isSet(usint value)
{
  const RLEVector& par = (const RLEVector&)(this->parent);

  if(value >= par.size) { return false; }

  this->valueLoop(value);

  return (this->val == value);
}

usint
RLEVector::Iterator::countRuns()
{
  const RLEVector& par = (const RLEVector&)(this->parent);

  if(par.items == 0) { return 0; }

  usint runs = 1;
  pair_type res = this->selectRun(0, par.items);
  usint last = res.first + res.second;

  while(last < par.size)
  {
    res = this->selectNextRun(par.items);
    if(res.first < par.size && res.first > last + 1) { runs++; }
    last = res.first + res.second;
  }

  return runs;
}
// FC> ---------------------------- END NOT STUDIED YET ----------------------------------






















//--------------------------------------------------------------------------

RLEEncoder::RLEEncoder(usint block_bytes, usint superblock_size) :
  VectorEncoder(block_bytes, superblock_size),
  run(EMPTY_PAIR)
{
}

RLEEncoder::~RLEEncoder()
{
}

void
RLEEncoder::setBit(usint value)
{
  this->setRun(value,1);
}

void
RLEEncoder::setRun(usint start, usint len)
{
  if(this->items == 0)
  {
    this->setFirstBit(start);
    if(len > 1)
    {
      this->RLEncode(1,len-1);
    }
    return;
  }
  if(start < this->size || len == 0) { return; }

  // Write as much into the buffer as possible.
  usint diff = start + 1 - this->size;
  usint free_bits = this->buffer->bitsLeft();
  usint code_bits = this->buffer->deltaCodeLength(diff);  // FC> distance from previous 1-bit
  
  // <FC>
  usint run_bits = this->buffer->deltaCodeLength(project(diff,len,1,(usint)atoi(getenv("THRESHOLD"))));
  if(free_bits >= code_bits+run_bits) // FC> The whole run must fit into the block now
  {
    free_bits -= code_bits+run_bits;
    this->RLEncode(diff,len);
  }
  else
  {
    this->size = start + 1;
    this->items++;
    // Didn't fit into the block. A new sample & block required.
    this->addNewBlock();
    if(len > 1)
    {
      this->RLEncode(1,len-1);
    }
  }
  // </FC>
}

void
RLEEncoder::addBit(usint value)
{
  this->addRun(value,1);
}

void
RLEEncoder::addRun(usint start, usint len)
{
  if(this->run.second == 0)
  {
    this->run = pair_type(start, len);
  }
  else if(start == this->run.first + this->run.second)
  {
    this->run.second += len;
  }
  else
  {
    this->setRun(this->run.first,this->run.second);
    this->run = pair_type(start, len);
  }
}

void
RLEEncoder::flush()
{
  this->setRun(this->run.first,this->run.second);
  this->run.second = 0;
}


} // namespace CSA
