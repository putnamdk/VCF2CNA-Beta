//------------------------------------------------------------------------------------
//
// genutil.h - module containing useful definitions for genomic processing
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2017 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#ifndef GENUTIL_H
#define GENUTIL_H

#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

const int MAX_EQUIV_INDEL_DISTANCE = 1000; // max base pairs between equivalent indels
const int MAX_POSITION        = 300000000; // max position within a chromosome
const int DEFAULT_BUFFER_SIZE =   1048576; // default buffer size for binary I/O

const int NUM_CHROMOSOMES = 24; // 1 to 22, 23=X, 24=Y
extern const std::string chrLongName [NUM_CHROMOSOMES + 1]; // "chr1" to "chrY"
extern const std::string chrShortName[NUM_CHROMOSOMES + 1]; // "1" to "Y"

uint8_t getChrNumber(const std::string& chrName);

int    stringToInt(const std::string& s);
double stringToDbl(const std::string& s);

inline int  roundit(double d)      { return static_cast<int>(d + 0.5); }
inline bool validPosition(int pos) { return (pos >= 1 && pos <= MAX_POSITION); }

bool isACGT (char ch);
bool isACGTN(char ch);

bool isAllACGT (const std::string& sequence);
bool isAllACGTN(const std::string& sequence);

std::string toupperSequence(const std::string& sequence);
std::string reverseSequence(const std::string& sequence);
std::string invertSequence (const std::string& sequence);

typedef std::vector<std::string> StringVector;

void getDelimitedStrings(const std::string& s, char delimiter, StringVector& v);

//------------------------------------------------------------------------------------

class Variant // represents an indel or SNV
{
public:
   Variant(uint8_t inChrNumber, uint32_t inPosition, const std::string& inSequence);

   Variant(const std::string& s);

   Variant(const Variant& v) // copy constructor
      : chrNumber(v.chrNumber), position(v.position), sequence(v.sequence) { }

   virtual ~Variant() { }

   virtual std::string toString() const;

   virtual bool isInsertion()     const { return (sequence[0] == 'I'); }
   virtual bool isDeletion()      const { return (sequence[0] == 'D'); }
   virtual bool isSubstitution()  const { return (sequence[0] == 'S'); }
   virtual bool isIndel()         const { return (isInsertion() || isDeletion()); }

   uint8_t     chrNumber; // 1 to 22, 23=X, 24=Y
   uint32_t    position;  // 1 to MAX_POSITION
   std::string sequence;  // "Ialt", "Dref", or "Srefalt"
};

typedef std::vector<Variant *> VariantVector;
typedef std::map<std::string, Variant *> VariantMap; // key is sequence

//------------------------------------------------------------------------------------

class Position // represents a position within a chromosome
{
public:
   Position(uint8_t inChrNumber, uint32_t inPosition);
   Position(const std::string& s);
   virtual ~Position();

   virtual std::string toString() const;

   uint8_t    chrNumber; // 1 to 22, 23=X, 24=Y
   uint32_t   position;  // 1 to MAX_POSITION
   VariantMap varmap;    // map of variants at this position
};

typedef std::vector<Position *> PositionVector;
typedef std::map<uint32_t, Position *> PositionMap; // key is position

Variant *saveVariantInPositionMap(PositionMap& pmap, Variant *v);

//------------------------------------------------------------------------------------

class Chromosome
{
public:
   Chromosome(uint8_t inChrNumber);
   Chromosome(const std::string& s);
   virtual ~Chromosome();

   virtual std::string toString() const { return chrLongName[chrNumber]; }

   uint8_t     chrNumber; // 1 to 22, 23=X, 24=Y
   PositionMap posmap;    // map of positions within this chromosome
};

//------------------------------------------------------------------------------------

class BinaryWriter // for writing a binary file
{
public:
   BinaryWriter(size_t bufferSize=DEFAULT_BUFFER_SIZE);
   virtual ~BinaryWriter() { delete[] buf; }

   virtual bool openFile(const char *filename, bool newFile);
   virtual void write_buffer(const void *buffer, size_t numBytes);
   virtual void write_string(const char *string);
   virtual void write_uint8 (uint8_t  value);
   virtual void write_uint16(uint16_t value);
   virtual void write_uint32(uint32_t value);
   virtual void write_uint64(uint64_t value);
   virtual void write_double(double value);
   virtual void flushBuffer();
   virtual void closeFile();

   virtual uint64_t bytesWritten() const { return bytesFlushed + offset; }

   int       fd;
   uint8_t  *buf;
   size_t    bufsize, offset;
   uint64_t  bytesFlushed;
};

//------------------------------------------------------------------------------------

class BinaryReader // for reading a binary file
{
public:
   BinaryReader(size_t bufferSize=DEFAULT_BUFFER_SIZE);
   virtual ~BinaryReader() { delete[] buf; }

   virtual bool openFile(const char *filename);
   virtual void seek(uint64_t byteOffset);
   virtual bool fillBuffer();
   virtual bool read_buffer(uint8_t *buffer, size_t numBytes);
   virtual bool read_string(char *string, size_t maxlen);
   virtual bool read_uint8 (uint8_t&  value);
   virtual bool read_uint16(uint16_t& value);
   virtual bool read_uint32(uint32_t& value);
   virtual bool read_uint64(uint64_t& value);
   virtual bool read_double(double& value);
   virtual bool skipBytes(size_t numBytes);
   virtual void closeFile();

   int      fd;
   uint8_t *buf;
   size_t   bufsize, buflen, offset;
};

//------------------------------------------------------------------------------------

class ReferenceGenome // for representing a reference genome and determining indel
                      // equivalence
{
public:
   ReferenceGenome(const std::string& twobit_filename, uint8_t chrNumber,
		   uint32_t beginpos, uint32_t endpos, std::string chrName="");
   virtual ~ReferenceGenome() { delete[] sequence; }

   virtual char getBase(uint32_t pos) const;

   virtual bool validDeletion(uint32_t pos, const std::string& seq) const;

   virtual bool equivalentInsertions(uint32_t pos1, const std::string& seq1,
		                     uint32_t pos2, const std::string& seq2) const;

   virtual bool equivalentDeletions (uint32_t pos1, const std::string& seq1,
		                     uint32_t pos2, const std::string& seq2) const;

   uint32_t  begin, end;
   char     *sequence;
};

StringVector *getChromosomeNames(const std::string& twobit_filename);

//------------------------------------------------------------------------------------

class BambinoParser // for parsing lines in a Bambino file
{
public:
   BambinoParser(const std::string& headingLine);
   virtual ~BambinoParser() { }

   virtual bool parseLine(const std::string& line, std::string& chrName,
		          int& position, std::string& variantType,
		          std::string& ref, std::string& alt,
		          int& refCount, int& altCount) const;

   // column numbers of columns of interest
   int chrCol, posCol, typeCol, refCol, altCol, refCountCol, altCountCol;

   int numColumns;
};

//------------------------------------------------------------------------------------

class BambinoParserTumor : public BambinoParser // includes tumor fields
{
public:
   BambinoParserTumor(const std::string& headingLine);
   virtual ~BambinoParserTumor() { }

   virtual bool parseLine(const std::string& line, std::string& chrName,
		          int& position, std::string& variantType,
		          std::string& ref, std::string& alt,
		          int& refCount, int& altCount,
			  int& refTumorCount, int& altTumorCount,
			  std::string& tumorSample) const;

   // column numbers of columns of interest
   int refTumorCountCol, altTumorCountCol, tumorSampleCol;
};

//------------------------------------------------------------------------------------

class SequenceTrie // represents a set of sequences as a trie
{
public:
   SequenceTrie();
   virtual ~SequenceTrie();

   virtual void addSequence (const std::string& sequence);
   virtual bool findSequence(const std::string& sequence) const;

   bool endOfSequence;
   SequenceTrie *link[4]; // for A, C, G, T
};

//------------------------------------------------------------------------------------

class NumberSet // aggregates a set of numbers
{
public:
   NumberSet();
   virtual ~NumberSet() { }

   void addNumber(double x);

   double average()  const;
   double variance() const;
   double stdev()    const;

   int n; // number of numbers in set
   double min, max, sum, sumsq;
};

//------------------------------------------------------------------------------------

class ObservationSet // used to compute Pearson correlation
{
public:
   ObservationSet();
   virtual ~ObservationSet() { }

   void addObservation(double x, double y);

   double pearsonCorrelationCoefficient() const;

   int n; // number of observations in set
   double sumx, sumy, sumxx, sumyy, sumxy;
};

//------------------------------------------------------------------------------------

class SpearmanObservation
{
public:
   SpearmanObservation(double inX, double inY)
      : x(inX), y(inY), xrank(-1), yrank(-1) { }

   virtual ~SpearmanObservation() { }

   double x, y, xrank, yrank;
};

typedef std::vector<SpearmanObservation> SpearmanObservationVector;

class SpearmanObservationSet // used to compute Spearman rank correlation
{
public:
   SpearmanObservationSet()
      : obs() { }

   virtual ~SpearmanObservationSet() { }

   void addObservation(double x, double y)
   { obs.push_back(SpearmanObservation(x, y)); }

   double rankCorrelationCoefficient();

   SpearmanObservationVector obs;
};

//------------------------------------------------------------------------------------
#endif
