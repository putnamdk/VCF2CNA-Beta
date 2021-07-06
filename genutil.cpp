//------------------------------------------------------------------------------------
//
// genutil.cpp - module containing useful definitions for genomic processing
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2017 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "genutil.h"

const std::string chrLongName[NUM_CHROMOSOMES + 1] =
{
   "",      "chr1",  "chr2",  "chr3",  "chr4",
   "chr5",  "chr6",  "chr7",  "chr8",  "chr9",
   "chr10", "chr11", "chr12", "chr13", "chr14",
   "chr15", "chr16", "chr17", "chr18", "chr19",
   "chr20", "chr21", "chr22", "chrX",  "chrY"
};

const std::string chrShortName[NUM_CHROMOSOMES + 1] =
{
   "",   "1",  "2",  "3",  "4",
   "5",  "6",  "7",  "8",  "9",
   "10", "11", "12", "13", "14",
   "15", "16", "17", "18", "19",
   "20", "21", "22", "X",  "Y"
};

//------------------------------------------------------------------------------------
// getChrNumber() returns the chromosome number (1-24) for a given chromosome name;
// zero is returned if the chromosome name is unrecognized

uint8_t getChrNumber(const std::string& chrName)
{
   if (chrName.length() > 3)
   {
      for (uint8_t i = 1; i <= NUM_CHROMOSOMES; i++)
         if (chrName == chrLongName[i])
            return i;
   }
   else
   {
      for (uint8_t i = 1; i <= NUM_CHROMOSOMES; i++)
         if (chrName == chrShortName[i])
            return i;
   }

   return 0; // unrecognized chromosome name
}

//------------------------------------------------------------------------------------
// stringToInt() converts a string to an integer; -1 is returned if the conversion
// cannot be performed

int stringToInt(const std::string& s)
{
   int n = s.length();

   if (n < 1 || n > 10)
      return -1;

   int value = 0;

   for (int i = 0; i < n; i++)
   {
      char c = s[i];

      if (std::isdigit(c))
         value = 10 * value + c - '0';
      else
         return -1;
   }

   return value;
}

//------------------------------------------------------------------------------------
// stringToDbl() converts a string to a double; -1.0 is returned if the conversion
// cannot be performed

double stringToDbl(const std::string& s)
{
   std::stringstream stream(s);
   double value = -1.0;

   stream >> value;
   return value;
}

//------------------------------------------------------------------------------------
// isACGT() returns true if the given character is A, C, G or T

bool isACGT(char ch)
{
   ch = std::toupper(ch);

   return (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T');
}

//------------------------------------------------------------------------------------
// isACGTN() returns true if the given character is A, C, G, T or N

bool isACGTN(char ch)
{
   ch = std::toupper(ch);

   return (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'N');
}

//------------------------------------------------------------------------------------
// isAllACGT() returns true if all characters in the sequence are A, C, G or T

bool isAllACGT(const std::string& sequence)
{
   int len = sequence.length();

   for (int i = 0; i < len; i++)
      if (!isACGT(sequence[i]))
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// isAllACGTN() returns true if all characters in the sequence are A, C, G, T or N

bool isAllACGTN(const std::string& sequence)
{
   int len = sequence.length();

   for (int i = 0; i < len; i++)
      if (!isACGTN(sequence[i]))
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// toupperSequence() converts all letters in a sequence to uppercase

std::string toupperSequence(const std::string& sequence)
{
   std::string upperSeq = "";

   int len = sequence.length();

   for (int i = 0; i < len; i++)
   {
      char ch   = std::toupper(sequence[i]);
      upperSeq += ch;
   }

   return upperSeq;
}

//------------------------------------------------------------------------------------
// reverseSequence() reverses the order of the characters in a sequence

std::string reverseSequence(const std::string& sequence)
{
   std::string reverseSeq = "";

   int len = sequence.length();

   for (int i = len - 1; i >= 0; i--)
      reverseSeq += sequence[i];

   return reverseSeq;
}

//------------------------------------------------------------------------------------
// invertSequence() swaps A and T, and C and G in the sequence

std::string invertSequence(const std::string& sequence)
{
   std::string inverseSeq = "";

   int len = sequence.length();

   for (int i = 0; i < len; i++)
   {
      char ch;

      switch (sequence[i])
      {
         case 'A': ch = 'T'; break;
	 case 'a': ch = 't'; break;

	 case 'T': ch = 'A'; break;
	 case 't': ch = 'a'; break;

         case 'C': ch = 'G'; break;
	 case 'c': ch = 'g'; break;

	 case 'G': ch = 'C'; break;
	 case 'g': ch = 'c'; break;

	 default : ch = sequence[i];
      }

      inverseSeq += ch;
   }

   return inverseSeq;
}

//------------------------------------------------------------------------------------
// getDelimitedStrings() extracts delimited string values from a string and appends
// them to a string vector

void getDelimitedStrings(const std::string& s, char delimiter, StringVector& v)
{
   int i = 0, len = s.length();

   do
   {
      std::string value = "";

      while (i < len && s[i] != delimiter)
	 value += s[i++];

      v.push_back(value);
   }
   while (i++ < len);
}

//------------------------------------------------------------------------------------
// Variant::Variant(uint8_t, uint32_t, const std::string&) validates the arguments
// before constructing a Variant object

Variant::Variant(uint8_t inChrNumber, uint32_t inPosition,
	         const std::string& inSequence)
   : chrNumber(inChrNumber), position(inPosition), sequence(inSequence)
{
   if (chrNumber >= 1 && chrNumber <= NUM_CHROMOSOMES && validPosition(position) &&
       sequence.length() >= 2)
   {
      sequence = toupperSequence(sequence);
      std::string refalt = sequence.substr(1);

      if (sequence[0] == 'I' && isAllACGT(refalt))
         return; // this is a valid insertion

      if (sequence[0] == 'D' && isAllACGTN(refalt))
         return; // this is a valid deletion

      if (sequence[0] == 'S' && isAllACGT(refalt) && refalt.length() == 2 &&
	  refalt[0] != refalt[1])
         return; // this is a valid SNV
   }

   throw std::runtime_error("invalid variant specification");
}

//------------------------------------------------------------------------------------
// Variant::Variant(const std::string&) parses and validates a variant string before
// constructing a Variant object

Variant::Variant(const std::string& s)
{
   int i = 0, len = s.length();

   while (i < len && s[i] != ':' && s[i] != '.') // look for colon or period
      i++;

   if (i > 0 && i < len - 5) // found separator between chromosome and position
   {
      std::string chrName = s.substr(0, i);
      chrNumber = getChrNumber(chrName);

      if (chrNumber > 0) // chromosome name is valid
      {
         int j = i + 1;

         while (j < len && s[j] != '.') // look for period
            j++;

	 if (j > i + 1 && j < len - 3) // found period between position and ref
	 {
            std::string posString = s.substr(i + 1, j - i - 1);
	    int pos = stringToInt(posString);

	    if (validPosition(pos)) // position is valid
	    {
	       position = pos;

               int k = j + 1;

	       while (k < len && s[k] != '.') // look for period
	          k++;

	       if (k > j + 1 && k < len - 1) // found period between ref and alt
	       {
	          std::string ref = s.substr(j + 1, k - j - 1);
		  ref = toupperSequence(ref);

		  std::string alt = s.substr(k + 1);
		  alt = toupperSequence(alt);

                  if (ref == "-" && isAllACGT(alt))
	          {
		     // found a valid insertion
		     sequence = "I" + alt;
		     return;
		  }

		  if (alt == "-" && isAllACGTN(ref))
		  {
	             // found a valid deletion
		     sequence = "D" + ref;
		     return;
		  }

		  if (ref.length() == 1 && isACGT(ref[0]) &&
		      alt.length() == 1 && isACGT(alt[0]) && ref[0] != alt[0])
		  {
	             // found a valid SNV
		     sequence = "S" + ref + alt;
		     return;
		  }
	       }
	    }
	 }
      }
   }

   throw std::runtime_error("invalid variant specification \"" + s + "\"");
}

//------------------------------------------------------------------------------------
// Variant::toString() returns the string representation of a variant

std::string Variant::toString() const
{
   std::stringstream stream;

   stream << chrLongName[chrNumber] << "." << position << ".";

   if (sequence[0] == 'I')
      stream << "-." << sequence.substr(1);
   else if (sequence[0] == 'D')
      stream << sequence.substr(1) << ".-";
   else if (sequence[0] == 'S')
      stream << sequence[1] << "." << sequence[2];

   return stream.str();
}

//------------------------------------------------------------------------------------
// Position::Position(uint8_t, uint32_t) validates the arguments before constructing a
// Position object

Position::Position(uint8_t inChrNumber, uint32_t inPosition)
   : chrNumber(inChrNumber), position(inPosition), varmap()
{
   if (chrNumber >= 1 && chrNumber <= NUM_CHROMOSOMES && validPosition(position))
      return; // this is a valid position

   throw std::runtime_error("invalid position specification");
}

//------------------------------------------------------------------------------------
// Position::Position(const std::string&) parses and validates a position string
// before constructing a Position object

Position::Position(const std::string& s)
   : varmap()
{
   int i = 0, len = s.length();

   while (i < len && s[i] != ':' && s[i] != '.') // look for colon or period
      i++;

   if (i > 0 && i < len - 1) // found separator between chromosome and position
   {
      std::string chrName = s.substr(0, i);
      chrNumber = getChrNumber(chrName);

      if (chrNumber > 0) // chromosome name is valid
      {
         std::string posString = s.substr(i + 1);
	 int pos = stringToInt(posString);

	 if (validPosition(pos)) // position is valid
	 {
            position = pos;
            return;
	 }
      }
   }

   throw std::runtime_error("invalid position specification \"" + s + "\"");
}

//------------------------------------------------------------------------------------
// Position::~Position() de-allocates all variants stored in the variant map

Position::~Position()
{
   for (VariantMap::iterator vpos = varmap.begin(); vpos != varmap.end(); ++vpos)
   {
      Variant *variant = vpos->second;
      delete variant;
      vpos->second = NULL;
   }
}

//------------------------------------------------------------------------------------
// Position::toString() returns the string representation of a position

std::string Position::toString() const
{
   std::stringstream stream;

   stream << chrLongName[chrNumber] << "." << position;

   return stream.str();
}

//------------------------------------------------------------------------------------
// saveVariantInPositionMap() saves the given variant in the specified position map;
// if the variant is already in the map, the Variant object passed to this function
// is a duplicate and is deleted; the given variant is returned by the function,
// except that when it is a duplicate, the original variant is returned instead;
// in all cases, the returned value is the variant stored in the map

Variant *saveVariantInPositionMap(PositionMap& pmap, Variant *v)
{
   PositionMap::iterator ppos = pmap.find(v->position);

   if (ppos == pmap.end())
   {
      // this position is not yet in the map
      Position *p = new Position(v->chrNumber, v->position);
      pmap.insert(std::make_pair(v->position, p));
      p->varmap.insert(std::make_pair(v->sequence, v));
      return v;
   }

   Position *p = ppos->second;

   VariantMap::iterator vpos = p->varmap.find(v->sequence);

   if (vpos == p->varmap.end())
   {
      // this variant is not yet in the map
      p->varmap.insert(std::make_pair(v->sequence, v));
      return v;
   }

   // this variant is already in the map
   delete v;            // discard duplicate variant
   return vpos->second; // return original variant
}

//------------------------------------------------------------------------------------
// Chromosome::Chromosome(uint8_t) validates the argument before constructing a
// Chromosome object

Chromosome::Chromosome(uint8_t inChrNumber)
   : chrNumber(inChrNumber), posmap()
{
   if (chrNumber >= 1 && chrNumber <= NUM_CHROMOSOMES)
      return; // this is a valid chromosome

   throw std::runtime_error("invalid chromosome specification");
}

//------------------------------------------------------------------------------------
// Chromosome::Chromosome(const std::string&) validates a chromosome name before
// constructing a Chromosome object

Chromosome::Chromosome(const std::string& s)
   : posmap()
{
   chrNumber = getChrNumber(s);

   if (chrNumber > 0)
      return; // this is a valid chromosome

   throw std::runtime_error("invalid chromosome specification \"" + s + "\"");
}

//------------------------------------------------------------------------------------
// Chromosome::~Chromosome() de-allocates all positions stored in the position map

Chromosome::~Chromosome()
{
   for (PositionMap::iterator ppos = posmap.begin(); ppos != posmap.end(); ++ppos)
   {
      Position *position = ppos->second;
      delete position;
      ppos->second = NULL;
   }
}

//------------------------------------------------------------------------------------
// BinaryWriter::BinaryWriter() allocates an internal buffer

BinaryWriter::BinaryWriter(size_t bufferSize)
   : fd(-1), bufsize(bufferSize), offset(0), bytesFlushed(0)
{
   buf = new uint8_t[bufsize];
}

//------------------------------------------------------------------------------------
// BinaryWriter::openFile() creates a new file for writing if newFile is true, and
// opens an existing file for writing if newFile is false; true is returned if
// successful

bool BinaryWriter::openFile(const char *filename, bool newFile)
{
   if (fd != -1) // file is already open
      return false;

   if (newFile)
      fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC,
	        S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
   else
      fd = open(filename, O_WRONLY);

   if (fd == -1) // error
      return false;

   offset       = 0;
   bytesFlushed = 0;
   return true;
}

//------------------------------------------------------------------------------------
// BinaryWriter::write_buffer() writes a buffer of bytes

void BinaryWriter::write_buffer(const void *buffer, size_t numBytes)
{
   if (numBytes > bufsize)
      throw std::runtime_error("binary write buffer is too small");

   if (offset + numBytes > bufsize)
      flushBuffer();

   std::memcpy(&buf[offset], buffer, numBytes);

   offset += numBytes;
}

//------------------------------------------------------------------------------------
// BinaryWriter::write_string() writes a C-style string, including the trailing null
// byte

void BinaryWriter::write_string(const char *string)
{
   write_buffer(string, std::strlen(string) + 1);
}

//------------------------------------------------------------------------------------
// BinaryWriter::write_uint8() writes a one-byte integer

void BinaryWriter::write_uint8(uint8_t value)
{
   write_buffer(&value, 1);
}

//------------------------------------------------------------------------------------
// BinaryWriter::write_uint16() writes a two-byte integer

void BinaryWriter::write_uint16(uint16_t value)
{
   uint8_t byte[2];
   byte[0] = static_cast<uint8_t>(value >> 8);
   byte[1] = static_cast<uint8_t>(value & 0xFF);

   write_buffer(byte, 2);
}

//------------------------------------------------------------------------------------
// BinaryWriter::write_uint32() writes a four-byte integer

void BinaryWriter::write_uint32(uint32_t value)
{
   uint8_t byte[4];

   for (int i = 0; i < 4; i++)
      byte[i] = static_cast<uint8_t>((value >> (24 - 8 * i)) & 0xFF);

   write_buffer(byte, 4);
}

//------------------------------------------------------------------------------------
// BinaryWriter::write_uint64() writes an eight-byte integer

void BinaryWriter::write_uint64(uint64_t value)
{
   uint8_t byte[8];

   for (int i = 0; i < 8; i++)
      byte[i] = static_cast<uint8_t>((value >> (56 - 8 * i)) & 0xFF);

   write_buffer(byte, 8);
}

//------------------------------------------------------------------------------------
// BinaryWriter::write_double() writes a double-precision floating-point value

void BinaryWriter::write_double(double value)
{
   write_uint64(*(uint64_t *)(&value));
}

//------------------------------------------------------------------------------------
// BinaryWriter::flushBuffer() writes the internal buffer to the file

void BinaryWriter::flushBuffer()
{
   if (fd == -1)
      throw std::runtime_error("binary file not open");

   if (offset == 0)
      return; // nothing to flush

   ssize_t bytes = write(fd, buf, offset);
   if (bytes != offset)
      throw std::runtime_error("binary file write error");

   offset = 0;
   bytesFlushed += bytes;
}

//------------------------------------------------------------------------------------
// BinaryWriter::closeFile() writes the internal buffer and closes the file

void BinaryWriter::closeFile()
{
   if (fd == -1) // no file is open
      return;

   flushBuffer();

   if (close(fd) == -1)
      throw std::runtime_error("binary file close error");

   fd           = -1;
   offset       =  0;
   bytesFlushed =  0;
}

//------------------------------------------------------------------------------------
// BinaryReader::BinaryReader() allocates an internal buffer

BinaryReader::BinaryReader(size_t bufferSize)
   : fd(-1), bufsize(bufferSize), buflen(0), offset(0)
{
   buf = new uint8_t[bufsize];
}

//------------------------------------------------------------------------------------
// BinaryReader::openFile() opens an existing file for reading; true is returned if
// successful

bool BinaryReader::openFile(const char *filename)
{
   if (fd != -1) // file is already open
      return false;

   fd = open(filename, O_RDONLY);
   if (fd == -1) // error
      return false;

   buflen = 0;
   offset = 0;
   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::seek() performs a seek operation to the specified byte offset from
// the beginning of the file

void BinaryReader::seek(uint64_t byteOffset)
{
   if (fd == -1)
      throw std::runtime_error("binary file not open");

   if (lseek(fd, byteOffset, SEEK_SET) == -1)
      throw std::runtime_error("binary file seek error");

   buflen = 0;
   offset = 0;
}

//------------------------------------------------------------------------------------
// BinaryReader::fillBuffer() reads from the file and puts the bytes read into the
// internal buffer; false is returned when EOF has been reached

bool BinaryReader::fillBuffer()
{
   if (fd == -1)
      throw std::runtime_error("binary file not open");

   ssize_t bytes = read(fd, buf, bufsize);
   if (bytes == -1)
      throw std::runtime_error("binary file read error");

   if (bytes == 0) // reached EOF
      return false;

   buflen = bytes;
   offset = 0;
   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::read_buffer() reads a buffer of bytes; false is returned if EOF is
// reached before the requested number of bytes has been read

bool BinaryReader::read_buffer(uint8_t *buffer, size_t numBytes)
{
   for (int i = 0; i < numBytes; i++)
   {
      if (offset >= buflen && !fillBuffer())
         return false;

      buffer[i] = buf[offset++];
   }

   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::read_string() reads a C-style string, copying the characters of the
// string up to a maximum length or up to and including the trailing null byte; false
// is returned if EOF is encountered

bool BinaryReader::read_string(char *string, size_t maxlen)
{
   for (int i = 0; i < maxlen; i++)
   {
      uint8_t byte;
      if (!read_buffer(&byte, 1))
         return false;

      string[i] = static_cast<char>(byte);
      if (byte == 0)
         break;
   }

   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::read_uint8() reads a one-byte integer; false is returned if EOF is
// encountered

bool BinaryReader::read_uint8(uint8_t& value)
{
   return read_buffer(&value, 1);
}

//------------------------------------------------------------------------------------
// BinaryReader::read_uint16() reads a two-byte integer; false is returned if EOF is
// encountered

bool BinaryReader::read_uint16(uint16_t& value)
{
   uint8_t byte[2];
   if (!read_buffer(byte, 2))
      return false;

   value = (byte[0] << 8) + byte[1];
   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::read_uint32() reads a four-byte integer; false is returned if EOF is
// encountered

bool BinaryReader::read_uint32(uint32_t& value)
{
   uint8_t byte[4];
   if (!read_buffer(byte, 4))
      return false;

   value = 0;
   for (int i = 0; i < 4; i++)
      value = (value << 8) + byte[i];

   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::read_uint64() reads an eight-byte integer; false is returned if EOF
// is encountered

bool BinaryReader::read_uint64(uint64_t& value)
{
   uint8_t byte[8];
   if (!read_buffer(byte, 8))
      return false;

   value = 0;
   for (int i = 0; i < 8; i++)
      value = (value << 8) + byte[i];

   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::read_double() reads a double-precision floating-point value; false is
// returned if EOF is encountered

bool BinaryReader::read_double(double& value)
{
   uint64_t ivalue;
   if (!read_uint64(ivalue))
      return false;

   value = *(double *)(&ivalue);

   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::skipBytes() skips the specified number of bytes in the input file;
// false is returned if EOF is encountered

bool BinaryReader::skipBytes(size_t numBytes)
{
   size_t bytesRemaining = (buflen > offset ? buflen - offset : 0);

   if (numBytes <= bytesRemaining) // all bytes to skip are already in the buffer
      offset += numBytes;
   else
   {
      numBytes -= bytesRemaining;

      while (true) // refill the buffer until enough bytes have been skipped
      {
         if (!fillBuffer())
            return false;

	 if (numBytes <= buflen)
	 {
            offset = numBytes;
            break;
	 }

	 numBytes -= buflen;
      }
   }

   return true;
}

//------------------------------------------------------------------------------------
// BinaryReader::closeFile() closes the file

void BinaryReader::closeFile()
{
   if (fd == -1) // no file is open
      return;

   if (close(fd) == -1)
      throw std::runtime_error("binary file close error");

   fd     = -1;
   buflen =  0;
   offset =  0;
}

//------------------------------------------------------------------------------------
// swap_uint32() swaps the byte ordering of a four-byte unsigned integer

static uint32_t swap_uint32(uint32_t value)
{
   return ((value << 24) + ((value & 0xFF00) << 8) + ((value & 0xFF0000) >> 8) +
	   (value >> 24));
}

//------------------------------------------------------------------------------------
// read_uint32() reads a four-byte unsigned integer from a 2bit file

static uint32_t read_uint32(BinaryReader& reader, bool swapBytes)
{
   uint32_t value;

   if (!reader.read_uint32(value))
      throw std::runtime_error("truncated 2bit file");

   return (swapBytes ? swap_uint32(value) : value);
}

//------------------------------------------------------------------------------------
// ReferenceGenome::ReferenceGenome() extracts a DNA sequence from a 2bit file and
// stores it in a ReferenceGenome object

ReferenceGenome::ReferenceGenome(const std::string& twobit_filename,
	                         uint8_t chrNumber, uint32_t beginpos,
				 uint32_t endpos, std::string chrName)
{
   // if no chromosome name specified, validate the specified chromosome number
   if (chrName == "" && (chrNumber == 0 || chrNumber > NUM_CHROMOSOMES))
      throw std::runtime_error("invalid chromosome specification");

   BinaryReader reader;

   if (!reader.openFile(twobit_filename.c_str()))
      throw std::runtime_error("unable to open " + twobit_filename);

   const uint32_t EXPECTED_SIGNATURE = 0x1A412743;
   uint32_t signature;

   if (!reader.read_uint32(signature) || signature != EXPECTED_SIGNATURE &&
       signature != swap_uint32(EXPECTED_SIGNATURE))
      throw std::runtime_error(twobit_filename + " is not a 2bit file");

   bool swapBytes = (signature != EXPECTED_SIGNATURE);

   uint32_t version  = read_uint32(reader, swapBytes);
   uint32_t chrCount = read_uint32(reader, swapBytes);
   uint32_t reserved = read_uint32(reader, swapBytes);

   // get the byte offset of our chromosome

   uint32_t chrOffset = 0;

   for (uint32_t i = 0; i < chrCount && chrOffset == 0; i++)
   {
      uint8_t nameLength, nameBuffer[255];
      if (!reader.read_uint8(nameLength) || 
	  !reader.read_buffer(nameBuffer, nameLength))
         throw std::runtime_error("truncated 2bit file " + twobit_filename);

      uint32_t offset = read_uint32(reader, swapBytes);

      std::string name = "";
      for (uint32_t j = 0; j < nameLength; j++)
         name += static_cast<char>(nameBuffer[j]);

      if (chrName == "" &&
          (name == chrShortName[chrNumber] || name == chrLongName[chrNumber]) ||
	  chrName != "" && name == chrName)
         chrOffset = offset;
   }

   if (chrOffset == 0)
   {
      if (chrName == "")
         chrName = chrShortName[chrNumber];

      throw std::runtime_error("chromosome " + chrName + " not found in " +
		               twobit_filename);
   }

   // jump to the chromosome header

   reader.seek(chrOffset);
   
   // read the chromosome header
   
   uint32_t numBases = read_uint32(reader, swapBytes);

   begin = beginpos;
   end   = (endpos > numBases ? numBases : endpos);

   if (begin == 0 || begin > end)
      throw std::runtime_error("invalid begin position");

   sequence = new char[end - begin + 1]; // allocate the internal buffer

   std::vector<uint32_t> nstart, nstop;

   uint32_t nBlockCount = read_uint32(reader, swapBytes);

   for (uint32_t i = 0; i < nBlockCount; i++)
      nstart.push_back(read_uint32(reader, swapBytes) + 1); // convert to 1-based

   for (uint32_t i = 0; i < nBlockCount; i++)
      nstop.push_back(nstart[i] + read_uint32(reader, swapBytes) - 1);

   uint32_t maskBlockCount = read_uint32(reader, swapBytes);

   uint32_t dnaOffset = chrOffset +
      sizeof(uint32_t) * (2 * nBlockCount + 2 * maskBlockCount + 4);

   // jump to the first DNA byte in the selected range

   reader.seek(dnaOffset + ((begin - 1) >> 2));

   // read the two-bit sequence data and store it in the internal buffer as characters

   const char symbol[4] = { 'T', 'C', 'A', 'G' };

   uint8_t byte;
   bool readByte = true;

   for (uint32_t pos = begin; pos <= end; pos++)
   {
      if (readByte && !reader.read_uint8(byte))
         throw std::runtime_error("truncated 2bit file " + twobit_filename);

      uint32_t shift = 2 * (3 - ((pos - 1) & 3));
      uint32_t index = (byte >> shift) & 3;

      sequence[pos - begin] = symbol[index];

      readByte = (shift == 0);
   }

   // now mark the unknown regions in the sequence

   for (uint32_t i = 0; i < nBlockCount; i++)
      if (nstart[i] <= end && nstop[i] >= begin)
      {
         uint32_t start = (nstart[i] < begin ? begin : nstart[i]);
	 uint32_t stop  = (nstop[i]  > end   ? end   : nstop[i]);

	 for (uint32_t pos = start; pos <= stop; pos++)
            sequence[pos - begin] = 'N'; // unknown position
      }

   reader.closeFile();
}

//------------------------------------------------------------------------------------
// ReferenceGenome::getBase() returns the nucleotide (A, C, G, T or N) at the
// specified position

char ReferenceGenome::getBase(uint32_t pos) const
{
   if (pos < begin || pos > end)
      return 'N';
   else
      return sequence[pos - begin];
}

//------------------------------------------------------------------------------------
// ReferenceGenome::validDeletion() returns true if the given deletion has a sequence
// that matches the reference genome

bool ReferenceGenome::validDeletion(uint32_t pos, const std::string& seq) const
{
   uint32_t n = seq.length();

   for (uint32_t i = 0; i < n; i++)
      if (seq[i] != getBase(pos + i))
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// ReferenceGenome::equivalentInsertions() returns true if two insertions are
// equivalent; it uses the algorithms given in S.V. Rice, "Determining Whether Two
// Indels Are Equivalent," St. Jude, 2015

bool ReferenceGenome::equivalentInsertions(uint32_t pos1, const std::string& seq1,
		                           uint32_t pos2, const std::string& seq2)
                                          const
{
   if (seq1.length() != seq2.length())
      return false;

   if (pos1 == pos2)
      return (seq1 == seq2);

   uint32_t j; // position before first  insertion
   uint32_t k; // position before second insertion

   if (pos1 < pos2)
   {
      j = pos1 - 1;
      k = pos2 - 1;
   }
   else // pos1 > pos2
   {
      j = pos2 - 1;
      k = pos1 - 1;
   }

   const std::string& v = (pos1 < pos2 ? seq1 : seq2);
   const std::string& w = (pos1 < pos2 ? seq2 : seq1);

   uint32_t m = k - j;
   uint32_t n = seq1.length();

   uint32_t vindex, windex, sindex;

   if (m < n)
   {
      for (vindex = 0, windex = n - m, sindex = j + 1;
	   vindex < m; vindex++, windex++, sindex++)
         if (v[vindex] != w[windex] || v[vindex] != getBase(sindex))
            return false;

      for (vindex = m, windex = 0; vindex < n; vindex++, windex++)
         if (v[vindex] != w[windex])
            return false;

      return true;
   }

   if (m == n)
   {
      for (vindex = 0, windex = 0, sindex = j + 1;
	   vindex < n; vindex++, windex++, sindex++)
         if (v[vindex] != w[windex] || v[vindex] != getBase(sindex))
            return false;

      return true;
   }

   // m > n

   for (vindex = 0, sindex = j + 1; vindex < n; vindex++, sindex++)
      if (v[vindex] != getBase(sindex))
         return false;

   for (windex = 0, sindex = k - n + 1; windex < n; windex++, sindex++)
      if (w[windex] != getBase(sindex))
         return false;

   uint32_t last = k - n;
   for (sindex = j + 1; sindex <= last; sindex++)
      if (getBase(sindex) != getBase(sindex + n))
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// ReferenceGenome::equivalentDeletions() returns true if two deletions are
// equivalent; it uses the algorithms given in S.V. Rice, "Determining Whether Two
// Indels Are Equivalent," St. Jude, 2015

bool ReferenceGenome::equivalentDeletions(uint32_t pos1, const std::string& seq1,
		                          uint32_t pos2, const std::string& seq2)
                                         const
{
   if (seq1.length() != seq2.length())
      return false;

   if (pos1 == pos2)
      return (seq1 == seq2);

   uint32_t j; // starting position of first  deletion
   uint32_t k; // starting position of second deletion

   if (pos1 < pos2)
   {
      j = pos1;
      k = pos2;
   }
   else // pos1 > pos2
   {
      j = pos2;
      k = pos1;
   }

   uint32_t n = seq1.length();

   for (uint32_t sindex = j; sindex < k; sindex++)
      if (getBase(sindex) != getBase(sindex + n))
         return false;

   return true;
}

//------------------------------------------------------------------------------------
// getChromosomeNames() returns a string vector containing all of the chromosome names
// in a 2bit file; the caller must delete the string vector when finished with it

StringVector *getChromosomeNames(const std::string& twobit_filename)
{
   BinaryReader reader;

   if (!reader.openFile(twobit_filename.c_str()))
      throw std::runtime_error("unable to open " + twobit_filename);

   const uint32_t EXPECTED_SIGNATURE = 0x1A412743;
   uint32_t signature;

   if (!reader.read_uint32(signature) || signature != EXPECTED_SIGNATURE &&
       signature != swap_uint32(EXPECTED_SIGNATURE))
      throw std::runtime_error(twobit_filename + " is not a 2bit file");

   bool swapBytes = (signature != EXPECTED_SIGNATURE);

   uint32_t version  = read_uint32(reader, swapBytes);
   uint32_t chrCount = read_uint32(reader, swapBytes);
   uint32_t reserved = read_uint32(reader, swapBytes);

   StringVector *chrName = new StringVector();

   for (uint32_t i = 0; i < chrCount; i++)
   {
      uint8_t nameLength, nameBuffer[255];
      if (!reader.read_uint8(nameLength) || 
	  !reader.read_buffer(nameBuffer, nameLength))
         throw std::runtime_error("truncated 2bit file " + twobit_filename);

      uint32_t offset = read_uint32(reader, swapBytes);

      std::string name = "";
      for (uint32_t j = 0; j < nameLength; j++)
         name += static_cast<char>(nameBuffer[j]);

      chrName->push_back(name);
   }

   reader.closeFile();

   return chrName;
}

//------------------------------------------------------------------------------------
// BambinoParser::BambinoParser() parses a heading line from a Bambino file and saves
// the number of columns in the file and the column numbers of columns of interest

BambinoParser::BambinoParser(const std::string& headingLine)
   : chrCol(-1), posCol(-1), typeCol(-1), refCol(-1), altCol(-1), refCountCol(-1),
     altCountCol(-1)
{
   StringVector heading;
   getDelimitedStrings(headingLine, '\t', heading);

   numColumns = heading.size();

   for (int i = 0; i < numColumns; i++)
      if (heading[i] == "Chr")                           chrCol      = i;
      else if (heading[i] == "Pos")                      posCol      = i;
      else if (heading[i] == "Type")                     typeCol     = i;
      else if (heading[i] == "Chr_Allele")               refCol      = i;
      else if (heading[i] == "Alternative_Allele")       altCol      = i;
      else if (heading[i] == "reference_normal_count")   refCountCol = i;
      else if (heading[i] == "alternative_normal_count") altCountCol = i;

   if (chrCol < 0 || posCol < 0 || typeCol < 0 || refCol < 0 || altCol < 0 ||
       refCountCol < 0 || altCountCol < 0)
      throw std::runtime_error("missing column(s) in Bambino file");
}

//------------------------------------------------------------------------------------
// BambinoParser::parseLine() parses a variant line read from a Bambino file; if it is
// valid, selected column values are passed back to the caller in reference parameters
// and true is returned by the function; otherwise, false is returned

bool BambinoParser::parseLine(const std::string& line, std::string& chrName,
		              int& position, std::string& variantType,
		              std::string& ref, std::string& alt,
		              int& refCount, int& altCount) const
{
   StringVector value;
   getDelimitedStrings(line, '\t', value);

   if (value.size() != numColumns)
      return false; // unexpected number of columns in line

   position = stringToInt(value[posCol]);
   refCount = stringToInt(value[refCountCol]);
   altCount = stringToInt(value[altCountCol]);

   if (position < 0 || refCount < 0 || altCount < 0)
      return false; // unable to convert strings to integers

   chrName     = value[chrCol];
   variantType = value[typeCol];
   ref         = value[refCol];
   alt         = value[altCol];

   return true;
}

//------------------------------------------------------------------------------------
// BambinoParserTumor::BambinoParserTumor() parses a heading line from a Bambino file
// and saves the number of columns in the file and the column numbers of columns of
// interest

BambinoParserTumor::BambinoParserTumor(const std::string& headingLine)
   : BambinoParser(headingLine), refTumorCountCol(-1), altTumorCountCol(-1),
     tumorSampleCol(-1)
{
   StringVector heading;
   getDelimitedStrings(headingLine, '\t', heading);

   for (int i = 0; i < numColumns; i++)
      if (heading[i] == "TumorSample")                  tumorSampleCol   = i;
      else if (heading[i] == "reference_tumor_count")   refTumorCountCol = i;
      else if (heading[i] == "alternative_tumor_count") altTumorCountCol = i;

   if (refTumorCountCol < 0 || altTumorCountCol < 0)
      throw std::runtime_error("missing column(s) in Bambino file");
}

//------------------------------------------------------------------------------------
// BambinoParserTumor::parseLine() parses a variant line read from a Bambino file; if
// it is valid, selected column values are passed back to the caller in reference
// parameters and true is returned by the function; otherwise, false is returned

bool BambinoParserTumor::parseLine(const std::string& line, std::string& chrName,
		                   int& position, std::string& variantType,
		                   std::string& ref, std::string& alt,
		                   int& refCount, int& altCount,
				   int& refTumorCount, int& altTumorCount,
				   std::string& tumorSample) const
{
   if (!BambinoParser::parseLine(line, chrName, position, variantType, ref, alt,
			         refCount, altCount))
      return false;

   StringVector value;
   getDelimitedStrings(line, '\t', value);

   refTumorCount = stringToInt(value[refTumorCountCol]);
   altTumorCount = stringToInt(value[altTumorCountCol]);

   if (refTumorCount < 0 || altTumorCount < 0)
      return false; // unable to convert strings to integers

   tumorSample = (tumorSampleCol >= 0 ? value[tumorSampleCol] : "");

   return true;
}

//------------------------------------------------------------------------------------
// SequenceTrie::SequenceTrie() initializes a sequence trie to represent an empty set
// of sequences

SequenceTrie::SequenceTrie()
   : endOfSequence(false)
{
   link[0] = NULL;
   link[1] = NULL;
   link[2] = NULL;
   link[3] = NULL;
}

//------------------------------------------------------------------------------------
// SequenceTrie::~SequenceTrie() de-allocates a sequence trie

SequenceTrie::~SequenceTrie()
{
   delete link[0];
   delete link[1];
   delete link[2];
   delete link[3];
}

//------------------------------------------------------------------------------------
// getLinkIndex() converts character A,C,G,T to 0,1,2,3; a value of -1 is returned if
// some other character is passed to this function

static inline int getLinkIndex(char base)
{
   switch (base)
   {
      case 'A': return  0;
      case 'C': return  1;
      case 'G': return  2;
      case 'T': return  3;
      default : return -1;
   }
}

//------------------------------------------------------------------------------------
// SequenceTrie::addSequence() adds a sequence to the set of sequences represented by
// the trie

void SequenceTrie::addSequence(const std::string& sequence)
{
   int seqlen = sequence.length();

   SequenceTrie *current = this;

   for (int i = 0; i < seqlen; i++)
   {
      int index = getLinkIndex(sequence[i]);

      if (index < 0)
         throw std::runtime_error("invalid sequence \"" + sequence + "\"");

      if (!current->link[index])
         current->link[index] = new SequenceTrie();

      current = current->link[index];
   }

   current->endOfSequence = true;
}

//------------------------------------------------------------------------------------
// SequenceTrie::findSequence() returns true if the given sequence is in the set of
// sequences represented by a trie

bool SequenceTrie::findSequence(const std::string& sequence) const
{
   int seqlen = sequence.length();

   const SequenceTrie *current = this;

   for (int i = 0; i < seqlen; i++)
   {
      int index = getLinkIndex(sequence[i]);

      if (index < 0 || !current->link[index])
         return false;

      current = current->link[index];
   }

   return current->endOfSequence;
}

//------------------------------------------------------------------------------------
// NumberSet::NumberSet() initializes a set of numbers to the empty set

NumberSet::NumberSet()
   : n(0), sum(0.0), sumsq(0.0)
{
}

//------------------------------------------------------------------------------------
// NumberSet::addNumber() adds the given number to the set

void NumberSet::addNumber(double x)
{
   if (++n == 1) // first number added to the set
      min = max = x;
   else if (x < min)
      min = x;
   else if (x > max)
      max = x;

   sum   += x;
   sumsq += x * x;
}

//------------------------------------------------------------------------------------
// NumberSet::average() computes and returns the average of the numbers in the set;
// a value of 0.0 is returned if the set is empty

double NumberSet::average() const
{
   return (n == 0 ? 0.0 : sum / n);
}

//------------------------------------------------------------------------------------
// NumberSet::variance() computes and returns the sample variance of the numbers in
// the set; a value of 0.0 is returned if the set contains fewer than two numbers

double NumberSet::variance() const
{
   if (n < 2)
      return 0.0;

   double avg = average();

   return ((sumsq - n * avg * avg) / (n - 1));
}

//------------------------------------------------------------------------------------
// NumberSet::stdev() computes and returns the sample standard deviation of the
// numbers in the set; a value of 0.0 is returned if the set contains fewer than two
// numbers

double NumberSet::stdev() const
{
   return std::sqrt(variance());
}

//------------------------------------------------------------------------------------
// ObservationSet::ObservationSet() initializes a set of observations to the empty set

ObservationSet::ObservationSet()
   : n(0), sumx(0.0), sumy(0.0), sumxx(0.0), sumyy(0.0), sumxy(0.0)
{
}

//------------------------------------------------------------------------------------
// ObservationSet::addObservation() adds the given observation to the set

void ObservationSet::addObservation(double x, double y)
{
   n++;

   sumx  += x;
   sumy  += y;

   sumxx += x * x;
   sumyy += y * y;
   sumxy += x * y;
}

//------------------------------------------------------------------------------------
// ObservationSet::pearsonCorrelationCoefficient() computes and returns Pearson's
// correlation coefficient for the set of observations; a value of 0.0 is returned if
// the coefficient is undefined

double ObservationSet::pearsonCorrelationCoefficient() const
{
   double d1 = n * sumxx - sumx * sumx;
   double d2 = n * sumyy - sumy * sumy;

   if (d1 <= 0.0 || d2 <= 0.0)
      return 0.0; // coefficient is undefined

   return (n * sumxy - sumx * sumy) / (std::sqrt(d1) * std::sqrt(d2));
}

//------------------------------------------------------------------------------------

struct SpearmanObservationCompareX
{
   bool operator()(const SpearmanObservation& a, const SpearmanObservation& b) const
   {
      return (a.x < b.x); // order by ascending x
   }
};

struct SpearmanObservationCompareY
{
   bool operator()(const SpearmanObservation& a, const SpearmanObservation& b) const
   {
      return (a.y < b.y); // order by ascending y
   }
};

//------------------------------------------------------------------------------------
// assignSpearmanRank() assigns rank values to xrank or yrank depending on the forX
// argument

static void assignSpearmanRank(SpearmanObservationVector& obs, bool forX)
{
   if (forX)
      std::sort(obs.begin(), obs.end(), SpearmanObservationCompareX());
   else
      std::sort(obs.begin(), obs.end(), SpearmanObservationCompareY());

   int n = obs.size();
   int i = 0;

   while (i < n)
   {
      double ranksum = i;
      int j = i + 1;

      while (j < n && (forX && obs[i].x == obs[j].x || !forX && obs[i].y == obs[j].y))
      {
         ranksum += j;
	 j++;
      }

      double rank = ranksum / (j - i);

      while (i < j)
         if (forX)
            obs[i++].xrank = rank;
         else
            obs[i++].yrank = rank;
   }
}

//------------------------------------------------------------------------------------
// SpearmanObservationSet::rankCorrelationCoefficient() computes and returns the
// Spearman rank correlation coefficient for the set of observations; a value of 0.0
// is returned if the coefficient is undefined

double SpearmanObservationSet::rankCorrelationCoefficient()
{
   int n = obs.size();
   if (n == 0)
      return 0.0; // coefficient is undefined

   assignSpearmanRank(obs, true);  // assign values to xrank
   assignSpearmanRank(obs, false); // assign values to yrank

   // Spearman rank correlation is Pearson correlation applied to the ranks

   ObservationSet pearsonObs;

   for (int i = 0; i < n; i++)
      pearsonObs.addObservation(obs[i].xrank, obs[i].yrank);

   return pearsonObs.pearsonCorrelationCoefficient();
}
