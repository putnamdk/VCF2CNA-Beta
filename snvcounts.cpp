//------------------------------------------------------------------------------------
//
// snvcounts.cpp - program that extracts mutant and total counts for SNVs in tumor and
//                 normal samples from a Bambino output file ("high_20") or a file in
//                 the Mutation Annotation Format (MAF); the output files written by
//                 this program are inputs to the consprep program
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2016 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "genutil.h"

const uint16_t MAX_COUNT = 65535;

uint64_t occurrences;                     // number of normal coverage values
uint64_t normalTotalCount[MAX_COUNT + 1]; // histogram of normal coverage values

void compressCounts(int inMutant, int inTotal, uint16_t& outMutant,
		    uint16_t& outTotal); // forward declaration

//------------------------------------------------------------------------------------

class PosCounts // concisely stores the counts for one sample
{
public:
   PosCounts(int inTumorMutant, int inTumorTotal, int inNormalMutant,
	     int inNormalTotal)
   {
      compressCounts(inTumorMutant,  inTumorTotal,  tumorMutant,  tumorTotal);
      compressCounts(inNormalMutant, inNormalTotal, normalMutant, normalTotal);
   }

   uint16_t tumorMutant, tumorTotal, normalMutant, normalTotal;
};

typedef std::map<int, PosCounts> PosMap; // key is position
PosMap posmap[NUM_CHROMOSOMES + 1];      // one map for each chromosome

//------------------------------------------------------------------------------------

class MAF_Parser // for parsing lines in Mutation Annotation Format (MAF)
{
public:
   MAF_Parser(const std::string& headingLine);
   virtual ~MAF_Parser() { }

   virtual bool parseLine(const std::string& line, std::string& chrName,
		          int& position, std::string& variantType,
			  int& tumorMutant, int& tumorTotal,
			  int& normalMutant, int& normalTotal) const;

   //column numbers of columns of interest
   int chrCol, posCol, typeCol, tumorMutantCol, tumorTotalCol, normalMutantCol,
       normalTotalCol;

   int numColumns;
};

//------------------------------------------------------------------------------------
// MAF_Parser::MAF_Parser() determines the column number of columns of interest by
// parsing a heading line read from a MAF file

MAF_Parser::MAF_Parser(const std::string& headingLine)
   : chrCol(-1), posCol(-1), typeCol(-1), tumorMutantCol(-1), tumorTotalCol(-1),
     normalMutantCol(-1), normalTotalCol(-1)
{
   StringVector heading;
   getDelimitedStrings(headingLine, '\t', heading);

   numColumns = heading.size();

   for (int i = 0; i < numColumns; i++)
   {
      std::string h = heading[i];

      if (h == "Chromosome")                                   chrCol          = i;
      else if (h == "Start_Position" || h == "Start_position") posCol          = i;
      else if (h == "Variant_Type"   || h == "VariantType")    typeCol         = i;
      else if (h == "Tumor_ReadCount_Alt")                     tumorMutantCol  = i;
      else if (h == "Tumor_ReadCount_Total")                   tumorTotalCol   = i;
      else if (h == "Normal_ReadCount_Alt")                    normalMutantCol = i;
      else if (h == "Normal_ReadCount_Total")                  normalTotalCol  = i;
   }

   if (chrCol < 0 || posCol < 0 || typeCol < 0 || tumorMutantCol < 0 ||
       tumorTotalCol < 0 || normalMutantCol < 0 || normalTotalCol < 0)
      throw std::runtime_error("missing column(s) in MAF file");
}

//------------------------------------------------------------------------------------
// MAF_Parser::parseLine() parses a non-heading line read from a MAF file

bool MAF_Parser::parseLine(const std::string& line, std::string& chrName,
		           int& position, std::string& variantType,
			   int& tumorMutant, int& tumorTotal,
			   int& normalMutant, int& normalTotal) const
{
   StringVector value;
   getDelimitedStrings(line, '\t', value);

   if (value.size() != numColumns)
      return false;

   position     = stringToInt(value[posCol]);
   tumorMutant  = stringToInt(value[tumorMutantCol]);
   tumorTotal   = stringToInt(value[tumorTotalCol]);
   normalMutant = stringToInt(value[normalMutantCol]);
   normalTotal  = stringToInt(value[normalTotalCol]);

   if (position < 0 || tumorMutant < 0 || tumorTotal < 0 || normalMutant < 0 ||
       normalTotal < 0)
      return false; // unable to convert strings to integers

   chrName      = value[chrCol];
   variantType  = value[typeCol];

   return true;
}

//------------------------------------------------------------------------------------
// compressCounts() converts counts from four-byte signed integers to two-byte
// unsigned integers

void compressCounts(int inMutant, int inTotal, uint16_t& outMutant,
		    uint16_t& outTotal)
{
   if (inMutant > inTotal)
      inMutant = inTotal;

   if (inTotal > MAX_COUNT) // count is too large to fit in two bytes
   {
      // proportionally reduce the mutant count so that the ratio of mutant/total
      // is preserved
      inMutant = roundit(MAX_COUNT * static_cast<double>(inMutant) / inTotal);
      inTotal  = MAX_COUNT;
   }

   outMutant = static_cast<uint16_t>(inMutant);
   outTotal  = static_cast<uint16_t>(inTotal);
}

//------------------------------------------------------------------------------------
// readFile() reads a Bambino output file or MAF file and stores the position data
// in an array of maps

void readFile(const std::string& filename)
{
   BambinoParserTumor *bp = NULL;
   MAF_Parser         *mp = NULL;

   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw std::runtime_error("unable to open " + filename);

   std::string line;

   if (!std::getline(infile, line))
      throw std::runtime_error("empty file " + filename);

   // examine the heading line to see what kind of file it is

   try
   {
      bp = new BambinoParserTumor(line);
   }
   catch (const std::runtime_error&) { }

   if (!bp) // it is not a Bambino output file
   {
      try
      {
         mp = new MAF_Parser(line);
      }
      catch (const std::runtime_error&) { }
   }

   if (!bp && !mp)
      throw std::runtime_error("unrecognized file format in " + filename);

   // now read the file

   while (std::getline(infile, line))
   {
      std::string chrName, type, ref, alt, tumorSample;
      int position, tumorRef,  tumorMutant,  tumorTotal,
	            normalRef, normalMutant, normalTotal;

      if (bp && bp->parseLine(line, chrName, position, type, ref, alt, normalRef,
			      normalMutant, tumorRef, tumorMutant, tumorSample))
      {
         tumorTotal  = tumorRef  + tumorMutant;
	 normalTotal = normalRef + normalMutant;
      }
      else if (mp && mp->parseLine(line, chrName, position, type, tumorMutant,
			           tumorTotal, normalMutant, normalTotal))
      {
      }
      else
         throw std::runtime_error("unable to parse line in " + filename + " \"" +
			          line + "\"");

      int chrnum;
      if (type != "SNP" || (chrnum = getChrNumber(chrName)) == 0)
         continue; // this is not an SNV or this is an unrecognized chromosome

      PosMap& pmap = posmap[chrnum];         // get the map for this chromosome
      if (pmap.find(position) == pmap.end()) // this position is not already in map
      {
         pmap.insert(std::make_pair(position,
	             PosCounts(tumorMutant, tumorTotal, normalMutant, normalTotal)));

	 occurrences++;

	 if (normalTotal > MAX_COUNT)
            normalTotalCount[MAX_COUNT]++;
	 else
	    normalTotalCount[normalTotal]++;
      }
   }

   if (infile.bad())
      throw std::runtime_error("read error in " + filename);

   infile.close();

   delete bp;
   delete mp;
}

//------------------------------------------------------------------------------------
// writeCounts() writes the counts in order by chromosome and position

void writeCounts(const std::string& filename)
{
   std::ofstream outfile(filename.c_str());
   if (!outfile.is_open())
      throw std::runtime_error("unable to open " + filename);

   outfile << "Chr"
	   << "\t" << "Pos"
	   << "\t" << "TumorMutant"
	   << "\t" << "TumorTotal"
	   << "\t" << "NormalMutant"
	   << "\t" << "NormalTotal"
	   << std::endl;

   for (int chrnum = 1; chrnum <= NUM_CHROMOSOMES; chrnum++)
   {
      PosMap& pmap = posmap[chrnum];

      for (PosMap::iterator ppos = pmap.begin(); ppos != pmap.end(); ++ppos)
         outfile << chrLongName[chrnum]
		 << "\t" << ppos->first
		 << "\t" << ppos->second.tumorMutant
		 << "\t" << ppos->second.tumorTotal
		 << "\t" << ppos->second.normalMutant
		 << "\t" << ppos->second.normalTotal
		 << "\n";
   }

   outfile.close();
}

//------------------------------------------------------------------------------------
// writeMedian() computes and writes the median normal coverage

void writeMedian(const std::string& filename)
{
   uint64_t half  = (occurrences + 1) / 2;
   uint64_t count = normalTotalCount[0];
   int i = 0;

   while (count < half)
      count += normalTotalCount[++i];

   std::ofstream outfile(filename.c_str());
   if (!outfile.is_open())
      throw std::runtime_error("unable to open " + filename);

   outfile << i << std::endl; // write the median normal coverage

   outfile.close();
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   if (argc != 4)
   {
      std::cout << "Usage: " << argv[0]
	        << " inputfile"
		<< " snvcounts_outputfile"
		<< " median_outputfile"
		<< std::endl;
      return 1;
   }

   try
   {
      std::string infilename  = argv[1];
      std::string cntfilename = argv[2];
      std::string medfilename = argv[3];

      readFile(infilename);
      writeCounts(cntfilename);
      writeMedian(medfilename);
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
