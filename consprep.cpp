//------------------------------------------------------------------------------------
//
// consprep.cpp - program that prepares data files needed by CONSERTING
//
// Author: Stephen V. Rice, Ph.D.
//
// Copyright 2016 St. Jude Children's Research Hospital
//
//------------------------------------------------------------------------------------

#include "genutil.h"

// command-line option variables and default values

const double DEFAULT_MEDIAN     = 30.00;
const double DEFAULT_MINFACTOR  =  0.50;
const double DEFAULT_MAXFACTOR  =  1.50;
const double DEFAULT_XMINFACTOR =  0.25;
const double DEFAULT_XMAXFACTOR =  1.50;

double median     = DEFAULT_MEDIAN;
double minfactor  = DEFAULT_MINFACTOR;
double maxfactor  = DEFAULT_MAXFACTOR;
double xminfactor = DEFAULT_XMINFACTOR;
double xmaxfactor = DEFAULT_XMAXFACTOR;

// one set for each chromosome holds the bad positions in that chromosome; these sets
// are initialized from data read from the goodbad_file
std::set<int> badlist[NUM_CHROMOSOMES + 1];

// one integer for each chromosome indicates the number of consecutive,
// non-overlapping 100-bp windows in that chromosome; the data in this array is
// read from the wincount_file
int numWindows[NUM_CHROMOSOMES + 1];

// filename suffixes and file handles for the output files

const char *AI_FILENAME_SUFFIX = ".ai";
std::ofstream *aifile; // allelic imbalance output file

const char *CHR_FILENAME_SUFFIX = "_%s_100";
std::ofstream *chrfile[NUM_CHROMOSOMES + 1]; // one output file for each chromosome

//------------------------------------------------------------------------------------

class PosData // data associated with a particular position within a chromosome
{
public:
   PosData(int inChrnum, int inPosition, int inTumorMutant, int inTumorTotal,
	   int inNormalMutant, int inNormalTotal)
      : chrnum(inChrnum), position(inPosition),
	tumorMutant(inTumorMutant), tumorTotal(inTumorTotal),
	normalMutant(inNormalMutant), normalTotal(inNormalTotal),
	window(inPosition / 100) { }

   int chrnum, position, tumorMutant, tumorTotal, normalMutant, normalTotal, window;
};

//------------------------------------------------------------------------------------
// processOptions() processes the command-line arguments; false is returned if any of
// the arguments are invalid

bool processOptions(int argc, char *argv[],
		    std::string& goodbad_filename,
		    std::string& wincount_filename,
		    std::string& output_filenamePrefix)
{
   int n = 0; // number of non-option arguments found

   for (int i = 1; i < argc; i++)
   {
      std::string s = argv[i];
      if (s.length() == 0)
         return false;

      if (s[0] == '-') // found an option
      {
         StringVector part;
	 getDelimitedStrings(s, '=', part);

	 if (!(part.size() == 2 &&
              (part[0] == "-median"     && (median     = stringToDbl(part[1])) >= 0 ||
               part[0] == "-minfactor"  && (minfactor  = stringToDbl(part[1])) >= 0 ||
               part[0] == "-maxfactor"  && (maxfactor  = stringToDbl(part[1])) >= 0 ||
               part[0] == "-xminfactor" && (xminfactor = stringToDbl(part[1])) >= 0 ||
               part[0] == "-xmaxfactor" && (xmaxfactor = stringToDbl(part[1])) >= 0)))
            return false;
      }
      else
         switch (++n)
	 {
            case  1: goodbad_filename      = s; break;
            case  2: wincount_filename     = s; break;
            case  3: output_filenamePrefix = s; break;
            default: return false;
	 }
   }

   return (n == 3 && minfactor <= maxfactor && xminfactor <= xmaxfactor);
}

//------------------------------------------------------------------------------------
// showOption() displays one command-line option

void showOption(std::string optname, std::string description, double defaultValue)
{
   std::printf("  %s\t%s, default is %4.2f\n", optname.c_str(), description.c_str(),
               defaultValue);
}

//------------------------------------------------------------------------------------
// showUsage() displays the command-line usage for this program

void showUsage(const char *progname)
{
   std::cout << "Usage: " << progname
	     << " [OPTION ...] goodbad_file wincount_file"
	     << " output_path_prefix < snvcounts_file"
	     << std::endl << std::endl;

   showOption("-median=N",     "median normal coverage",         DEFAULT_MEDIAN);
   showOption("-minfactor=N",  "minimum scale factor, non-chrX", DEFAULT_MINFACTOR);
   showOption("-maxfactor=N",  "maximum scale factor, non-chrX", DEFAULT_MAXFACTOR);
   showOption("-xminfactor=N", "minimum scale factor, chrX",     DEFAULT_XMINFACTOR);
   showOption("-xmaxfactor=N", "maximum scale factor, chrX",     DEFAULT_XMAXFACTOR);
}

//------------------------------------------------------------------------------------
// readGoodBadList() reads a file containing SNVs that have been designated as
// SuperGood or SuperBad; the positions of bad SNVs are saved in a data structure

void readGoodBadList(const std::string& filename)
{
   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw std::runtime_error("unable to open " + filename);

   std::string line;

   while (std::getline(infile, line))
   {
      StringVector column;
      getDelimitedStrings(line, '\t', column);

      if (column.size() != 2)
         throw std::runtime_error("unexpected #columns in line of " + filename +
			          " \"" + line + "\"");

      if (column[1] != "SuperBad")
         continue; // ignore heading line and SuperGood lines

      try
      {
         Variant variant(column[0]); // parse the SNV string

	 // the SNV string is valid
	 int chrnum   = variant.chrNumber;
	 int position = variant.position;

	 // add this bad position if it has not already been saved
	 if (badlist[chrnum].find(position) == badlist[chrnum].end())
	    badlist[chrnum].insert(position);
      }
      catch (const std::runtime_error&)
      {
         throw std::runtime_error("invalid variant specification in " + filename +
			          " \"" + column[0] + "\"");
      }
   }

   if (infile.bad())
      throw std::runtime_error("read error in " + filename);

   infile.close();
}

//------------------------------------------------------------------------------------
// readNumWindows() reads a file containing the number of 100-bp windows in each
// chromosome; this data is saved in the numWindows array

void readNumWindows(const std::string& filename)
{
   std::ifstream infile(filename.c_str());
   if (!infile.is_open())
      throw std::runtime_error("unable to open " + filename);

   std::string line;

   while (std::getline(infile, line))
   {
      StringVector column;
      getDelimitedStrings(line, '\t', column);

      int chrnum = getChrNumber(column[0]);
      if (chrnum == 0)
         continue; // ignore heading line and unrecognized chromosomes

      if (column.size() != 2)
         throw std::runtime_error("unexpected #columns in line of " + filename +
			          " \"" + line + "\"");

      numWindows[chrnum] = stringToInt(column[1]);
   }

   if (infile.bad())
      throw std::runtime_error("read error in " + filename);

   infile.close();

   for (int chrnum = 1; chrnum <= NUM_CHROMOSOMES; chrnum++)
      if (numWindows[chrnum] <= 0)
         throw std::runtime_error("invalid or missing #windows for " +
			          chrLongName[chrnum] + " in " + filename);
}

//------------------------------------------------------------------------------------
// createOutputFiles() creates the output files and writes a heading line to each

void createOutputFiles(const std::string& filenamePrefix)
{
   std::string filename = filenamePrefix + AI_FILENAME_SUFFIX;

   aifile = new std::ofstream(filename.c_str());
   if (!aifile->is_open())
      throw std::runtime_error("unable to open " + filename);

   *aifile << "Chr"
	   << "\t" << "Pos"
	   << "\t" << "AIDiff"
	   << "\t" << "BAFT"
	   << "\t" << "BAFN"
	   << std::endl;

   for (int chrnum = 1; chrnum <= NUM_CHROMOSOMES; chrnum++)
   {
      char suffix[100];
      std::sprintf(suffix, CHR_FILENAME_SUFFIX, chrLongName[chrnum].c_str());

      filename = filenamePrefix + suffix;

      chrfile[chrnum] = new std::ofstream(filename.c_str());
      if (!chrfile[chrnum]->is_open())
         throw std::runtime_error("unable to open " + filename);

      *chrfile[chrnum] << "Dcvg"
	               << "\t" << "Gcvg"
		       << std::endl;
   }
}

//------------------------------------------------------------------------------------
// closeOutputFiles() closes all of the output files

void closeOutputFiles()
{
   aifile->close();

   for (int chrnum = 1; chrnum <= NUM_CHROMOSOMES; chrnum++)
      chrfile[chrnum]->close();
}

//------------------------------------------------------------------------------------
// readNextPosition() reads the next line from stdin and returns a pointer to a newly
// allocated PosData object containing the data in the line, or returns NULL if EOF
// has been reached

PosData *readNextPosition()
{
   std::string line;

   while (std::getline(std::cin, line))
   {
      StringVector column;
      getDelimitedStrings(line, '\t', column);

      if (column.size() != 6)
         throw std::runtime_error("unexpected #columns in line read from stdin \"" +
			          line + "\"");

      int chrnum = getChrNumber(column[0]);
      if (chrnum == 0)
         continue; // ignore heading line and unrecognized chromosomes

      int position     = stringToInt(column[1]);
      int tumorMutant  = stringToInt(column[2]);
      int tumorTotal   = stringToInt(column[3]);
      int normalMutant = stringToInt(column[4]);
      int normalTotal  = stringToInt(column[5]);

      if (position < 0 || tumorMutant < 0 || tumorTotal < 0 || normalMutant < 0 ||
	  normalTotal < 0)
         throw std::runtime_error("invalid data in line read from stdin \"" + line +
			          "\"");

      if (tumorMutant > tumorTotal)
         tumorMutant = tumorTotal;

      if (normalMutant > normalTotal)
         normalMutant = normalTotal;

      return new PosData(chrnum, position, tumorMutant, tumorTotal, normalMutant,
		         normalTotal);
   }

   return NULL; // reached EOF
}

//------------------------------------------------------------------------------------
// processWindow() processes positions that fall in a particular window and writes to
// the chromosome file the average tumor coverage and average normal coverage of these
// positions; note that positions not in chrX that have a bad SNV are excluded from
// the computation of average coverage; positions with normal coverage below the
// minimum or above the maximum are also excluded; this function returns the first
// position not in the current window, or returns NULL if EOF has been reached

PosData *processWindow(PosData *pd)
{
   const int chrX = 23;
   const double EPSILON = 0.0001; // to avoid division by zero

   int count          = 0;
   int sumTumorTotal  = 0;
   int sumNormalTotal = 0;

   int chrnum = pd->chrnum;
   int window = pd->window;

   double minCoverage = median * (chrnum == chrX ? xminfactor : minfactor);
   double maxCoverage = median * (chrnum == chrX ? xmaxfactor : maxfactor);

   while (pd && chrnum == pd->chrnum && window == pd->window)
   {
      if (chrnum == chrX ||
          badlist[chrnum].find(pd->position) == badlist[chrnum].end())
      {
         double normalMAF = pd->normalMutant / (pd->normalTotal + EPSILON);

	 if (pd->tumorTotal > 15 && pd->normalTotal > 15 &&
	     normalMAF > 0.4 && normalMAF < 0.6)
	 {
            double tumorMAF = pd->tumorMutant / (pd->tumorTotal + EPSILON);

	    char buffer[100];

	    std::sprintf(buffer, "%s\t%d\t%.2f\t%.2f\t%.2f\n",
			 chrLongName[chrnum].c_str(), pd->position,
			 std::abs(tumorMAF - normalMAF), tumorMAF, normalMAF);

	    *aifile << buffer;
	 }

	 if (pd->normalTotal >= minCoverage && pd->normalTotal <= maxCoverage)
	 {
            count++;
	    sumTumorTotal  += pd->tumorTotal;
	    sumNormalTotal += pd->normalTotal;
	 }
      }

      delete pd;
      pd = readNextPosition();
   }

   *chrfile[chrnum] << roundit(sumTumorTotal  / (count + EPSILON)) << "\t"
		    << roundit(sumNormalTotal / (count + EPSILON)) << "\n";

   return pd;
}

//------------------------------------------------------------------------------------
// processAllChromosomes() reads position data from stdin and writes one line for each
// window in each chromosome giving the average tumor coverage and average normal
// coverage of positions in that window

void processAllChromosomes()
{
   PosData *pd = readNextPosition(); // read first position

   for (int chrnum = 1; chrnum <= NUM_CHROMOSOMES; chrnum++)
   {
      int wincount = numWindows[chrnum];

      for (int window = 0; window < wincount; window++)
         if (pd && chrnum == pd->chrnum && window == pd->window)
            pd = processWindow(pd); // process one or more positions in this window
         else // no positions in this window
	    *chrfile[chrnum] << "0\t0\n"; // average coverage is zero
   }

   if (pd)
      throw std::runtime_error("lines read from stdin are invalid or unsorted");
}

//------------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
   std::string goodbad_filename, wincount_filename, output_filenamePrefix;

   if (!processOptions(argc, argv, goodbad_filename, wincount_filename,
		       output_filenamePrefix))
   {
      showUsage(argv[0]);
      return 1;
   }

   try
   {
      readGoodBadList(goodbad_filename);
      readNumWindows(wincount_filename);

      createOutputFiles(output_filenamePrefix);
      processAllChromosomes();
      closeOutputFiles();
   }
   catch (const std::runtime_error& error)
   {
      std::cerr << argv[0] << ": " << error.what() << std::endl;
      return 1;
   }

   return 0;
}
