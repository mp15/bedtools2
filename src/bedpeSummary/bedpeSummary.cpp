/*****************************************************************************
  bedpeIntersect.cpp

  (c) 2025 - Martin Pollard
  Campbell Group
  Genome Research Ltd t/a Wellcome Sanger Institute
  mp15@sanger.ac.uk

  Licenced under the GNU General Public License 2.0 license.
******************************************************************************/
#include "lineFileUtilities.h"
#include "bedFilePE.h"
#include "GenomeFile.h"
#include "version.h"

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <experimental/iterator>
#include <cstdlib>
#include <cmath>


// define our program name
#define PROGRAM_NAME "bedpesummary"

// define our parameter checking macro
#define PARAMETER_CHECK(param, paramLen, actualLen) (strncmp(argv[i], param, min(actualLen, paramLen))== 0) && (actualLen == paramLen)

// function declarations
static void bedpesummary_help(void);
static void ProcessBedPE_summary(BedFilePE *bedpe);

int bedpesummary_main(int argc, char* argv[]) {
    // our configuration variables
    bool showHelp = false;

    // input files
    std::string bedpeFileA = "stdin";

    for(int i = 1; i < argc; i++) {
        int parameterLength = (int)strlen(argv[i]);

        if((PARAMETER_CHECK("-h", 2, parameterLength)) ||
        (PARAMETER_CHECK("--help", 5, parameterLength))) {
            showHelp = true;
        }
    }

    if(showHelp) bedpesummary_help();

    // do some parsing (all of these parameters require 2 strings)
    for(int i = 1; i < argc; i++) {

        int parameterLength = (int)strlen(argv[i]);

        if(PARAMETER_CHECK("-i", 2, parameterLength)) {
            if ((i+1) < argc) {
                bedpeFileA = argv[i + 1];
                i++;
            }
        }
        else {
            std::cerr << std::endl << "*****ERROR: Unrecognized parameter: " << argv[i] << " *****" << std::endl << std::endl;
            showHelp = true;
        }
    }

    if (!showHelp) {
        BedFilePE *bedpe= new BedFilePE(bedpeFileA);

        ProcessBedPE_summary(bedpe);
    }
    else {
        bedpesummary_help();
    }
    return 0;
}


static void bedpesummary_help(void) {
    std::cerr << "\nTool:    bedtools bedpesummary (aka bedpeSummary)" << std::endl;
    std::cerr << "Version: " << VERSION << std::endl;
    std::cerr << "Summary: Summarises a BEDPE file." << std::endl << std::endl;

    std::cerr << "Usage:   " << PROGRAM_NAME << " [OPTIONS] -i <bedpe>" << endl << endl;

    std::cerr << "Options: " << endl;

    // end the program here
    exit(1);
}

static inline CHRPOS calc_median_chrposv(std::vector<CHRPOS>& accum)
{
    std::sort(accum.begin(), accum.end());
    size_t n = accum.size();
    return n % 2 ? accum[n / 2] : (accum[n / 2 - 1] + accum[n / 2]) / 2;
}

class HistogramData
{
public:
    HistogramData(const std::vector<CHRPOS>& data, const int num_bins) :
        n_bins(num_bins), min_val(0), max_val(0), bin_width(0), bin_data(num_bins, 0)
    {
        if (data.empty() || num_bins <= 0) return;

        min_val = *std::min_element(data.begin(), data.end());
        max_val = *std::max_element(data.begin(), data.end());
        bin_width = (max_val - min_val) / num_bins;
    
        // Protect against zero-width bins
        if (bin_width == 0) return;
   
        for (CHRPOS value : data) {
            int bin_idx = static_cast<int>((value - min_val) / bin_width);
            // Handle edge case: if value == max, put it in the last bin
            if (bin_idx == num_bins) bin_idx--;
            bin_data[bin_idx]++;
        }
    }
public:
    int n_bins;
    CHRPOS min_val;
    CHRPOS max_val;
    CHRPOS bin_width;

    std::vector<int> bin_data;
};

static void ProcessBedPE_summary(BedFilePE *bedpe) {
    BEDPE bedpeEntry_a;
    BedLineStatus bedpeStatus_a;
    int lineNum_a = 0;
    // totals
    int n_intrachrom = 0, n_interchrom = 0, inversion = 0, insertion = 0, deletion = 0;
    CHRPOS total_distance = 0;
    std::vector<CHRPOS> accum;

    // open the BEDPE file for reading.
    bedpe->Open();
    // Prime the pump
    bedpeStatus_a = bedpe->GetNextBedPE(bedpeEntry_a, lineNum_a);
    if (bedpeStatus_a == BED_INVALID) {bedpe->Close(); return;}

    do {
        // Ensure we have a record and not a header line etc.
        if (bedpeStatus_a == BED_VALID) {
            if (bedpeEntry_a.chrom1 == bedpeEntry_a.chrom2) {
                n_interchrom++;
            } else {
                n_intrachrom++;
                CHRPOS distance = llabs(bedpeEntry_a.start2-bedpeEntry_a.start1);
                accum.push_back(distance);
                total_distance += distance;
                // Intrachromosomal
                if (bedpeEntry_a.strand1 == bedpeEntry_a.strand2) {
                    inversion++;
                } else if (bedpeEntry_a.strand1 == "+" && bedpeEntry_a.strand2 == "-") {
                    deletion++;
                } else if (bedpeEntry_a.strand1 == "-" && bedpeEntry_a.strand2 == "+") {
                    insertion++;
                }
            }
        }
    } while ((bedpeStatus_a = bedpe->GetNextBedPE(bedpeEntry_a, lineNum_a)) != BED_INVALID);
    // calc mean len
    CHRPOS mean_len = ((n_intrachrom != 0) ? total_distance / n_intrachrom : std::nan(""));
    // calc median len
    CHRPOS median_len = calc_median_chrposv(accum);
    std::cout << "{";
    std::cout << "\"inversion\" : " << inversion << ", \"insertion\" : " << insertion
              << ", \"deletion\" : " << deletion  << ", "<< std::endl;
    std::cout << "\"n_interchrom\" : " << n_interchrom << ", \"n_intrachrom\" : " << n_intrachrom << ", ";

    std::cout << "\"mean intrachromasomal sv length\" : " << mean_len  << ", "<< std::endl;
    std::cout << "\"median intrachromasomal sv length\" : " << median_len << ", "<< std::endl;
    HistogramData hist(accum, 10);
    std::cout << "\"histogram\" : { \"min_val\" : " << hist.min_val << ", ";
    std::cout << "\"bin_width\" : " << hist.bin_width << ", \"bin_counts\": ["<< std::endl;

    std::copy(std::begin(hist.bin_data),
          std::end(hist.bin_data),
          std::experimental::make_ostream_joiner(std::cout, ", "));

    std::cout << "]}}" << std::endl;

    bedpe->Close();
}
