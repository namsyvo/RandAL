#include <cstring>
#include <algorithm>
#include <sstream>
#include <math.h>

#include "ReadAlignment.hpp"

using namespace std;

class Parameter
{
    public:
        string infoRead;
        string nameRead;
        string CIGAR;
        string qualRead;
        int flag;
};

void usage();
void printParameters();
void parse_parameters (int argc, char **argv);

string reverseRead(string& rea);
template <class T> string to_string1(T t);
void process(string& r, Parameter& pa);
void alignReads(const string& r, const Parameter& p, uint64_t pos);

int length;
float error_prob = -1;
int att_para = -1;
int th_para = -1;

int mode = -1;
int th = -1;
int att = -1;
int rep = -1;
int minws = -1;
int maxcdt = -1;

char *backwardIndexFile = "";
char *forwardIndexFile = "";
char *refFile = "";
char *readFile = "";
char *mapFile = "";
uint32_t num = 0;
ReadAlignment ra;

ifstream infile;
ofstream outfile;

int main(int argc, char **argv)
{
    if (argc < 4)
    {
        usage();
    }
    else
    {
        parse_parameters(argc, argv);
    }

    srand(time(NULL)); //initialize the random generator

    infile.open(readFile);     //Read the read file (*.fq)
    outfile.open(mapFile);     //Open the map file (*.sam)
    string r = "";
    Parameter para;
    getline(infile, r);     //The 1st read
    para.infoRead = r.erase(0, 1);
    para.nameRead = r.substr(0, r.find("_"));
    string read1st = "";
    getline(infile, read1st);

    //Set parameters
    length = read1st.length();
    if (error_prob == -1) error_prob = 0.02;
    if (att_para == -1) att_para = 1;
    if (th_para == -1) th_para = 4;

    if (mode == -1) mode = 3;
    if (rep == -1) rep = 10;
    if (maxcdt == -1) maxcdt = 2;

    if (th == -1) th = int(floor(error_prob * length + th_para * sqrt(length * error_prob * (1 - error_prob))));
    if (att == -1) att = att_para * (th + 1);
    if (minws == -1) minws = int(ceil(length/(th + 1))) + 3;

    printParameters();
    ra.Init(mode, rep, th, att, minws, backwardIndexFile, forwardIndexFile, refFile);

    para.CIGAR = to_string1(read1st.length()) + "M";
    getline(infile, r); 
    getline(infile, para.qualRead);

    process(read1st, para);

    //For each line (read) find the match positions and store those positions to output file
    while (getline(infile, r))
    {       
        if (r[0] == '@')
        {
            para.infoRead = r.erase(0, 1);  
            getline(infile, r);                                 
            process(r, para);
        }           
    }
    infile.close();
    outfile.close();
    cout << "Finish aligning reads, number of mapped reads: " << num << endl;
}

//
void process(string& r, Parameter& pa)
{
    vector<size_t> positions;
    vector<uint16_t> distances; 
    ra.AlignRead(r, positions, distances);    
    string rr;
    size_t pos=0;
    int dis;

    if (positions.empty())
    {
        rr = reverseRead(r);
        ra.AlignRead(rr, positions, distances);
        if (positions.empty())
        {
            pa.flag = 4;
            alignReads(rr, pa, pos);
        } 
        else
        {
            dis = distances[0];
            pos = positions[0];
			if (positions.size() > maxcdt)
				return;
            for (unsigned int i = 0 ; i < positions.size() ; i++)
            {
                if (distances[i] == 0)
                {
                    pos = positions[i];
                    break;
                }
                if (dis > distances[i])
                {
                    dis = distances[i];
                    pos = positions[i];
                }
            }
            pa.flag = 16;
            alignReads(rr, pa, pos);
            num++;
        }
    }
    else
    {
        dis = distances[0];
        pos = positions[0];
	    if (positions.size() > maxcdt)
			return;
        for (unsigned int i=0; i<positions.size(); i++)
        {
            if (distances[i] == 0)
            {
                pos = positions[i];
                break;
            }
            if (dis > distances[i])
            {
                dis = distances[i];
                pos = positions[i];
            }
        }
        pa.flag = 0;
        alignReads(r, pa, pos);
        num++;
     }
}

//
template <class T> string to_string1(T t)
{
    stringstream ss;
    ss << t;
    return ss.str();
}

//
inline string reverseRead(string& rea)
{
    string t = "";  
    for (size_t i=0; i<rea.length(); i++)
    {
        if (rea[i] =='A') t = t+'T'; else
        if (rea[i] =='T') t = t+'A'; else
        if (rea[i] =='G') t = t+'C'; else
        if (rea[i] =='C') t = t+'G'; else
        t = t+rea[i];
    }
    reverse(t.begin(), t.end());
    return t;
}

//
inline void alignReads(const string& r, const Parameter& p, uint64_t pos)
{
    if (p.flag == 4)
    {
        outfile << p.infoRead << "\t" << p.flag << "\t*\t0\t0\t*\t*\t0\t0\t";
        outfile << r << "\t" << p.qualRead << endl;
    }
    else
    {
        outfile << p.infoRead << "\t" << p.flag << "\t" << p.nameRead << "\t" << pos << "\t" << "255\t" << p.CIGAR << " \t*\t0\t0\t";
        outfile << r << "\t" << p.qualRead << endl;
    }
}

//
void usage()
{
    cerr << endl
    << "Usage:  randal IndexFileName ReadFileName MapFileName [OPTION]" << endl
    << "        IndexFileName           : index file name (from preprocessing step with rand-build)" << endl
    << "        ReadFileName            : read file name (in FASTQ format)" << endl
    << "        MapFileName             : file to store alignment result (in SAM format)" << endl
    << "        [OPTION]            : a list of zero or more optional arguments" << endl
    << "                -mode number    : distance mode (2:hamming distance; 3:edit distance; default 3)" << endl
    << "                -th number      : threshould for distance (integer, >= 0)" << endl
    << "                -att number     : number of attempts for the randomization (integer, > 0)" << endl
    << "                -minws number   : minimal length of common substring between read and reference (integer, > 0, < read length)" << endl
    << "                -maxcdt number  : maximal number of candidates for alignment (integer, > 0)" << endl;
    exit(0);
}

//
void printParameters()
{
    cout << "\nParameters:\n"
         << "Distance mode: " << mode << ", Distance threshold: " << th
         << ", Number of atempts: " << att << ", Minimal LCS: " << minws << ", Maximal alignment: " << maxcdt
         << ", Backward Index File: " << backwardIndexFile << ", Forward Index File: " << forwardIndexFile
         << ", Reference File: " << refFile << ", Read File: " << readFile << ", Map File:" << mapFile << endl;
}

//
void parse_parameters (int argc, char **argv)
{
    char *indexFile;
    // the first three arguments shouldn't start with minus
    if(argv[1][0] == '-' ||argv[2][0] == '-'||argv[3][0] == '-')
    {
        usage();
        return;
    }
    indexFile   = argv[1];
    readFile    = argv[2];
    mapFile     = argv[3];

    char *ref_ext = ".ref";
    char *bw_ext = ".bw";
    char *fw_ext = ".fw";
    int len = strlen(indexFile);
    refFile = new char[len + 3];
    backwardIndexFile = new char[len + 2];
    forwardIndexFile = new char[len + 2];
    strcpy(refFile, indexFile);
    strcat(refFile, ref_ext);
    strcpy(backwardIndexFile, indexFile);
    strcat(backwardIndexFile, bw_ext);
    strcpy(forwardIndexFile, indexFile);
    strcat(forwardIndexFile, fw_ext);

    int argno;
    for (argno = 4; argno < argc - 1; argno++)
    {
        if (argv[argno][0] == '-')
        {
            if (!strcmp (argv[argno], "-mode"))
            {
                if (argno == argc - 1) std::cerr << "Must specify searching mode after -mode" << std::endl;
                mode = atoi(argv[++argno]);
                if ((mode != 1) && (mode != 2) && (mode != 3))
                {
                    usage();
                    return;
                }
            }
            else if (!strcmp (argv[argno], "-th"))
            {
                if (argno == argc - 1) std::cerr << "Must specify threshold after -th" << std::endl;
                th = atoi(argv[++argno]);
                if (th < 0)
                {
                    usage();
                    return;
                }
            }
            else if (!strcmp (argv[argno], "-att"))
            {
                if (argno == argc - 1) std::cerr << "Must specify number of attempts after -att" << std::endl;
                att = atoi(argv[++argno]);
                if (att <= 0)
                {
                    usage();
                    return;
                }
            }
            else if (!strcmp (argv[argno], "-minws"))
            {
                if (argno == argc - 1) std::cerr << "Must specify minimal length of common substring between read and reference after -minws" << std::endl;
                minws = atoi(argv[++argno]);
                if (minws <= 0)
                {
                    usage();
                    return;
                }
            }
            else if (!strcmp (argv[argno], "-maxcdt"))
            {
                if (argno == argc - 1) std::cerr << "Must specify maximal number of candidates for alignment after -maxcdt" << std::endl;
                maxcdt = atoi(argv[++argno]);
                if (maxcdt <= 0)
                {
                    usage();
                    return;
                }
            }
        }
    } 
}
