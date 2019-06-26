// C++ csv/tsv data generation companion
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
//#include <functional>
//#include <utility>
#include <map>
#include <math.h>
#include <algorithm>
#include <vector>
using namespace std;


int barcodefinder(string FileName, int length)
{
  //string FileName;
  int barcodelines;
  int barcodelengthlim;
  int barcodelengththousand = 0;
  int barcode;

  map<char, char> rev;
  rev['A']='T';
  rev['T']='A';
  rev['C']='G';
  rev['G']='C';

  ifstream infile(FileName);
  //infile.open(FileName);
  string line;
  if (infile.is_open())
  {
    while (getline(infile, line))
    {
      if (line[0] != '>' && line.length() == length)
      {
        while (barcodelines < 1000)
        {
          for (barcodelengthlim=0; barcodelengthlim<8; barcodelengthlim++)
          {
            if (line[barcodelengthlim] == rev[line[line.length()-(barcodelengthlim+1)]])
            {
              barcodelengththousand += 1;
            }
            else{break;}
          }
          ++barcodelines;
        }
      }
    }
  }
  //cout << barcodelengththousand << '\n';
  barcode = round(barcodelengththousand/1000);
  return barcode;
  //cout << barcode;
}


int main()
{
  string FileName;
  cout << "File name: ";
  cin >> FileName;
  int length;
  cout << "Read lengths: ";
  cin >> length;
  int mink;
  cout << "Min k: ";
  cin >> mink;
  int maxk;
  cout << "Max k: ";
  cin >> maxk;

  int barcode;
  int numoflines = 0;
  int pos;

  cout << FileName << ", " << length << ", " << mink << ", " << maxk << endl;

  char bases[]{'A', 'T', 'C', 'G'};
  string basesstr = "ATCG";

  map<char, char> rev;
  rev['A']='T';
  rev['T']='A';
  rev['C']='G';
  rev['G']='C';



  barcode = barcodefinder(FileName, length);
  cout << "Barcode length: " << barcode << endl;



  map<string, int> kmers;

  int l = mink;
  int ll = l;
  int h = maxk;
  int i;
  int numkmers = 0;

  int logmake;
  int logcheck;
  int logvals[h];

  for (ll=l; ll<(h+1); ll++)
  {
    numkmers += pow(4,ll);
    cout << "Total kmer kombos: " << numkmers << endl;
  }

  for (logmake = 0; logmake<(h+1); logmake++)
  {
    logvals[logmake] = pow(4, logmake);
  }

  string ckmer;

  for (l=l; l<(h+1); l++)
  {
    ckmer = string(l, bases[0]);
    //cout << ckmer << endl;
    kmers[ckmer] = 0;
    //cout << " FIRST DONE " << endl;
    for (i=0; i<(pow(4,l)); i++)
    {
      ckmer.replace((l-1), 1, string(1, bases[i%4]));
      //cout << ckmer << " i-val: " << i << endl;
      logcheck = 1;
      if (floor(float(i)/4) != ceil(float(i)/4))
      {
        kmers[ckmer] = 0;
        //cout << " 123 DONE " << endl;
      }
      for (logcheck = 1; logcheck<l; logcheck++)
      {
        if (floor(float(i)/logvals[logcheck])==ceil(float(i)/logvals[logcheck])  && (float(i)/logvals[logcheck])!=0)
        {
          ckmer.replace(l-logcheck-1, 1, string(1, bases[(i/logvals[logcheck])%4]));
          //cout << ckmer << "   lc: " << logcheck << " lv: " << logvals[logcheck] << " i-val: " << i << " div: " << (float(i)/float(logvals[logcheck])) << " base: " << bases[(i/logvals[logcheck])%4] << " position: " << l-1-logcheck << ckmer[l-1-logcheck] << endl;
          if(floor(float(i)/4) == ceil(float(i)/4) && bases[(i/logvals[logcheck])%4] != 'A')
          {
            kmers[ckmer] = 0;
            //cout << " BIGGEST 4 DONE " << endl;
          }
        }
      }
    }
  }
  //kmers = kmergen(3,4);

  int lk = mink;
  int hk = maxk;

  string line;

  ifstream infile(FileName);
  if (infile.is_open())
  {
    while (getline(infile, line))
    {
      if (line[0] != '>' && line.length() == length)
      {
        ++numoflines;
        cout << "Line: " << line << endl;
        lk = mink;
        for (lk=lk; lk<(hk+1); lk++)
        {
          pos = 0;
          for (pos = barcode; pos < (length-barcode-lk+1); pos++)
          {
            cout << "Pos: " << pos << endl;
            cout << "Kmer: " << line.substr(pos, lk) << endl;
            kmers[line.substr(pos, lk)]++;
          }
        }
      }
    }
  }

  infile.close();

  cout << "Total lines: " << numoflines << endl;

  map<string, int>::iterator it;

  for (it=kmers.begin(); it!=kmers.end(); ++it)
    cout << it->first << " => " << it->second << '\n';

  cout << "Map size: " << kmers.size() << endl;


  return 0;
}
