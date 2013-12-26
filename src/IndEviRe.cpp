#include <algorithm>
#include <vector>
#include <map>
#include <assert.h>
#include <errno.h>
#include <ext/algorithm>
#include <ext/functional>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <math.h>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "PTNode.hh"
#include "TreeWriter.hh"
#include "BitDb.hh"
#include "CartesianProductDb.hh"
#include "IndEviStatic.hh"


struct tipair{
	Transactionset::transaction_t row;//actually row+1 in the array
	Itemset::item_t col;//actually col+1 in the array
};


struct tipaircomp{
	bool operator() (const tipair& lhs, const tipair& rhs) const
	{
		if (lhs.row<rhs.row)
			return true;
		
		if ((lhs.row==rhs.row) && (lhs.col<rhs.col))
			return true;
		return false;
	}
};



int main(int argc, char **argv)
{
  
  if (argc < 6)
  {
    std::cerr << "Usage: " << argv[0] << " fi|fci|mfi dataset support trans_index item_index\n";
    std::cerr << "  fi => (normal) frequent itemsets, fci => frequent closed itemsets, mfi => maximal frequent itemsets\n";
    exit(1);
  }
  
  if (!(strcmp(argv[1], "fi") == 0) && !(strcmp(argv[1], "fci") == 0) && !(strcmp(argv[1], "mfi") == 0))
  {
    std::cerr << "Error: invalid frequent itemset type " << argv[1] << ".  Valid types are fi, fci, and mfi.\n";
    exit(1);
  }
  const char* fi_type = argv[1];

  struct timeval tv;
  gettimeofday(&tv, 0);
  srand48(tv.tv_sec + tv.tv_usec);

  IndEviStatic::init_timer();

  errno = 0;
  double support = strtod(argv[3], NULL);
  if (errno != 0)
  {
    std::cerr << "Error converting support " << argv[3] << " to double.\n";
    exit(1);
  }

  int trans_index = strtol(argv[4], NULL, 10);
  
  int item_index=strtol(argv[5], NULL, 10);


  IndEviStatic::ds_fpath = argv[2];
  std::string::size_type last_slash_pos = IndEviStatic::ds_fpath.rfind('/');
  std::string ds_dir(IndEviStatic::ds_fpath.substr(0, last_slash_pos + 1));
  std::string ds_fname(IndEviStatic::ds_fpath.substr(last_slash_pos + 1));
  std::string ds_fbase(ds_fname.substr(0, ds_fname.find('.')));

  std::string fi_dir(ds_dir + "freq_itemsets");
  if (mkdir(fi_dir.c_str(), 0775) != 0)
  {
    if (errno != EEXIST)
    {
      perror("mkdir freq_itemsets");
      exit(1);
    }
  }

  std::ostringstream fi_fname_str;
  fi_fname_str << fi_dir << "/" << ds_fbase << "_" << support << "_" << fi_type << ".mafia";
  const std::string fi_fname(fi_fname_str.str());

  // regenerate frequent itemset file if necessary
  struct stat ds_stat, fi_stat;
  if (stat(IndEviStatic::ds_fpath.c_str(), &ds_stat) == -1)
  {
    perror("stat dataset file");
    exit(1);
  }
  if (stat(fi_fname.c_str(), &fi_stat) != -1)
  {
    if (fi_stat.st_mtime <= ds_stat.st_mtime)
    {
      std::cout << "frequent itemset file " << fi_fname << " is not newer than dataset file " << IndEviStatic::ds_fpath << "; deleting frequent itemset file\n";
      if (unlink(fi_fname.c_str()) == -1)
      {
        perror("unlink fi_fname");
        exit(1);
      }
    }
  }

  FILE *fi_file;
  if ((fi_file = fopen(fi_fname.c_str(), "rb")) == NULL)
  {
    if (errno != ENOENT)
    {
      perror("fopen fi_fname (try 1)");
      exit(1);
    }
    std::cout << "freq itemset file " << fi_fname << " does not exist; creating.\n";

    std::ostringstream mk_fi_cmd;
    //mk_fi_cmd << "bin/mafia -fi " << support / 100.0 << " -ascii " << IndEviStatic::ds_fpath << " " << fi_fname;
    mk_fi_cmd << "bin/mafia -" << fi_type << " " << support / 100.0 << " -ascii " << IndEviStatic::ds_fpath << " " << fi_fname;
    std::cout << "Executing " << mk_fi_cmd.str() << "\n";
    int fi_cmd_st;
    if ((fi_cmd_st = system(mk_fi_cmd.str().c_str())) == -1)
    {
      std::cerr << "call to `" << mk_fi_cmd.str() << "` failed\n";
      perror("make fi file");
      exit(1);
    }

    if (WEXITSTATUS(fi_cmd_st) != 0)
    {
       std::cerr << "command `" << mk_fi_cmd.str() << "` returned " << WEXITSTATUS(fi_cmd_st) << "\n";
       exit(1);
    }

    if ((fi_file = fopen(fi_fname.c_str(), "rb")) == NULL)
    {
      perror("fopen fi_fname (try 2)");
      exit(1);
    }
  }

  std::cout << "using freq itemset file " << fi_fname << "\n";

  PTNode* root = new PTNode();

  size_t line_size = 0;
  char *fi_line = NULL;

  while (getline(&fi_line, &line_size, fi_file) != -1)
  {
    //std::cout << "Reading line " << fi_line;
    int item_digit_len, item_int;
    std::set<Itemset::item_t> items;
    char *fi_line_scan = fi_line;
    while (sscanf(fi_line_scan, "%d%n", &item_int, &item_digit_len) == 1)
    {
      //std::cout << "  item " << item_int << "\n";
      if (item_int < std::numeric_limits<Itemset::item_t>::min() || item_int > std::numeric_limits<Itemset::item_t>::max())
        {
        std::cout << "Item out of range: " << item_int << "\n";
          exit(1);
        }
        const Itemset::item_t item = static_cast<Itemset::item_t>(item_int);
      items.insert(item);
      fi_line_scan += item_digit_len;
    }
    long abs_support;
    // TODO: adapt to mafia output format
    if (sscanf(fi_line_scan, " (%ld)", &abs_support) != 1)
    {
      std::cerr << "Error parsing absolute support from '" << fi_line_scan << "'\n";
      std::cerr << "Expected a string like \" (abs_support)\"\n";
      exit(1);
    }
    //std::cout << "  absolute support: " << abs_support << "\n";

    if (items.size() > 0)
    {
  
      PTNode *ptNode = root->build_node(items);

      ptNode->coverable_itemset = true;
      
  
      ptNode->gamma = (float)(items.size() + abs_support) / (float)(items.size() * abs_support);
      //std::cout << "adding node " << ptNode << " with calculated gamma " << ptNode->gamma << " from items.size() " << items.size() << ", abs_support " << abs_support << "\n";
    }
  }

  if (NULL != fi_line)
  {
    free(fi_line);
  }

  if (fclose(fi_file) == -1)
  {
    perror("fclose fi_file");
    exit(1);
  }


  FILE *ds_file;
  if ((ds_file = fopen(IndEviStatic::ds_fpath.c_str(), "rb")) == NULL)
  {
    std::cerr << IndEviStatic::ds_fpath << "\n";
    perror("fopen ds_fpath");
    exit(1);
  }

  char *ds_line;
  if ((ds_line = (char*)malloc(line_size)) == NULL)
  {
    perror("malloc ds_line");
    exit(1);
  }

  // populate first level of nodes with transactions by scanning the dataset
  IndEviStatic::max_item = std::numeric_limits<Itemset::item_t>::min();
  IndEviStatic::min_item = std::numeric_limits<Itemset::item_t>::max();

  uint64_t trans_num_int64 = 0, tot_cells_num = 0;
  int item_int, item_digit_len;
  char *ds_line_scan;
  PTNode::PTChildMap::iterator item_child;
  std::set<Itemset::item_t> added_single_items;//for removal single items added as non-MBC;    
  while (getline(&ds_line, &line_size, ds_file) != -1)
  {
    ++trans_num_int64;
    if (trans_num_int64 < std::numeric_limits<Transactionset::transaction_t>::min() || trans_num_int64 > std::numeric_limits<Transactionset::transaction_t>::max())
    {
      std::cout << "Transaction id out of range: " << trans_num_int64 << "\n";
      exit(1);
    }
    ds_line_scan = ds_line;
    while (sscanf(ds_line_scan, "%d%n", &item_int, &item_digit_len) != EOF)
    {
      if (item_int < std::numeric_limits<Itemset::item_t>::min() || item_int > std::numeric_limits<Itemset::item_t>::max())
      {
        std::cout << "Error: item out of range: " << item_int << "\n";
        exit(1);
      }
      const Itemset::item_t item = static_cast<Itemset::item_t>(item_int);
      //std::cout << "Processing transaction " << trans_num_int64 << ", item " << item << ", item_digit_len " << item_digit_len << "\n";
      IndEviStatic::max_item = std::max(IndEviStatic::max_item, item);
      IndEviStatic::min_item = std::min(IndEviStatic::min_item, item);

      PTNode *child;
      if ((item_child = root->children.find(item)) == root->children.end())
      {
        std::set<Itemset::item_t> single_item;
        single_item.insert(item);
        child = root->build_node(single_item);
		added_single_items.insert(item);
        child->coverable_itemset = true;
        //std::cout << "built single node child with item " << child->itemset << " (size " << child->itemset.size() << ")\n";
      } else {
        child = item_child->second;
        //std::cout << "retreived single node " << child << "\n";
      }
      //std::cout << "  child " << item << " matches\n";
      child->transactionset.push_back(static_cast<Transactionset::transaction_t>(trans_num_int64));
      ds_line_scan += item_digit_len;
      ++tot_cells_num;
    }
  }

  IndEviStatic::num_transactions = static_cast<Transactionset::transaction_t>(trans_num_int64);
  
  // populate with a bit for every zero
  BitDb uncov_zero_db(IndEviStatic::num_transactions, IndEviStatic::min_item, IndEviStatic::max_item);

  if (fseek(ds_file, 0, SEEK_SET) != 0)
  {
    perror("fseek ds_file");
    exit(1);
  }

  Transactionset::transaction_t cur_trans = 0;
  while (getline(&ds_line, &line_size, ds_file) != -1)
  {
    ++cur_trans;
    //std::cout << "populating uncov_zero_db for line " << ds_line;
    std::vector<bool> line_items(IndEviStatic::max_item);
    ds_line_scan = ds_line;
    while (sscanf(ds_line_scan, "%d%n", &item_int, &item_digit_len) != EOF)
    {
      const Itemset::item_t item = static_cast<Itemset::item_t>(item_int);
      line_items[item] = true;
      ds_line_scan += item_digit_len;
    }
    for (Itemset::item_t i=IndEviStatic::min_item; i<=IndEviStatic::max_item; ++i)
    {
      if (!line_items[i])
      {
        //std::cout << "  inserting uncov bit for trans " << cur_trans << ", item " << i << "\n";
        uncov_zero_db.insert(cur_trans, i);
      }
    }
  }

  free(ds_line);

  if (fclose(ds_file) == -1)
  {
    perror("fclose ds_fpath");
    exit(1);
  }

  //std::cout << "uncov_zero_db:\n" <<  uncov_zero_db << "\n";

  for (PTNode::PTChildMap::iterator rcIt=root->children.begin(); rcIt!=root->children.end(); ++rcIt)
  {
    PTNode *child = rcIt->second;
    //std::cout << "checking root child with itemset: " << child->itemset << "\n";
    if (child->gamma == std::numeric_limits<float>::max())
    {
      child->gamma = (float)(child->transactionset.size() + child->itemset.size()) / (float)((double)child->transactionset.size() * (double)child->itemset.size());
      //std::cout << "updated root child's init_support gamma to " << child->gamma << "\n";
    }
  }

  std::cout << "orig db size_horiz: " << tot_cells_num + IndEviStatic::num_transactions << "\n";
  std::cout << "orig db size_min: " << IndEviStatic::num_transactions + (IndEviStatic::max_item - IndEviStatic::min_item + 1) << "\n";


  BitDb *covered = new BitDb(IndEviStatic::num_transactions, IndEviStatic::min_item, IndEviStatic::max_item);

  typedef PTNode::GONodeSet GONodeSet;

  GONodeSet* gons = root->gamma_ordered_nodeset();

  CartesianProductDb bicliques;

    
 //below are code specifically for IndEviRe

  
  
  for(GONodeSet::iterator ists=gons->begin(); ists!=gons->end(); ists++)
  {
	Transactionset* trans_chosen = new Transactionset;
	(*ists)->recalc_gamma(*covered, trans_chosen, root);
	if(((*ists)->itemset).size()==1){
		 Itemset::item_t tmp_single_item=*(((*ists)->itemset).begin());
		 if (added_single_items.find(tmp_single_item)!=added_single_items.end())
			continue;
	}
	
	CartesianProduct cp((*ists)->itemset, *trans_chosen);
    bicliques.push_back(cp);
	delete trans_chosen;
    trans_chosen = NULL;
  }
  
  std::ostringstream in_GeneList;
  in_GeneList << ds_dir;
  in_GeneList << ds_fbase+"_GeneList.txt";

  std::ostringstream in_PhenotypeList;
  in_PhenotypeList << ds_dir;
  in_PhenotypeList << ds_fbase+"_PhenotypeList.txt";
 
  std::ifstream GeneList_inp(in_GeneList.str().c_str());
  std::ifstream PhenotypeList_inp(in_PhenotypeList.str().c_str());
  if ((!GeneList_inp)||(!PhenotypeList_inp)) {//
		std::cout << "Error: Cannot open " << in_GeneList.str() <<" or "<< in_PhenotypeList.str() << std::endl;//
		return -1;//
	}//    
  
  std::map<int, std::string> gList;
  std::map<int, std::string> phList;
  
int id=1;
  std::string gene;
  while (getline(GeneList_inp, gene))
  {
	gList.insert(std::pair<int, std::string>(id,gene));
	id++;	
  }
  std::cout<< "Gene list size: "<<gList.size()<<std::endl;
  GeneList_inp.close();
  
  id=IndEviStatic::min_item;
  std::string phenotype;
  while (getline(PhenotypeList_inp,phenotype))
  {
	phList.insert(std::pair<int, std::string>(id,phenotype));
	id++;	
  }
  std::cout<< "phenotype list size: "<<phList.size()<<std::endl;
  PhenotypeList_inp.close();
  
  std::ostringstream in_maxCliques;
  in_maxCliques << IndEviStatic::ds_fpath.substr(0, IndEviStatic::ds_fpath.length() - 4);
  in_maxCliques << "_support-" << support<<"_"<<fi_type;
  in_maxCliques << ".max_Cliques";
  
  std::cout << "Reading max cliques from " << in_maxCliques.str() << "\n";
  
  std::ifstream maxCliques_inp(in_maxCliques.str().c_str()); 
  if (!maxCliques_inp) {//
		std::cout << "Error: Cannot open " << in_maxCliques.str() << std::endl;//
		return -1;//
	}//    
  
  int tr_num=-1;
  int itm_num=-1;
  int bicliques_num=-1;
  
  
  std::string max_clique;
  bool clique_found=false;
  while(getline(maxCliques_inp, max_clique))
  {
	size_t left_para=max_clique.find("(");
	size_t mid_comma=max_clique.find(",");
	size_t right_para=max_clique.find(")");
	size_t sep_colon=max_clique.find(":");
	tr_num=atoi((max_clique.substr(left_para+1, mid_comma-left_para-1)).c_str());
	itm_num=atoi((max_clique.substr(mid_comma+1, right_para-mid_comma-1)).c_str());
	bicliques_num=atoi((max_clique.substr(sep_colon+2)).c_str());
	if ((trans_index==tr_num)&&(item_index==itm_num)){
		clique_found=true;
		break;
	}
  }
  
  if (!clique_found)
  {
	 std::cout << "No supporting clique detected. Exit."<<"\n";
	 exit(1);
  }
  else
  {
	std::cout<< "reading successful: tr_num "<<tr_num<<" itm_num "<<itm_num<< " bicliques_num "<<bicliques_num<<"\n";
  }  
  maxCliques_inp.close();
  if (!uncov_zero_db.exists(tr_num, itm_num))
  {
	std::cout<<"("<<tr_num<<","<<itm_num<<")=1"<<std::endl;
	std::cout<<"("<<(gList.find(tr_num))->second<<","<<(phList.find(itm_num))->second<<")=1"<<std::endl;
  }
  else
  {
	std::cout<<"("<<tr_num<<","<<itm_num<<")=0"<<std::endl;
	std::cout<<"("<<(gList.find(tr_num))->second<<","<<(phList.find(itm_num))->second<<")=0"<<std::endl;
  }
	
  int bicliques_counter=0;
  
  //starting reconstruct evidence
  for (CartesianProductDb::iterator ists=bicliques.begin(); ists!=bicliques.end(); ists++)
  {
	bicliques_counter++;
	if (bicliques_counter<bicliques_num)
		continue;
	
	std::cout<<"The maximal biclique index: "<<bicliques_counter<<std::endl<<"The maximal biclique: "<<*ists<<std::endl;
	std::cout<<"The phenotype-gene relationship corresponding to the maximal biclique:"<<std::endl;

	for(Itemset::iterator itit=ists->itemset.begin(); itit!=ists->itemset.end(); itit++)
		std::cout<<(phList.find(*itit))->second<<", ";			
		
	std::cout<<"; ";
		
	for(Transactionset::iterator trit=ists->transactionset.begin(); trit!=ists->transactionset.end(); trit++)
		std::cout<<(gList.find(*trit))->second<<", ";	

	std::cout<<std::endl;	
	
	
	Transactionset::transaction_t match_t[1];
	match_t[0]=tr_num;
	Itemset::item_t match_i[1];
	match_i[0]=itm_num;
	
	if (!uncov_zero_db.exists(tr_num, itm_num))
	{
		Transactionset::iterator Setri=std::search((ists->transactionset).begin(), (ists->transactionset).end(), match_t, match_t+1);
		Itemset::iterator Seiti=std::search((ists->itemset).begin(), (ists->itemset).end(), match_i, match_i+1);
		if ((Setri!=ists->transactionset.end()) && (Seiti!=ists->itemset.end()))
		{
			std::cout<<"The biclique covers the entry."<<std::endl;
			break;
		}			
	}


	
	std::cout<<"The biclique does NOT cover the entry. The maximal supporting biclique is: "<<std::endl;

	int effect_itm=0;	
	for(Itemset::iterator itit=ists->itemset.begin(); itit!=ists->itemset.end(); itit++)
	{
		if (!uncov_zero_db.exists(tr_num, *itit))
		{
			std::cout<<(phList.find(*itit))->second<<", ";
			effect_itm++;
		}
	}
	
	std::cout<<"; ";
	
	int effect_tr=0;	
	
	for(Transactionset::iterator trit=ists->transactionset.begin(); trit!=ists->transactionset.end(); trit++)
	{
		if (!uncov_zero_db.exists(*trit, itm_num))
		{
			std::cout<<(gList.find(*trit))->second<<", ";
			effect_tr++;
		}
	}	
	
	std::cout<<"("<<effect_itm<<"*"<<effect_tr<<")"<<std::endl;
	
	break;

  }
  //ending reconstruct evidence
  
   bicliques.clear();
  

  delete covered;

  delete gons;
 
  delete root;

  return 0;
} 
