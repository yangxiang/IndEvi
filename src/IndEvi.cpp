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
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include "PTNode.hh"
#include "TreeWriter.hh"
#include "CartesianProductDb.hh"
#include "IndEviStatic.hh"


struct tipair{// This is a pair of transaction and item
	Transactionset::transaction_t row;//actually row+1 in the array
	Itemset::item_t col;//actually col+1 in the array
};

struct bicliqueSize{
	int transaction_width;
	int item_length;
	int area; //area=transaction_width*item_length;
};

struct edgeinfo{
	bicliqueSize bestclique;
	int cliquecovercount;
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

struct edgecomp0{//by area
	bool operator() (const edgeinfo& lhs, const edgeinfo& rhs) const
	{
		
		if (lhs.bestclique.area>rhs.bestclique.area)
			return true;
		
		if ((lhs.bestclique.area==rhs.bestclique.area) && (lhs.cliquecovercount>rhs.cliquecovercount))
			return true;
		
		return false;
	}
};

bool ecomp0 (const edgeinfo& lhs, const edgeinfo& rhs) //by area
{
	if (lhs.bestclique.area>rhs.bestclique.area)
		return true;

	if ((lhs.bestclique.area==rhs.bestclique.area) && (lhs.cliquecovercount>rhs.cliquecovercount))
		return true;		
		
	return false;
}

bool bicomp0(const bicliqueSize& lhs, const bicliqueSize& rhs)//by area
{
		
		if (lhs.area>rhs.area)
			return true;
		
		return false;
}

struct edgecomp1{//by shortest edge, then by area
	bool operator() (const edgeinfo& lhs, const edgeinfo& rhs) const
	{
		int l_min=lhs.bestclique.transaction_width<lhs.bestclique.item_length?lhs.bestclique.transaction_width:lhs.bestclique.item_length;
		int r_min=rhs.bestclique.transaction_width<rhs.bestclique.item_length?rhs.bestclique.transaction_width:rhs.bestclique.item_length;
		
		if (l_min>r_min)
			return true;
		
		if ((l_min==r_min) && (lhs.bestclique.area>rhs.bestclique.area))
			return true;
			
		if ((l_min==r_min) && (lhs.bestclique.area==rhs.bestclique.area) && (lhs.cliquecovercount>rhs.cliquecovercount))
			return true;
			
		return false;
	}
};

bool ecomp1 (const edgeinfo& lhs, const edgeinfo& rhs)  //by shortest edge, then by area
{
	int l_min=lhs.bestclique.transaction_width<lhs.bestclique.item_length?lhs.bestclique.transaction_width:lhs.bestclique.item_length;
	int r_min=rhs.bestclique.transaction_width<rhs.bestclique.item_length?rhs.bestclique.transaction_width:rhs.bestclique.item_length;
		
	if (l_min>r_min)
		return true;
		
	if ((l_min==r_min) && (lhs.bestclique.area>rhs.bestclique.area))
		return true;

	if ((l_min==r_min) && (lhs.bestclique.area==rhs.bestclique.area) && (lhs.cliquecovercount>rhs.cliquecovercount))
		return true;
			
	return false;
}
	
bool bicomp1(const bicliqueSize& lhs, const bicliqueSize& rhs)//by shortest edge, then by area
{
		int l_min=lhs.transaction_width<lhs.item_length?lhs.transaction_width:lhs.item_length;
		int r_min=rhs.transaction_width<rhs.item_length?rhs.transaction_width:rhs.item_length;
		
		if (l_min>r_min)
			return true;
		
		if ((l_min==r_min) && (lhs.area>rhs.area))
			return true;
		return false;
}



int main(int argc, char **argv)
{

  if (argc < 5)
  {
    std::cerr << "Usage: " << argv[0] << " fi|fci|mfi dataset support [query_phenotype] \n";
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
  

  
  

  int query_phenotype=-1;
  if (argc >= 5)
  {
    query_phenotype = strtol(argv[4], NULL, 10);
  }
  

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

  
  
 //below are code specifically for IndEvi
 
  
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
  
  
  //starting transactional data completion
  std::vector<edgeinfo> score_matrix[IndEviStatic::num_transactions];
  
  bool trans_visit_record[IndEviStatic::num_transactions];
  for (uint32_t i=1; i<=IndEviStatic::num_transactions; i++)
  {  
     trans_visit_record[i-1]=false;
	 score_matrix[i-1].resize(IndEviStatic::max_item-IndEviStatic::min_item+1);
	 
  }
  for (uint32_t i=1; i<=IndEviStatic::num_transactions; i++)
	for (uint16_t j=IndEviStatic::min_item; j<=IndEviStatic::max_item; j++)
	{
		score_matrix[i-1][j-IndEviStatic::min_item].bestclique.transaction_width=-1;
		score_matrix[i-1][j-IndEviStatic::min_item].bestclique.item_length=-1;
		score_matrix[i-1][j-IndEviStatic::min_item].bestclique.area=-1;
		score_matrix[i-1][j-IndEviStatic::min_item].cliquecovercount=0;
	}
  
  std::map<tipair, int, tipaircomp> maxCliques;
  
  int bicliques_counter=0;
  for (CartesianProductDb::iterator ists=bicliques.begin(); ists!=bicliques.end(); ists++)
  {
	 
	 bicliques_counter++;
	 std::vector<int> row_projection;
	 std::vector<int> col_projection;
	 
	 row_projection.resize(IndEviStatic::num_transactions);
	 col_projection.resize(IndEviStatic::max_item-IndEviStatic::min_item+1);
	 
	 for (Transactionset::transaction_t ti=1; ti<=IndEviStatic::num_transactions; ti++)
			row_projection[ti-1]=0;
	 
	 for (Itemset::item_t isi=IndEviStatic::min_item; isi!=IndEviStatic::max_item; isi++)
			col_projection[isi-IndEviStatic::min_item]=0;
 
	 
	 for (Transactionset::transaction_t ti=1; ti<=IndEviStatic::num_transactions; ti++)
	 {
		for (Itemset::iterator isi=(ists->itemset).begin(); isi!=(ists->itemset).end(); isi++)
		{
			if (!uncov_zero_db.exists(ti, *isi))
				row_projection[ti-1]++;
		}
	 }

	 for (Itemset::item_t isi=IndEviStatic::min_item; isi<=IndEviStatic::max_item; isi++)
	 {
		for (Transactionset::iterator ti=(ists->transactionset).begin(); ti!=(ists->transactionset).end(); ti++)
		{
			if (!uncov_zero_db.exists(*ti, isi))
				col_projection[isi-IndEviStatic::min_item]++;
		}
	 }

	 


	 for (Transactionset::transaction_t ti=1; ti<=IndEviStatic::num_transactions; ti++)
	 {
		if (row_projection[ti-1]<=0)
			continue;	

			
		for (Itemset::item_t isi=IndEviStatic::min_item; isi<=IndEviStatic::max_item; isi++)
		{
			if (col_projection[isi-IndEviStatic::min_item]<=0)
				continue;
			int support_T=row_projection[ti-1];
			int support_I=col_projection[isi-IndEviStatic::min_item];
			
			
			if (!uncov_zero_db.exists(ti, isi))
			{
				Transactionset::transaction_t match_t[1];
				match_t[0]=ti;
				Transactionset::iterator it_proj_t=std::search((ists->transactionset).begin(),(ists->transactionset).end(), match_t, match_t+1);
				if (it_proj_t!=(ists->transactionset).end()) 
					support_T--;				
				if (support_T<=0)
					continue;

				Itemset::item_t match_i[1];
				match_i[0]=isi;	
				Itemset::iterator it_proj_i=std::search((ists->itemset).begin(),(ists->itemset).end(), match_i, match_i+1);
				if(it_proj_i!=(ists->itemset).end())
					support_I--;
					
				if (support_I<=0)
					continue;				
			}
			
			
			bicliqueSize score;
			score.transaction_width=support_T;
			score.item_length=support_I;
			score.area=score.transaction_width*score.item_length;
				
			score_matrix[ti-1][isi-IndEviStatic::min_item].cliquecovercount++;
			if (bicomp0(score, score_matrix[ti-1][isi-IndEviStatic::min_item].bestclique))
			{
				score_matrix[ti-1][isi-IndEviStatic::min_item].bestclique=score;
				tipair tmp_Clique;
				tmp_Clique.row=ti;
				tmp_Clique.col=isi;
				std::map<tipair, int, tipaircomp>::iterator Clit=maxCliques.find(tmp_Clique);
				if (Clit==maxCliques.end())
					maxCliques.insert(Clit, std::pair<tipair, int>(tmp_Clique, bicliques_counter));
				else
					Clit->second=bicliques_counter;					
			}				
		}
	 }
  }
  //ending transactional data completion
  
  
  
  //Output Cliques Begin
  std::ostringstream out_maxCliques;
  out_maxCliques << IndEviStatic::ds_fpath.substr(0, IndEviStatic::ds_fpath.length() - 4);
  out_maxCliques << "_support-" << support<<"_"<<fi_type;
  out_maxCliques << ".max_Cliques";
  
  std::cout << "Writing max cliques to " << out_maxCliques.str() << "\n";
  std::cout << "(Each line records the index of a maximal biclique of M)\n";
  
  std::ofstream maxCliques_outp(out_maxCliques.str().c_str());
  
  for (std::map<tipair, int, tipaircomp>::iterator Clit=maxCliques.begin(); Clit!=maxCliques.end(); Clit++)
	 maxCliques_outp<<"("<<(Clit->first).row<<","<<(Clit->first).col<<") : "<<Clit->second<<"\n";
	 
  maxCliques_outp.close();
  //Output Cliques End  
  
  //Output Score Matrix Begin
  std::ostringstream out_matrix_trwidth, out_matrix_itmlength, out_matrix_covercount;
  
  //Output trwidth
  out_matrix_trwidth << IndEviStatic::ds_fpath.substr(0, IndEviStatic::ds_fpath.length() - 4);
  out_matrix_trwidth << "_support-" << support<<"_"<<fi_type;
  out_matrix_trwidth << ".score_matrix_trwidth";

  std::cout << "Writing score matrix to " << out_matrix_trwidth.str() << "\n";
  std::cout << "(Each line corresponds to a transaction)\n";

  std::ofstream matrix_outp_trwidth(out_matrix_trwidth.str().c_str());
  
  for (uint32_t i=1; i<=IndEviStatic::num_transactions; i++)
  {  
	 for (uint16_t j=0; j<=IndEviStatic::max_item-IndEviStatic::min_item; j++)
		matrix_outp_trwidth << score_matrix[i-1][j].bestclique.transaction_width <<" ";
	 matrix_outp_trwidth<<"\n";
  }
  
  matrix_outp_trwidth.close();
  
  //Output itmlength
  out_matrix_itmlength << IndEviStatic::ds_fpath.substr(0, IndEviStatic::ds_fpath.length() - 4);
  out_matrix_itmlength << "_support-" << support<<"_"<<fi_type;
  out_matrix_itmlength << ".score_matrix_itmlength";

  std::cout << "Writing score matrix to " << out_matrix_itmlength.str() << "\n";
  std::cout << "(Each line corresponds to a transaction)\n";

  std::ofstream matrix_outp_itmlength(out_matrix_itmlength.str().c_str());
  
  for (uint32_t i=1; i<=IndEviStatic::num_transactions; i++)
  {  
	 for (uint16_t j=0; j<=IndEviStatic::max_item-IndEviStatic::min_item; j++)
		matrix_outp_itmlength << score_matrix[i-1][j].bestclique.item_length<<" ";
	 matrix_outp_itmlength<<"\n";
  }
  
  matrix_outp_itmlength.close();
  
  //Output covercount
  out_matrix_covercount << IndEviStatic::ds_fpath.substr(0, IndEviStatic::ds_fpath.length() - 4);
  out_matrix_covercount << "_support-" << support<<"_"<<fi_type;
  out_matrix_covercount << ".score_matrix_covercount";

  std::cout << "Writing score matrix to " << out_matrix_covercount.str() << "\n";
  std::cout << "(Each line corresponds to a transaction)\n";

  std::ofstream matrix_outp_covercount(out_matrix_covercount.str().c_str());
  
  for (uint32_t i=1; i<=IndEviStatic::num_transactions; i++)
  {  
	 for (uint16_t j=0; j<=IndEviStatic::max_item-IndEviStatic::min_item; j++)
		matrix_outp_covercount << score_matrix[i-1][j].cliquecovercount<<" ";
	 matrix_outp_covercount<<"\n";
  }
  
  matrix_outp_covercount.close();  
  
  //Output Score Matrix End
  
  //Begin outputting gene ranks and known-gene ranks
  //Begin calculating foldenrich
    std::multimap<float, tipair> foldenrich;
	
	 
	    std::ostringstream in_GeneList;
        in_GeneList << ds_dir;
        in_GeneList << ds_fbase+"_GeneList.txt";
	    std::ifstream GeneList_inp(in_GeneList.str().c_str());
		if (!GeneList_inp) {//
			std::cout << "Error: Cannot open " << in_GeneList.str() << std::endl;//
			return -1;//
		}//		
	    std::map<int, std::string> gList;
	    int id=1;
		std::string gene;
	    while (getline(GeneList_inp,gene))
	    {
		   gList.insert(std::pair<int, std::string>(id,gene));
		   id++;	
	    }	 
	    GeneList_inp.close();
		
		for (int p=IndEviStatic::min_item; p<=IndEviStatic::max_item; p++)
		{
		   std::multimap<edgeinfo, int, edgecomp0> gene_scores;
		   for (int i=0; i<=(int)IndEviStatic::num_transactions-1; i++)
		   {
			  gene_scores.insert(std::pair<edgeinfo, int>(score_matrix[i][p-IndEviStatic::min_item],i));//second is gene number -1;
		   }

		   std::vector<float> gene_pecentile;
		   gene_pecentile.resize(IndEviStatic::num_transactions);
		   int rank=1;
		   int samerank=0;
		   std::multimap<edgeinfo, int, edgecomp0>::iterator it=gene_scores.begin();
		   std::multimap<edgeinfo, int, edgecomp0>::iterator it_prev=it;
		   it++;
		   std::vector<int> samerankindex;//
		   for (; it!=gene_scores.end(); it++)
		   {
		       if ((!ecomp0(it_prev->first, it->first)) && (!ecomp0(it->first, it_prev->first)))
			   {
				  samerankindex.push_back(it->second);//
				  samerank++;
			   }
			   else
			   {
			      if(samerank>0)//
				  {//
					  for(int j=1; j<=samerank; j++)//
					  {//
						if(samerank-j<0){
							std::cout<<"error in samerank reassign."<<std::endl;
							std::cout<<"rank:"<<rank<<std::endl;
							std::cout<<"samerank:"<<samerank<<std::endl;
							std::cout<<"j:"<<j<<std::endl;
						}						
					  }//
				  }//
				  rank=rank+samerank+1;
				  for (int j=0; j<=samerank; j++)
				  {
					gene_pecentile[it_prev->second]=((float)(rank-1))/(IndEviStatic::num_transactions);
					if (!uncov_zero_db.exists(it_prev->second+1,p))
					{
						tipair tmp_ti;
						tmp_ti.col=p;
						tmp_ti.row=it_prev->second+1;
						foldenrich.insert(std::pair<float,tipair>(gene_pecentile[it_prev->second], tmp_ti));
					}					
					it_prev++;
				  }
 				  samerank=0;
				  samerankindex.clear();//
			   }
		   }
		   if(samerank>0)//
		   {//
			  for(int j=1; j<=samerank; j++)//
			  {//
				if(samerank-j<0){
					std::cout<<"error in samerank reassign."<<std::endl;
					std::cout<<"rank:"<<rank<<std::endl;
					std::cout<<"samerank:"<<samerank<<std::endl;
					std::cout<<"j:"<<j<<std::endl;
				}						
			  }//
		  }//		   
		   for (; it_prev!=it; it_prev++)
		   {
				gene_pecentile[it_prev->second]=((float)(rank+samerank))/(IndEviStatic::num_transactions);
				if (!uncov_zero_db.exists(it_prev->second+1,p))
				{
					tipair tmp_ti;
					tmp_ti.col=p;
					tmp_ti.row=it_prev->second+1;
					foldenrich.insert(std::pair<float,tipair>(gene_pecentile[it_prev->second], tmp_ti));	
				}				
		   }

			if (query_phenotype==p)
			{
			   std::ostringstream out_TopGenes;
			   out_TopGenes << IndEviStatic::ds_fpath.substr(0, IndEviStatic::ds_fpath.length() - 4);
			   out_TopGenes << "_support-" << support<<"_"<<fi_type;
			   out_TopGenes << ".TopGenes_" <<argv[4];
			  
			   std::cout << "Writing Gene ranks to " << out_TopGenes.str() << "\n";
			   std::cout << "(Each line corresponds to a gene)\n";
			  
			   std::ofstream TopGenes_outp(out_TopGenes.str().c_str());
				if (!TopGenes_outp) {//
					std::cout << "Error: Cannot open " << out_TopGenes.str() << std::endl;//
					return -1;//
				}//				   
			  
			   for (std::multimap<edgeinfo, int, edgecomp0>::iterator GSit=gene_scores.begin(); GSit!=gene_scores.end(); GSit++)
			   {
				  TopGenes_outp<<"Gene No: "<<GSit->second+1<<"; Gene Name: "<<(gList.find(GSit->second+1))->second;
				  TopGenes_outp<<"; transaction_width (genes): "<<(GSit->first).bestclique.transaction_width<<"; item_length (phenotypes): "<<(GSit->first).bestclique.item_length;
				  TopGenes_outp<<"; number of supporting bicliques: "<<(GSit->first).cliquecovercount<<"; pecentile: "<<gene_pecentile[GSit->second]<<std::endl;
			   }
			   TopGenes_outp<<std::endl<<"Below are ranks of known genes:"<<std::endl;
				
			   for (std::multimap<edgeinfo, int, edgecomp0>::iterator GSit=gene_scores.begin(); GSit!=gene_scores.end(); GSit++)
			   {
				  if (uncov_zero_db.exists(GSit->second+1,query_phenotype))
					continue;
				  
				  TopGenes_outp<<"Gene No: "<<GSit->second+1<<"; Gene Name: "<<(gList.find(GSit->second+1))->second;
				  TopGenes_outp<<"; transaction_width: "<<(GSit->first).bestclique.transaction_width<<"; item_length: "<<(GSit->first).bestclique.item_length;
				  TopGenes_outp<<"; number of supporting bicliques: "<<(GSit->first).cliquecovercount<<"; pecentile: "<<gene_pecentile[GSit->second]<<std::endl;
			   }		   
			   
			   TopGenes_outp.close();	
		    }
	
		}
	 
	 
  
  bicliques.clear();
  

  delete covered;

  delete gons;
 
  delete root;

  return 0;
}
