#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <stdio.h>
#include <algorithm>    // std::find
#include <vector>       // std::vector
#include <set>			//std::set
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <sstream>
#include <numeric>
#include <random>
using namespace std;

#define pi          3.141592653589793238462643383279502884197   
// Mersenne Twister generate random number 
int SEED = std::time(NULL);//Will be used to obtain a seed for the random number engine
std::mt19937 engine(SEED);    
std::uniform_real_distribution<double> Rnd(0.0, 1.0); //random number in [0, 1)

double rand_num;



vector<vector<int> > edge_list;
vector<double> theta;
vector<double> kappa;
void Read_parameters(string filename, int &realN, double &beta, double & mu) //should check the file name.
{
	// # N	 beta	mu	a	eta	 alpha_opt	epsilon2	chi2 
	//filename+=".txt";	
	string line;
	ifstream file(filename.c_str());
	if (file.is_open())
	{		
		while (getline(file, line))
		{
			if (line[0] != '#' ) //without comment
			{
				std::istringstream iss(line);
				//int xN;
				double chi2;
				iss >> realN>> beta>> mu;	
				iss.clear();
			}
		}
	}
	else 
	{
		cout << "Unable to open file: "<<filename<<endl; 
		exit (EXIT_FAILURE);
	}
	file.close();	

}
vector<vector<int> > Read_Graph(string filename, int &N) //should check the file name.
{
		//filename+=".txt";
	vector<vector<int> > adj_list;		 
	string line;
	ifstream file(filename.c_str());
	if (file.is_open())
	{		
		while (getline(file, line))
		{
			if (line[0] != '#' ) //without comment
			{
				std::istringstream iss(line);
				int i,j;
				double wij;	// 		
				while(iss >> i >> j )
				{	
					while(i>=adj_list.size())
					{			
						adj_list.push_back(vector<int>());							
						
					}
					while(j>=adj_list.size())
					{			
						adj_list.push_back(vector<int>());						
						
					}
					if ( (i!=j) && (find(adj_list[i].begin(), adj_list[i].end(), j) == adj_list[i].end()) )// can't self and multi-connection!
					{
						adj_list[i].push_back(j);	
						adj_list[j].push_back(i);												
					}
					//skip the last
					iss.clear();
				}
			}
		}
	}
	else 
	{
		cout << "Unable to open file: "<<filename<<endl; 
		exit (EXIT_FAILURE);
	}
	file.close();	
	N=adj_list.size();
	//cout<<N<<endl;
	return adj_list;
}
void Read_kappa_theta(string filename, vector<double> &kappa, vector<double> &theta) //should check the file name.
{
	//# i        kappa            theta  
	//filename+=".txt";	
	kappa.clear();
	theta.clear();
	string line;
	ifstream file(filename.c_str());
	if (file.is_open())
	{		
		while (getline(file, line))
		{
			if (line[0] != '#' ) //without comment
			{
				std::istringstream iss(line);
				int i;
				double xkappa,xtheta;	// 		
				while(iss >> i >> xkappa >> xtheta)
				{	
					while(i>=kappa.size())
					{			
						kappa.push_back(0);	
						theta.push_back(0);												
					}					
					kappa[i]=xkappa;	
					theta[i]=xtheta;								
				}
			}
		}
	}
	else 
	{
		cout << "Unable to open file: "<<filename<<endl; 
		exit (EXIT_FAILURE);
	}
	file.close();	
}
void sort_index(vector<vector<int> > &edge_list, vector<double> &kappa, vector<double> &theta, int &N) //re-asign the id for each node
{
    // Declaring vector of pairs
    vector< pair <double, int> > vect;
 
    // Initializing 1st and 2nd element of pairs of theta
    // Entering values in vector of pairs
    for (int i=0; i<N; i++)
	{
		if(edge_list[i].size()>0)
		{
			vect.push_back( make_pair(theta[i], i) );
		}
	}  

    // Using simple sort() function to sort
    sort(vect.begin(), vect.end());
 
     // Printing the sorted vector(after using sort())
    vector<double> new_theta;
	vector<double> new_kappa;
	vector<int> new_id(N,0);
    for (int i=0; i<vect.size(); i++)
    {
        // vect[i].first-->theta, vect[i].second-->old index
        // 1st and 2nd element of pair respectively
		new_theta.push_back(vect[i].first);
		new_kappa.push_back( kappa[vect[i].second] );
        new_id[vect[i].second]=i;
    } 
	new_theta.swap(theta);
	new_kappa.swap(kappa);
	vector<vector<int> > new_edge_list(vect.size());
	for(int i=0;i<N;i++)
	{
		for(int j=0;j<edge_list[i].size();j++)
		{
			int a=new_id[i];
			int b=new_id[edge_list[i][j]];
			new_edge_list[a].push_back(b);
		}
	
	}
	new_edge_list.swap(edge_list);
	N=vect.size();
	//cout<<"N is"<<N<<endl;
   // return new_id;
}
vector<vector<int> > Renormalization(vector<vector<int> > edge_list, int Nn, vector<double> &kappa, vector<double> &theta, vector<int> &node_in_which_supernode, int r,double beta)
{	
	
	int Nl1; //Number of nodes in layer 1;
	if(Nn%r==0)
		Nl1=Nn/r;
	else
		Nl1=Nn/r+1;
	vector<double> theta_l1(Nl1, 0.0);
	vector<double> kappa_l1(Nl1, 0.0);
	//node i in wich supernode 
		
	vector<vector<int> > aij(Nl1);  //aij edge list in layer 1
	//calculate kappa, theta and supernode set
 	for(int i=0;i<Nl1;i++)
	{
		double sumkappa=0.0;
		double sumtheta=0.0;
		for(int j=i*r; j<(i+1)*r; j++)
		{
			if(j<Nn)
			{				
				node_in_which_supernode[j]=i;//node j in supernode i.
				sumkappa+=pow(kappa[j],beta);				
				sumtheta+=pow((theta[j]*kappa[j]),beta);
			}
		}
		kappa_l1[i]=pow(sumkappa,1.0/beta);
		theta_l1[i]=pow(sumtheta/sumkappa,1.0/beta);
	}
	// renormalization the edge list
	int supernode1, supernode2;
	for(int i=0;i<Nn;i++)
	{
		supernode1=node_in_which_supernode[i];
		for(int l=0;l<edge_list[i].size();l++)
		{
			int j=edge_list[i][l];
			if(j>i)
			{
				supernode2=node_in_which_supernode[j];
				if ( (supernode1!=supernode2) && (find(aij[supernode1].begin(), aij[supernode1].end(), supernode2) == aij[supernode1].end()) )// can't self and multi-connection!
				{
					aij[supernode1].push_back(supernode2);
					aij[supernode2].push_back(supernode1);
				}
			}
		
		}
	
	
	}
	
	//reture kappa, theta and node_in_which_supernode in layer 1;
	kappa.swap(kappa_l1);
	theta.swap(theta_l1);
	//return aij
	//cout<<"finish Renormalization"<<endl;
	return aij;
}
vector<double > cluster_coefficence(vector<vector<int> > edge_list, int Nn, double &meanC) //The clustering coefficient of each node
{
	int i,j,k,u,v;
	vector<double > Cluster_coeff;
	double C,ci;	
	C=0.0;
	double sum=0; //sum is the number of degree larger than 1
	for(i=0;i<Nn;i++)
	{
		ci=0.0;
		if(edge_list[i].size()>1)//degree larger than 1
		{
			for(j=0;j<edge_list[i].size();j++)
			{
				for(k=j+1;k<edge_list[i].size();k++)
				{
					u=edge_list[i][j];
					v=edge_list[i][k];
					if(find(edge_list[u].begin(),edge_list[u].end(),v) !=edge_list[u].end())
					{
						ci+=1.0;
					}
				}
			}		
			ci=2.0*ci/(edge_list[i].size()* (edge_list[i].size()-1.0));
			sum++;
		}
		else
			ci=0.0;
		Cluster_coeff.push_back(ci);
		C+=ci;
	}
	//C/=(double)(Nn);
	C/=sum; //sum is the number of degree larger than 1
	meanC=C;
	//cout<<"the mean clustering coefficient is:"<<'\t'<<meanC<<endl;
	return Cluster_coeff;
}
double meank_maxk(vector<vector<int> > edge_list, int Nn)
{
	double meank=0;
	int maxk=0;
	double sumn = 0;
	for(int i=0;i<Nn;i++)
	{
		if(edge_list[i].size()>maxk)
			maxk=edge_list[i].size();
		if (edge_list[i].size() > 0)
			sumn++;
		meank+=edge_list[i].size();
	}
	meank /= sumn;
	//cout<<"maxk="<<maxk<<'\t'<<"mean_degree="<<meank<<endl;
	return meank;
}
void output_edgelist(string filename, vector<vector<int> > edge_list,int Nn)
{
	string file1,file2;
	file1=filename+"_edgelist.txt";	
	file2=filename+"_coordinates.txt";	
	ofstream outfile10(file1.c_str(),ios::out); 	
	ofstream outfile11(file2.c_str(),ios::out); 
	outfile10<<setprecision(10)<<left<<" "<<setw(10)<<"# i"<<" "<<setw(10)<<"j"<<endl;
	outfile11<<setprecision(10)<<left<<" "<<setw(10)<<"# i"<<" "<<setw(16)<<"kappa"<<" "<<setw(16)<<"theta"<<endl;
	for(int i=0;i<Nn;i++)
	{
		for(int k=0;k<edge_list[i].size();k++)
		{
			if(edge_list[i][k]>i)				 
				outfile10<<setprecision(10)<<left<<" "<<setw(10)<<i<<" "<<setw(10)<<edge_list[i][k]<<endl;
		}
		if(edge_list[i].size()>0)		
		outfile11<<setprecision(10)<<left<<" "<<setw(10)<<i<<" "<<setw(16)<<kappa[i]<<" "<<setw(16)<<theta[i]<<endl;
	}
	outfile10.close();
	outfile11.close();
}
void output_supernode(string filename, vector<int>  node_in_which_supernode,int Nn)
{
	string file1;
	file1=filename+"_supernode.txt";		 
	ofstream outfile10(file1.c_str(),ios::out); 	
	 
	for(int i=0;i<Nn;i++)
	{		
		outfile10<<left<<setw(15)<<i<<" "<<node_in_which_supernode[i]<<endl;
	}
	outfile10.close();	
}
void output_parameters(string filename, int Nn, double beta, double mu)
{
	string file1;
	file1=filename+"_parameters.txt";		
	ofstream outfile10(file1.c_str(),ios::out); 
	outfile10<<setprecision(10)<<left<<" "<<setw(20)<<"# N"<<" "<<setw(20)<<"beta"<<" "<<setw(20)<<"mu"<<endl;		
	outfile10<<setprecision(10)<<left<<" "<<setw(20)<<Nn<<" "<<setw(20)<<beta<<" "<<setw(20)<<mu<<endl;		
	outfile10.close();
	
}
int Real_N(vector<vector<int> > edge_list, int Nn)
{
	
	int realN=0;	
	for(int i=0;i<Nn;i++)
	{		
		if (edge_list[i].size() > 0)
		{
			realN++;
		}
		
	}
	return realN;
}
void renormalization_unweighted_net(string file_edge, string file_parameters, string file_coordinate, string file_outname,int Layers)
{
	int realN;	
	int r=2;
	double beta;
	double mu;
	int N;         //Node number 
	double meanC=0.0;
	double meanK=0.0;
		//change the file name dynamical 
		
		Read_parameters(file_parameters, realN, beta, mu);		
		edge_list=Read_Graph(file_edge,N);		
		Read_kappa_theta(file_coordinate, kappa, theta);
		 

		string file_renormalization= file_outname+"_layer_C_K_N.txt";	
		ofstream outfile10(file_renormalization.c_str(),ios::out);
		
		vector<double > Cc;
		Cc=cluster_coefficence(edge_list,N,meanC);
		meanK=meank_maxk(edge_list,N);
		realN=Real_N(edge_list,N);
		cout<<left<<" "<<setw(10)<<0<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<mu<<endl;
		
		outfile10<<left<<" "<<setw(10)<<"# layer"<<" "<<setw(20)<<"meanC"<<" "<<setw(20)<<"meanK"<<" "<<setw(20)<<"realN"<<" "<<setw(20)<<"N"<<" "<<setw(20)<<"mu"<<endl;
		outfile10<<left<<" "<<setw(10)<<0<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<mu<<endl;
		
		//outfile10<<left<<setw(10)<<0<<" "<<setw(15)<<meanC<<" "<<setw(15)<<meanK<<" "<<setw(15)<<realN<<" "<<setw(15)<<N<<" "<<setw(15)<<mu<<endl;
		int layer=0;
		ostringstream file5;
		file5 <<file_outname<<"_RG_layer_"<<layer;
		string file_layer_0= file5.str();		
		output_edgelist(file_layer_0, edge_list, N);
		output_parameters(file_layer_0, N, beta, mu);
		
					
		//renormalization
		for(int layer=1;layer<Layers;layer++)
		{

			vector<int> node_in_which_supernode(N, -1);
			edge_list=Renormalization(edge_list, N, kappa, theta, node_in_which_supernode, r,beta);//return edge_list kappa theta and node_in_which_supernode.
			//new network size is rescaled to N/r;
			int Noriginal=N;
			if(N%r==0)
				N=N/r;
			else
				N=N/r+1;
			mu=mu/r;
			//cal mean C and mean degree 
			//cout<<"N "<<N<<endl;
			Cc=cluster_coefficence(edge_list, N, meanC);
			meanK=meank_maxk(edge_list,N);
			realN=Real_N(edge_list,N);

			cout<<left<<" "<<setw(10)<<layer<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<mu<<endl;
		
			outfile10<<left<<" "<<setw(10)<<layer<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<mu<<endl;
		
			//change the file name dynamical 
			 
			ostringstream fileNameStream;
			fileNameStream <<file_outname<<"_RG_layer_"<<layer;
			string fileName = fileNameStream.str();	
			
			 	
			output_edgelist(fileName, edge_list, N);
			output_parameters(fileName, N, beta, mu);	
			output_supernode(fileName, node_in_which_supernode,Noriginal);
			
		}
		outfile10.close();
	 
}
/*
int main()
{
	double time_begin=time(NULL);
	
	string edgelist_filename[12]={"data/cargoshipsBB_edgelist.txt", \
								"data/Ecoli_edgelist.txt", \								
								"data/UScommuteBB_edgelist.txt",\
								"data/TW_103_weighted_edgelist.txt",\								
								"data/TW_80_weighted_edgelist.txt",\
								"data/TW_65_weighted_edgelist.txt",\
								"data/Facebook_like_Social_Network_edgelist.txt",\
								"data/Newman_scientific_collaboration_network_sum_edgelist.txt",\
								"data/Openflights_edgelist.txt",\								
								"data/NZ_edgelist.txt",\
								"data/poppy_1_edgelist.txt",\
								"data/foxglove_2_edgelist.txt" 
								};
	
	string coordinate_filename[12]={"data/cargoshipsBB_coordinates.txt", \
								"data/Ecoli_coordinates.txt", \								
								"data/UScommuteBB_coordinates.txt",\
								"data/TW_103_coordinates.txt",\								
								"data/TW_80_coordinates.txt",\
								"data/TW_65_coordinates.txt",\
								"data/Facebook_like_Social_Network_coordinates.txt",\								
								"data/Newman_scientific_collaboration_network_sum_coordinates.txt",\
								"data/Openflights_coordinates.txt",\								
								"data/NZ_coordinates.txt",\
								"data/poppy_1_coordinates.txt",\
								"data/foxglove_2_coordinates.txt" 
								};
	
    string parameters_filename[12]={"data/cargoshipsBB_Final_parameters.txt", \
								"data/Ecoli_Final_parameters.txt", \								
								"data/UScommuteBB_Final_parameters.txt",\
								"data/TW_103_Final_parameters.txt",\								
								"data/TW_80_Final_parameters.txt",\
								"data/TW_65_Final_parameters.txt",\
								"data/Facebook_like_Social_Network_Final_parameters.txt",\								
								"data/Newman_scientific_collaboration_network_sum_Final_parameters.txt",\
								"data/Openflights_Final_parameters.txt",\								
								"data/NZ_Final_parameters.txt",\
								"data/poppy_1_Final_parameters.txt",\
								"data/foxglove_2_Final_parameters.txt" 
								};
	
	string outfile_rootname[12]={"data/cargoshipsBB", \
								"data/Ecoli", \								
								"data/UScommuteBB",\
								"data/TW_103",\								
								"data/TW_80",\
								"data/TW_65",\
								"data/Facebook_like_Social_Network",\								
								"data/Newman_scientific_collaboration_network_sum",\
								"data/Openflights",\								
								"data/NZ",\
								"data/poppy_1",\
								"data/foxglove_2"
								};
	 
	int Layers[12]={3, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
	for(int n=0;n<1; n++)
	{	 
		renormalization_unweighted_net(edgelist_filename[n], parameters_filename[n], coordinate_filename[n], outfile_rootname[n], Layers[n]);
	}
	double time_end=time(NULL);
    cout<<"RUN TIME: "<<time_end-time_begin<<endl;	
	return 0;
}
*/
int main(int argc, char *argv[])
{	
	
	//./GR data/cargoshipsBB_edgelist.txt data/cargoshipsBB_Final_parameters.txt data/cargoshipsBB_coordinates.txt data/cargoshipsBB 3
	string edgelist_filename=argv[1];	 
	string parameters_filename=argv[2];
	string coordinate_filename=argv[3];
	string outfile_rootname=argv[4];
	int Layers=stoi(argv[5]); //total_layer	
 	renormalization_unweighted_net(edgelist_filename, parameters_filename, coordinate_filename, outfile_rootname, Layers);
	return 0;
}
