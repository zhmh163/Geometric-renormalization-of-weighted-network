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
double rand_num;
// Mersenne Twister generate random number 
int SEED = std::time(NULL);//Will be used to obtain a seed for the random number engine
std::mt19937 engine(SEED);    
std::uniform_real_distribution<double> Rnd(0.0, 1.0); //random number in [0, 1)

//*********************************************************************//
void Read_parameters(string filename, int &N, double &beta, double &mu, double &a, double &eta, double &alpha,double &epsilon2)
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
				iss >> N>> beta>> mu>> a >> eta>> alpha >>epsilon2>>chi2;				
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
//*********************************************************************//
vector<vector<int> > Read_Weight_Graph(string filename, vector<vector<double> > &weight, int &N) //should check the file name.
{
	//filename+=".txt";
	vector<vector<int> > adj_list;	
	vector<vector<double> > weight_list;	
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
				while(iss >> i >> j >> wij )
				{	
					while(i>=adj_list.size())
					{			
						adj_list.push_back(vector<int>());	
						weight_list.push_back(vector<double>());
						
					}
					while(j>=adj_list.size())
					{			
						adj_list.push_back(vector<int>());
						weight_list.push_back(vector<double>());
						
					}
					if ( (i!=j) && (find(adj_list[i].begin(), adj_list[i].end(), j) == adj_list[i].end()) )// can't self and multi-connection!
					{
						adj_list[i].push_back(j);	
						adj_list[j].push_back(i);	
						weight_list[i].push_back(wij);	
						weight_list[j].push_back(wij);							
					}
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
	weight=weight_list;	
	N=adj_list.size();
	//cout<<N<<endl;
	return adj_list;
}
//*********************************************************************//
void Read_coordinates(string filename, vector<double> &kappa,vector<double> &theta)	
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
void Cal_sigma(vector<double> kappa, vector<double> &sigma, double a, double eta)
{	
	sigma.clear();	
	for(int i=0;i<kappa.size();i++)
	{
		sigma.push_back(a*pow(kappa[i],eta));		
	}
	
}
double Cal_nu(vector<double> sigma, double alpha, double beta, double mu, double I2)
{
	double ave_sigma = std::accumulate( sigma.begin(), sigma.end(), 0.0) / (double)(sigma.size()); 
	double I3=pi/(beta*sin((pi-pi*alpha)/beta));
	double nu=1.0/(2*pow(mu,(1-alpha))*I2*I3*ave_sigma);
	return nu;
}
//*********************************************************************//
vector<vector<double> > Cal_Yl(vector<vector<int> > edge_list, vector<vector<double> > weight, \
    vector<vector<int> > edge_super, vector<vector<int> > supernode_has_nodes, int Nl1,int r)
{
	vector<vector<double> > Yl(r*r+1);
	for(int supernode1=0;supernode1<Nl1;supernode1++)
	{
		for(int l=0;l<edge_super[supernode1].size();l++)
		{
			int supernode2=edge_super[supernode1][l];
			if(supernode2>supernode1)
			{
				vector<double> wk;
				wk.clear();
				for(int s1=0;s1<supernode_has_nodes[supernode1].size();s1++)
				{
					int m=supernode_has_nodes[supernode1][s1];
					for(int s2=0;s2<supernode_has_nodes[supernode2].size();s2++)
					{
						int n=supernode_has_nodes[supernode2][s2];
						for(int ix=0;ix<edge_list[m].size();ix++)
						{
							if(edge_list[m][ix]==n) // link exists 
							{
								wk.push_back(weight[m][ix]);
								break;
							}
						}
					}
				}
				double sum_wk=std::accumulate( wk.begin(), wk.end(), 0.0) ;
				double temp_yl=0;
				for(int k=0;k<wk.size();k++)
				{
					temp_yl+=pow(wk[k]/sum_wk,2.0);					
				}
				int k=(int)(wk.size());
				Yl[k].push_back(temp_yl);

			}

		}

	}
    return Yl;
}
//*********************************************************************//
void Renormalization(vector<vector<int> > &edge_list, vector<vector<double> > &weight, vector<int> &node_in_which_supernode,\
	vector<double> &kappa, vector<double> &theta, vector<double> &sigma, double alpha, double mu, double nu, double beta, \
	double eta, double a, int r,double D,double W0,vector<vector<double> > &Yl)
{	
	
	int Nl1; //Number of nodes in super layer;
	int Nn=edge_list.size();//Number of nodes in sub-node layer;
	if(Nn%r==0)
		Nl1=Nn/r;
	else
		Nl1=Nn/r+1;
	// initialization and define vector in super layer
	vector<vector<int> > edge_super(Nl1);
	vector<vector<double> > weight_super(Nl1);
	//vector<vector<double> > noise_super(Nl1);
	
	vector<double> kappa_super(Nl1, 0.0);
	vector<double> theta_super(Nl1, 0.0);
	vector<double> sigma_super(Nl1, 0.0);
	
	vector<int> temp_node_in_which_supernode(Nn,-1);
	node_in_which_supernode=temp_node_in_which_supernode;
	
	vector<vector<int> > supernode_has_nodes(Nl1);
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
				supernode_has_nodes[i].push_back(j);
				sumkappa+=pow(kappa[j],beta);				
				sumtheta+=pow((theta[j]*kappa[j]),beta);
			}
		}
		kappa_super[i]=pow(sumkappa,1.0/beta);
		theta_super[i]=pow(sumtheta/sumkappa,1.0/beta);
		sigma_super[i]=a*pow(kappa_super[i],eta);
	}
	// renormalization the edge list
	double Phi=beta/(D*(eta-1)+alpha);
	//double C=pow(r,alpha/D);
	double C=1.0;
	
	
	int supernode1, supernode2;
	for(int m=0;m<edge_list.size();m++)
	{
		supernode1=node_in_which_supernode[m];
		for(int l=0;l<edge_list[m].size();l++)
		{
			int n=edge_list[m][l];
			if(n>m)
			{
				supernode2=node_in_which_supernode[n];
				int flag=0;
				for(int k=0;k<edge_super[supernode1].size();k++)
				{
					if(edge_super[supernode1][k]==supernode2)
					{
						flag=1;
						weight_super[supernode1][k]+=pow(weight[m][l],Phi);
						for(int k2=0;k2<edge_super[supernode2].size();k2++)
						{
							if(edge_super[supernode2][k2]==supernode1)
							{
								weight_super[supernode2][k2]+=pow(weight[m][l],Phi);
								break;
							}
						}
						break;
					}
				}			
				if (flag==0 && supernode1!=supernode2 )
				{
					edge_super[supernode1].push_back(supernode2);
					edge_super[supernode2].push_back(supernode1);
					double wij=pow(weight[m][l],Phi);
					weight_super[supernode1].push_back(wij);	
					weight_super[supernode2].push_back(wij);	

					//try no noise 
					//noise_super[supernode1].push_back(1.0);	
					//noise_super[supernode2].push_back(1.0);	
				}
			}		
		}	
		
	}
	 
	//normalizate C so that the average weight=W0 
	vector<double> w_list;
	for(int i=0;i<weight_super.size();i++)
	{
		for(int l=0;l<weight_super[i].size();l++)
		{
			weight_super[i][l]=pow(weight_super[i][l],1.0/Phi);
			w_list.push_back(weight_super[i][l]);
		}
	}		
	double W=accumulate( w_list.begin(), w_list.end(), 0.0) / w_list.size();
	C=W0/W;
	for(int i=0;i<weight_super.size();i++)
	{
		for(int l=0;l<weight_super[i].size();l++)
		{
			weight_super[i][l]=C*weight_super[i][l];			 
		}
	}	
	
	// Cal Yl in each group
	Yl=Cal_Yl(edge_list, weight, edge_super, supernode_has_nodes, Nl1, r);
	//return everything including node_in_which_supernode in super layer
	edge_list=edge_super;
	weight=weight_super;
	//noise=noise_super;
	
	kappa=kappa_super;
	theta=theta_super;
	sigma=sigma_super;	
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
vector<double > Cal_Yi(vector<vector<int> > edge_list, vector<vector<double> > weight_list)
{
	int Nn=edge_list.size();
	vector<double> Yi;	 
	for (int i=0; i<Nn; i++)
	{		 
		double si=0;
		for(int j=0;j<weight_list[i].size();j++)
		{	
			double w=weight_list[i][j];
			si+=w;
		}
		double sum2=0;
		for(int j=0;j<weight_list[i].size();j++)
		{	
			double w=weight_list[i][j];
			sum2+=pow(w/si,2.0);
		}
		Yi.push_back(sum2);		 
	}
	return Yi;
	 
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
void output_edgelist(string filename, vector<vector<int> > edge_list, vector<vector<double> > weight, \
	vector<double> kappa, vector<double> theta, vector<double> sigma, vector<double> strength, vector<double> Yi)
{
	
	string file1,file2;
	file1=filename+"_weight_edgelist.txt";	
	file2=filename+"_coordinates.txt";	
	ofstream outfile10(file1.c_str(),ios::out); 	
	ofstream outfile11(file2.c_str(),ios::out); 
	
	outfile10<<setprecision(10)<<left<<" "<<setw(10)<<"# i"<<" "<<setw(10)<<"j"<<" "<<setw(16)<<"wij"<<endl;
	outfile11<<setprecision(10)<<left<<" "<<setw(10)<<"# i"<<" "<<setw(16)<<"kappa"<<" "<<setw(16)<<"theta"<<" "<<setw(16)<<"sigma"<<" "<<setw(16)<<"strength"<<" "<<setw(16)<<"degree"<<" "<<setw(16)<<"disparity"<<endl;
	int N=edge_list.size();
	for(int i=0;i<N;i++)
	{
		for(int k=0;k<edge_list[i].size();k++)
		{
			if(edge_list[i][k]>i)
				outfile10<<setprecision(10)<<" "<<left<<setw(10)<<i<<" "<<setw(10)<<edge_list[i][k]<<" "<<setw(16)<<weight[i][k]<<endl;		
		}
		if(edge_list[i].size()>0)
		outfile11<<setprecision(10)<<left<<" "<<setw(10)<<i<<" "<<setw(16)<<kappa[i]<<" "<<setw(16)<<theta[i]<<" "<<setw(16)<<sigma[i]<<" "<<setw(16)<<strength[i]<<" "<<setw(16)<<(int)(edge_list[i].size())<<" "<<setw(16)<<Yi[i]<<endl;
	}
	outfile10.close();
	outfile11.close();
}

void output_parameters(string filename, int Nn, double beta, double mu,double nu, double alpha, double eta,double a)
{
	string file1;
	file1=filename+"_parameters.txt";		
	ofstream outfile10(file1.c_str(),ios::out); 
	outfile10<<setprecision(10)<<left<<" "<<setw(20)<<"# N"<<" "<<setw(20)<<"beta"<<" "<<setw(20)<<"mu"<<" "<<setw(20)<<"nu";
	outfile10<<" "<<setw(20)<<"alpha"<<" "<<setw(20)<<"eta"<<" "<<setw(20)<<"a"<<endl;	
	
	outfile10<<setprecision(10)<<left<<" "<<setw(20)<<Nn<<" "<<setw(20)<<beta<<" "<<setw(20)<<mu<<" "<<setw(20)<<nu;
	outfile10<<" "<<setw(20)<<alpha<<" "<<setw(20)<<eta<<" "<<setw(20)<<a<<endl;		
	outfile10.close();
	
}
void output_supernode(string filename, vector<int>  node_in_which_supernode)
{
	string file1;
	file1=filename+"_supernode.txt";		 
	ofstream outfile10(file1.c_str(),ios::out); 	
	 
	for(int i=0;i<node_in_which_supernode.size();i++)
	{		
		outfile10<<left<<setw(15)<<i<<" "<<node_in_which_supernode[i]<<endl;
	}
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
void Cal_mean_weight_strength(vector<vector<double> > weight_list, double &ave_w, double &ave_s, vector<double> &strength)
{
	//cal strength 
	
	vector<double> s;
	vector<double> w_list;
	for(int i=0;i<weight_list.size();i++)
	{
		double si=0;
		for(int j=0;j<weight_list[i].size();j++)
		{
			si+=weight_list[i][j];
			w_list.push_back(weight_list[i][j]);
		}
		s.push_back(si);	
	}
	ave_w=accumulate( w_list.begin(), w_list.end(), 0.0) / w_list.size();
	ave_s=accumulate( s.begin(), s.end(), 0.0) / s.size();	
	strength=s;
}
void Cal_ave_std(vector<double> X, double &ave, double &std)
{
	ave=0;
	std=0;
	ave=std::accumulate( X.begin(), X.end(), 0.0) / (double)(X.size()); 
	double sum=0;
	for(int i=0;i<X.size();i++)
	{
		sum+=(X[i]-ave)*(X[i]-ave);
	}
	sum/= (double)(X.size());
	std = sqrt(sum);	
}
void output_group_disparity(string filename, vector<vector<double> > Yl)
{
	
	string file1 ;
	file1=filename+"_Yl.txt";	 
	ofstream outfile10(file1.c_str(),ios::out); 	
		
	outfile10<<setprecision(10)<<left<<" "<<setw(10)<<"# l"<<" "<<setw(16)<<"ave.Yl"<<" "<<setw(16)<<"std.Yl"<<" "<<setw(16)<<"size.Yl"<<endl;
	 
	for(int l=1;l<Yl.size();l++)
	{
		double ave_Yl;
		double std_Yl;
		Cal_ave_std(Yl[l], ave_Yl, std_Yl);
		outfile10<<setprecision(10)<<left<<" "<<setw(10)<<l<<" "<<setw(16)<<ave_Yl<<" "<<setw(16)<<std_Yl<<" "<<setw(16)<<Yl[l].size()<<endl;
	}
	outfile10.close();
	 
}
void renormalization_weight_net(string file_edge, string file_parameters, string file_coordinate, string file_outname,int Layers)
{
	
		//topological parameters
		int N;	
		int realN;
		double beta;	
		double mu;
		//double gamma;
		//weight parameters	
		double alpha;
		double nu;	
		double eta;	
		double a;
		double epsilon2;

		// vector for edge and weight 
		vector<vector<int> > edge_list;	
		vector<vector<double> > weight;
		 
		vector<double> theta;
		vector<double> kappa;
		vector<double> sigma;
		vector<double> strength;
		
		// RG parameters
		int r=2;
		double D=1;
		vector<int> node_in_which_supernode;
		double meanC=0.0;
		double meanK=0.0;
		double ave_w=0;
		double ave_s=0;
				
		//change the file name dynamical
		Read_parameters(file_parameters, N, beta, mu, a, eta, alpha,epsilon2); 	
		edge_list=Read_Weight_Graph(file_edge, weight, N);
		Read_coordinates(file_coordinate, kappa, theta); 	
		Cal_sigma(kappa, sigma, a, eta);
		nu=Cal_nu(sigma, alpha, beta, mu, 1.0);//I2=average noise, here is 1.0
		//start RG for weight network
		string file_renormalization= file_outname+"_layer_C_K_N.txt";
		ofstream outfile10(file_renormalization,ios::out);

		vector<double > Cc;
		vector<double> Yi;	
		Cc=cluster_coefficence(edge_list,N,meanC);
		meanK=meank_maxk(edge_list,N);
		realN=Real_N(edge_list,N);
		Cal_mean_weight_strength(weight, ave_w, ave_s, strength);		 
		Yi=Cal_Yi(edge_list, weight);
		double W0=ave_w;
	
		int layer=0;
		cout<<left<<left<<" "<<setw(10)<<layer<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<ave_w<<" "<<setw(20)<<ave_s<<" "<<setw(20)<<mu<<endl;
		
		outfile10<<left<<" "<<setw(10)<<"# layer"<<" "<<setw(20)<<"meanC"<<" "<<setw(20)<<"meanK"<<" "<<setw(20)<<"realN"<<" "<<setw(20)<<"N"<<" "<<setw(20)<<"ave_weight"<<" "<<setw(20)<<"ave.strength"<<" "<<setw(20)<<"mu"<<endl;
		outfile10<<left<<" "<<setw(10)<<layer<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<ave_w<<" "<<setw(20)<<ave_s<<" "<<setw(20)<<mu<<endl;

		ostringstream file5;
		file5 <<file_outname<<"_RG_layer_"<<layer;
		string file_layer_0= file5.str();		
		output_edgelist(file_layer_0, edge_list, weight, kappa, theta, sigma, strength, Yi);	
		output_parameters(file_layer_0, N, beta, mu, nu, alpha, eta, a);
 
		//renormalization
		for(int layer=1;layer<Layers;layer++)
		{		
 			vector<vector<double> > Yl; 
			Renormalization(edge_list, weight, node_in_which_supernode, kappa, theta, sigma, alpha, mu, nu, beta, eta, a, r, D, W0,Yl);
			 
			//cal strength for super nodes
			Cal_mean_weight_strength(weight, ave_w, ave_s, strength);//return average weight, ave strength and strength for each nodes. 
			// Cal disparity in each node i and group l
			Yi=Cal_Yi(edge_list, weight);
			
			//new network size is rescaled to N/r;
			N=edge_list.size();
			mu=mu/r;
			//nu=nu*pow(r,-alpha/D);
			//update nu each layer.
			nu=Cal_nu(sigma, alpha, beta, mu, 1.0);//I2=average noise, here is 1.0
			
			//cal mean C and mean degree			
			Cc=cluster_coefficence(edge_list, N, meanC);
			meanK=meank_maxk(edge_list,N);
			realN=Real_N(edge_list,N);
            			
			
			cout<<left<<" "<<setw(10)<<layer<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<ave_w<<" "<<setw(20)<<ave_s<<" "<<setw(20)<<mu<<endl;
			outfile10<<left<<" "<<setw(10)<<layer<<" "<<setw(20)<<meanC<<" "<<setw(20)<<meanK<<" "<<setw(20)<<realN<<" "<<setw(20)<<N<<" "<<setw(20)<<ave_w<<" "<<setw(20)<<ave_s<<" "<<setw(20)<<mu<<endl;

			//change the file name dynamical 
			ostringstream fileNameStream;
			fileNameStream <<file_outname<<"_RG_layer_"<<layer;
			string fileName = fileNameStream.str();	
		
			output_edgelist(fileName, edge_list, weight, kappa, theta, sigma, strength,Yi);	
			output_parameters(fileName, N, beta, mu, nu, alpha, eta, a);			
			output_supernode(fileName, node_in_which_supernode);
			output_group_disparity(fileName, Yl);
		}
		outfile10.close();
	 
	
}
/*
int main()
{
 	 
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
	for(int n=0;n<12; n++)
	{	 
		renormalization_weight_net(edgelist_filename[n], parameters_filename[n], coordinate_filename[n], outfile_rootname[n], Layers[n]);
	}
 
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
 	renormalization_weight_net(edgelist_filename, parameters_filename, coordinate_filename, outfile_rootname, Layers);
	return 0;
}