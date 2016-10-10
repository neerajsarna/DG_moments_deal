namespace Ordering_Values
{
	using namespace std;
	struct Ordering
	{
  		Ordering( VectorXd &Array );
  		bool operator() (int i,int j);    
  		VectorXi index;  
  		VectorXd array;
	};

/******************************************************************/
/******************************************************************/
/******************************************************************/


Ordering::Ordering( VectorXd &Array )
{
  int n = Array.size();
  array = Array;
  
  vector<int> sorted(n);
  for( int i = 0; i<n; i++ ) sorted[i] = i;
  sort(sorted.begin(),sorted.end(),*this);
  index = Map<VectorXi>(sorted.data(),n);
}

bool Ordering::operator() (int i,int j) 
{ 
  return( array[i]<array[j] );
}
}
