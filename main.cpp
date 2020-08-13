// This main diagonalizes 1 boson hopping in 2 sites, with the 1D Bose-Hubbard model.

//#include <algorithm>
//#include <cmath>
//#include <iostream>
//#include <istream>
//#include <iterator>
//#include <map>
//#include <ostream>
//#include <string>
//#include <utility>
//#include <vector>
//#include <cstdlib>
//#include <Eigen/Core>
//#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <istream>
#include <iterator>
#include <map>
#include <ostream>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <cstdlib>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>




/**
	Basis_Representation is a class, which should do the following:
		prompt && store user defined parameters about the second quantized Hamiltonian, and
			use make_Vect_n to define a Basis(list) of all allowed states, and
				define a hash list[or Lin Table] from basis vector to basis vector, and
					output.

			Assumptions: LatDim==1[implicit] statistics s==+1 [implicit,bosons] spin
			*/

class Basis_Representation
{
	typedef std::vector<std::tuple<int,int>> Coord_List;
	typedef std::map<long long int, int>     Initial_States;
	typedef Eigen::Triplet<double>           T;
	typedef std::vector<int>                 Occupation_Number;
	typedef std::tuple<int, int>             coordinate;

	public:

		/**
			- CONSTRUCTORS -
			Basis_Representation	empty constructor
			Basis_Representation	accepts total number of particle && total lattice sites
			Basis_Representation	accepts total number of particle && total lattice sites && coupling parameters for Hamiltonian(t && J)

			- OUTPUTS -
			print_Basis				terminal print vector of tuple basis_vector_n
			print_All				terminal print all essential info about the invoked object.
		 */


			Basis_Representation();
			Basis_Representation(int, int);
			Basis_Representation(int, int, int, int);
			Basis_Representation(int, int, int, int, int);
			Basis_Representation(int, int, int, int, int, int, int);

			void set_Num_States();
			int  get_Num_States();
			int  get_Ini_States_size();
			int  getT();
			int  getJ();


			int  choose(int,int);

			void append(long long int, long long int);
			void append(long long int, std::string, std::map<long long int, std::string>);

			long long int
			     retrive(long long int);

			long long int  		to_Decimal(Occupation_Number);
			Occupation_Number   to_ONum(long long int);

			void print();
			void print(int);
			void print(unsigned int);
			void print(std::string);
			void print(Basis_Representation::Occupation_Number);
			void print(Basis_Representation::Occupation_Number, unsigned int);
			void print(long long int, bool);
			void print(long long int);
			void print_Ini_States();
			void print_Item(long long int);
			void print_Map(std::map<long long int, std::string>);
			void print_TList(std::vector<Eigen::Triplet<double>>);
			void TListtoFStates(std::vector<Eigen::Triplet<double>>);

			void 				iterate();

			long long int 		iterator();
			void 				reset();

			void 				make_Initial_States();

			long long int 	 	double_translate(long long int);

			Basis_Representation::Occupation_Number
								double_translate (Basis_Representation::Occupation_Number);

		/**
		    - HASH -
		 */
			std::map<long long int, std::string> create_at (int, std::map<long long int, std::string>);
			std::map<long long int, std::string> annihilate_at (int);
			std::map<long long int, std::string> annihilate_at (int, std::map<long long int, std::string>);

			int get_Particle_at(int, long long int);

		/**
		    - HAMILTONIAN MATRIX -
		 */
			void make_1PT_Operator();
			void make_2PT_Operator();
			void make_Chem_Pot_Ops();
			std::map<long long int, std::string> make_Another_Ops();

			std::map<long long int, std::string> hop(int, int);
			std::map<long long int, std::string> number_Operator(int);

			std::map<long long int, std::string> number_Operator
					(int, std::map<long long int, std::string>);

			std::map<long long int, std::string> add_constant
				 (int, std::map<long long int, std::string>);

			std::vector<std::string> split(std::string, char);

			void set_Final_States (std::map<long long int, std::string>);
			std::map<long long int, std::string> get_Final_States ();

			void set_TList(int);
			void add_TList(std::vector< Eigen::Triplet<double> >);

			Eigen::SparseMatrix<double> make_Hamiltonian ();

			std::vector< Eigen::Triplet<double> > to_Triplet_List(int);

			/**
			- COMPOSITE INDICES -
			 */

			int make_Coordinate(int, int);

			coordinate i_to_coord(int);
			int coord_to_i(coordinate);


	private:

		/**
			- USER INPUT VARIABLES -
			Particle_Num			Total number of particles
			Lattice_Len				Total number of lattice sites
			t						coupling parameter for 1pt term
			J						coupling parameter for 2pt term

			- INTERNAL VARIABLES -
			basis_vect_n			a list of tuples (0/1 0/1)
			lin_value				a list of decimal numbers encoding elements in basis_vect_n
			lin_key					a list of decimals encoding the keys for lin_value
			num_states				the dimension of the allowed occupation number basis

			- CONSTANTS -
			Lattice_Len				number of lattice sites

			- FUNCTIONS -
			hop 					translate particle  from site to site
			interact 				couples 2 particles between sites
			create		 			apply creation operator on the 		#th element on vector_n
			annihilate				apply annihilation operator on the  #th element on vector_n
			count_Num				apply number operator

			- MAKERS -
			make_Basis_Vect_n		write down the occupation number set allowed by the Hamiltonian
			make_Vect_n				appends tuple<int int> to basis_vect_n

			- SETTERS -
			set_Param				initialize all user input variables
			set_Int_Var				initialize all internal variables a±=s zeros

			- GETTERS -
			get_Part_Num			returns Particle_Num
			get_Lattice_Len			returns Lattice_Len
			get_t					returns coupling t
			get_J					returns coupling J

			get_Num_States			returns the total number of states

			get_Vect_n				retrieves an element of basis_vect_n i.e vector of tuples
			get_Vect_n_Elem			retrieves an element of tuple<int int> vector_n

			- HELPERS -
			put_Zero				initializes vector_n[ i ] to zero
			put_One					adds the integer value 1 in vector_n [ i ]
			Look_Up					uses either Lin Table or Hashing to retrieve elements in vector_n

			print_String			basically print a sentence in the terminal.

			- LOOK UP OPERATIONS -
			set_lin_key				converting basis_vector_n into lin_key
			set_lin_value			converting basis_vector_n into lin_value

		 */

			int Particle_Num,
				Lattice_Len,

				t,
				J,
				mu,
				rowLen,
				colLen;

			int     num_states;
			int	    state_size;
			int states_created;


			// map ( i , initial_states )
			std::map<long long int, int>
											initial_states;
			std::map<long long int,
					        std::string>
											final_states;

			std::vector< Eigen::Triplet<double>	>
											TList;

			Coord_List CList;

			void set_Param(int, int, int, int, int);
			void set_Param(int, int, int, int, int, int, int);
			void set_Int_Var();
			void set_State_Size();
};

Basis_Representation::Basis_Representation()
{
	set_Param(0,0,0,0,0);
}

Basis_Representation::Basis_Representation(int A, int B)
/**
	 @param
 */
{
	set_Param(A,B,0,0,0);
	set_Int_Var();

	this->make_Initial_States();
}

Basis_Representation::Basis_Representation(int A, int B, int C, int D)
{
	set_Param(A,B,C,D,0);
	set_Int_Var();

	this->make_Initial_States();

}

Basis_Representation::Basis_Representation(int A, int B, int C, int D, int E)
{
	this->set_Param(A,B,C,D,E);
	set_Int_Var();

	this->make_Initial_States();
}

Basis_Representation::Basis_Representation(int A, int B, int C, int D, int E, int F, int G)
{
	this->set_Param(A,B,C,D,E,F,G);
	set_Int_Var();

	this->make_Initial_States();
}

void Basis_Representation::set_Param(int A, int B, int C, int D, int E)
/**
   @ (Particle_Num, Lattice_Len, t, J)
 */
{
		Particle_Num = A;
		Lattice_Len	 = B;

				  t	 = C;
				  J  = D;
				  mu = E;
}

void Basis_Representation::set_Param(int A, int B, int C, int D, int E, int F, int G)
{
	Particle_Num = A;
	Lattice_Len	 = B;

			  t	 = C;
			  J  = D;
			  mu = E;
	      rowLen = F;
		  colLen = G;
}

void Basis_Representation::set_Int_Var()
/**
   @param initialize basis_vect_n && num_states
 	 uses void make_Basis_Vect_n
 */
{
	this->set_Num_States();
	this->set_State_Size();

	states_created = 0; // iterator
}

int  Basis_Representation::choose(int n, int k)
{
    int res = 1;

    // Since C(n, k) = C(n, n-k)
    if ( k > n - k )
        k = n - k;

    // Calculate value of
    // [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1]
    for (int i = 0; i < k; ++i)
    {
        res *= (n - i);
        res /= (i + 1);
    }

    return res;
}

void Basis_Representation::set_Num_States()
{
	int n = Particle_Num + Lattice_Len - 1 ,
		k = Particle_Num;

		num_states = choose(n,k);
}

int  Basis_Representation::get_Num_States()
{
	return num_states;
}

void Basis_Representation::set_State_Size()
{
	state_size = Lattice_Len;
}

void Basis_Representation::print()
{
	std::cout << '*';
}

void Basis_Representation::append(long long int key, long long int value)
{
	this->initial_states.insert
		(std::pair<long long int, long long int> (key, value));
}

long long int Basis_Representation::retrive(long long int key)
{
	return initial_states[key];
}

Basis_Representation::Occupation_Number Basis_Representation::to_ONum(long long int input)
{
	int L = this->Lattice_Len;
	Occupation_Number tempONum(L,0);

	int quotient = 0,
	   remainder = 0;
	lldiv_t     temp;

	long long int dividend =         input,
		           divisor =  (long long int)
				   	   	   	   (1 + this->Particle_Num);

	bool IsEnd   =  false;

	int    digit = 0;

	while(IsEnd  == false)
	{
			 temp = lldiv (dividend, divisor);

		 quotient = temp.quot;
		remainder =  temp.rem;

		tempONum.at(digit) = remainder;
		digit++;

//		this->print(tempONum.at(digit));

		dividend  = quotient;
		quotient  = 0;
		remainder = 0;

		if (dividend == 0 || digit > L)
			IsEnd = true;
	}

	return tempONum;
}

long long int Basis_Representation::to_Decimal(Basis_Representation::Occupation_Number input)
{
	int base   = 1 + this -> Particle_Num;
	long long int output = 0;

	for (int power = 0; power < this->Lattice_Len; power ++)
		output = output + input[power] * std::pow(base, power);

	return output;
}

void Basis_Representation::print(Basis_Representation::Occupation_Number input)
{
	int Len = this->Lattice_Len;

	std::cout << "(";

	for (int index = 0; index < Len; index++)
	{
		std::cout << input.at(index);

		if (index != Len-1)
			std::cout << " ";
		if (index == Len-1)
			std::cout << ") ";
	}
}

void Basis_Representation::print(int input)
{
	std::cout << input << std::endl;
}

void Basis_Representation::print(unsigned int input)
{
	std::cout << input << std::endl;
}

void Basis_Representation::print(std::string input)
{
//	std::cout << "*" << input << std::endl;
}

void Basis_Representation::print(Basis_Representation::Occupation_Number input, unsigned int count)
/**
	@param  accept ONum && Pos Int counter
 	@return prints ONum if counter is 0.
 */
{
	if (count == 0)
	{
		this->iterate();
		int Len = this->Lattice_Len;

		for (int index = 0; index < Len; index++)
		{
			std::cout  << input[index] << ' ';

			if (index == Len-1)
			{
				std::cout << std::endl;
			}
		}
	}

	else if (count > 0)
	{
//		std::cout << "Err Calling Function print(ONum, unassigned int)." << std::endl;
	}
}

void Basis_Representation::print(long long int key, bool temp)
/**
	@return print: key -> value(initial_states)
 */
{
	long long int temp_item = this->retrive(key);
	this->print(this->to_ONum(temp_item));
}

void Basis_Representation::print(long long int input)
{
	std::cout << input << std::endl;
}

void Basis_Representation::print_Item(long long int key)
{
	std::cout << key << ' ';
	this->print( this->to_ONum( this->retrive(key) ) );
}

void Basis_Representation::iterate()
/**
	@param iterate iterator states_created
	from class Basis_Representation
 */
{
	states_created++;
}

void Basis_Representation::reset()
/**
   assign iterator the value 0
*/
{
	this->states_created = 0;
}

long long int Basis_Representation::iterator()
{
	return states_created;
}

void Basis_Representation::make_Initial_States()
/**
	uses append(key value)
	assigns a key:(decimal)
	&& a ONum value: (int) interator w/ lexical ordering
 */
{
//	this->print(this->states_created);

	this->reset();

//	std::map<long long int, int>::iterator temp_ite
//		= initial_states.begin();

	int N = this->Particle_Num,
	    L = this->Lattice_Len;

	int j = 0;

	bool NotEnd = true;
	bool NotLast = true;

	std::vector<int> This(L,0);
	std::vector<int> Next(L,0);

	This.at(0) = N;

	this->iterate();
//		initial_states.insert(temp_ite,
//			std::pair<long long int, int>
//			(this->to_Decimal(This), -1+ this->states_created)
//						 );
		initial_states.emplace(this->to_Decimal(This), -1+ this->states_created);

//	std::cout << initial_states.at(this->to_Decimal(This));
//			  << std::endl;

	while (NotLast)
	// begin while loop: Basis_Generation
	{
		while (NotEnd)
		{

			std::copy(This.begin(), This.begin()+j, Next.begin());
			Next.at(j)   = This.at(j)   -1;

			Next.at(j+1) = Next.at(j+1) +1;

			std::copy( Next.begin(), Next.end(), This.begin() );
			std::fill( Next.begin(), Next.end(), 0 );

			this->iterate();

//			initial_states.insert(temp_ite,
//					std::pair<long long int, int>
//					(this->to_Decimal(This), -1+ this->states_created)
//						);
			initial_states.emplace(this->to_Decimal(This), -1+ this->states_created);

			j++;

			if ( j == L-1 )
				NotEnd = false;
		}

		if (j == L-1)
		/**
			If_loop: set index j to the second to last non-zero entry of This,
		 */
		{
			j--;

			while (This.at(j)==0 && j != 0 )
			{
				j--;
			}
		}

		if (This.at(j) > 0)
		{
			std::copy(This.begin(), This.begin()+j, Next.begin());
				Next.at(j)   = This.at(j) -1;
				Next.at(j+1) = Next.at(j+1) +1 + This.at(L-1);
				std::copy( Next.begin(), Next.end(), This.begin() );
				std::fill( Next.begin(), Next.end(), 0 );

			this -> iterate();
//				initial_states.insert(  temp_ite, std::pair<long long int, int>
//					( this->to_Decimal(This), -1+this->states_created )  );
			initial_states.emplace(this->to_Decimal(This), -1+ this->states_created);
		}


		if (This.at(L-1) == N)
			NotLast = false;

		j++;
			if (j < L-1)
			{
				NotEnd = true;

			}
	} // ends while loop: Basis_Generation.

	std::cout << "States_Created: "
			  << this->iterator()
			  << " ?= "
			  << this->get_Num_States()
			  << std::endl;

//	this ->print_Ini_States();

	this->reset();
}

long long int Basis_Representation::double_translate(long long int input)
{
	int output = 0;
	int temp = (int) input;
	Basis_Representation::Occupation_Number temp_ONum = this->to_ONum(temp);

	output = this->to_Decimal(temp_ONum);

	return output;
}

Basis_Representation::Occupation_Number Basis_Representation::double_translate (Basis_Representation::Occupation_Number input)
{
	long long int temp_int;
	temp_int = this -> to_Decimal(input);
	return this->to_ONum(temp_int);
}

void Basis_Representation::print_Ini_States()
{
	std::cout << "printing initial states: " << std::endl;
	std::map<long long int, int>::iterator temp_ite
	                 = this->initial_states.begin();

	int temp_count = 0;

	while (temp_ite != this->initial_states.end())
	{
		this->print( this->to_ONum(temp_ite->first) );
		this->print( temp_ite->second);
		std::cout << std::endl;
		temp_ite++;
		temp_count++;
	}

	std::cout << "states printed: "
			  << temp_count
			  << " ?= "
			  << this->get_Num_States()
			  << std::endl;
}

int Basis_Representation::get_Particle_at(int site, long long int input)
{
//	int L = this->Lattice_Len;

	long long int output = 0;

//	int quotient = 0,
//	   remainder = 0;
//	lldiv_t     temp;
//
//	long long int dividend = input,
//		           divisor = (long long int) (1 + this->Particle_Num);
//
//		      int    digit = 0;
//
//	while (digit != L-1)
//	{
//			 temp = lldiv (dividend, divisor);
//
//		 quotient = temp.quot;
//		remainder =  temp.rem;
//
//		if (digit == site)
//			output = remainder;
//
//		digit++;
//
////		std::cout << digit << std::endl;
//
//		dividend  = quotient;
//		quotient  = 0;
//		remainder = 0;
//	}

	Basis_Representation::Occupation_Number ONum_input;

	ONum_input = this->to_ONum(input);

	output = ONum_input.at(site);

	return output;
}

std::vector<std::string> Basis_Representation::split(std::string strToSplit, char delimeter)
{
    std::stringstream ss(strToSplit);
    std::string item;
	std::vector<std::string> splittedStrings;
    while (std::getline(ss, item, delimeter))
	{
		splittedStrings.push_back(item);
    }
	return splittedStrings;
}

void Basis_Representation::print_Map(std::map<long long int, std::string> map_to_print)
{
	std::string value_Str, first_token, second_token;
		std::vector <std::string> vec_token;

	bool IsEmpty = map_to_print.empty();
		if(IsEmpty == true)
		std::cout << "@print_Map input empty?(1) " << IsEmpty << std::endl;


	for	(std::pair <long long int, std::string> element : map_to_print)
	{
		vec_token = this->split(element.second, ',');
			first_token  = vec_token.at(0);
			second_token = vec_token.at(1);

		std::cout << "initial state: ";
//			this->print( this->to_ONum(element.first) );
			std::cout << "final state: ";
//			this->print(  this->to_ONum( std::stoi(first_token) )  );
			std::cout << "coefficient: "
					  << std::stod(second_token)
					  << std::endl;
	}
}

std::map<long long int, std::string> Basis_Representation::annihilate_at (int entry, std::map<long long int, std::string> input_states)
{

		int     N = this->Particle_Num;
		int      		  base = 1 + N;

		int   		 temp_ptle_num = 0;
	long long int        temp_key1 = 0;
	long long int        temp_key2 = 0;
	 double       		temp_coeff = 0;

	std::string   out_Str, in_Str, first_token, second_token;

	std::vector <std::string> vec_token;

	std::map <long long int, std::string> output;

	for	(
		std::pair <long long int, std::string>
			element : input_states
	)
	{
		std::cout << "@annihilate at @begin for loop" << std::endl;
			in_Str = element.second;
			vec_token   = this->split(in_Str, ',');
			first_token = vec_token[0];
			temp_ptle_num = this->get_Particle_at ( entry, std::stoi(first_token) );
//			std::cout << "particle number at " << entry << " => " << temp_ptle_num << std::endl;
		if (temp_ptle_num < (N+1))
		{
			if (temp_ptle_num > 0)
			{
				second_token = vec_token[1];
					temp_coeff = std::stod(second_token);
					temp_coeff = temp_coeff * std::sqrt( temp_ptle_num );

				temp_key1 = element.first;
					temp_key2 = std::stoi(first_token) - (long long int) std::pow(base, entry);

				out_Str = std::to_string(temp_key2)
					+ "," + std::to_string(temp_coeff);
					output.insert(std::pair<long long int, std::string> (temp_key1, out_Str));

//				std::cout << "@annihilate_at [overload] "
//						  << "key: " << temp_key1
//						  << "value: " << out_Str
//						  << std::endl;

				temp_coeff = 0;
					temp_key1=0;
					temp_key2=0;
					second_token.clear();
					out_Str.clear();
			}
		}

		first_token.clear();
			in_Str.clear();
			vec_token.clear();
			temp_ptle_num=0;
	} // for loop ends.

	return output;
}

std::map<long long int, std::string> Basis_Representation::annihilate_at (int entry) // @suppress("No return")
{
//	Initial_States::iterator temp_ite
//    = this->initial_states.begin();

	          int          N = this->Particle_Num;
	          int       base = 1+ N;

	          int temp_ptle_num = 0;
	long long int     temp_key1 = 0;
	long long int     temp_key2 = 0;
		   double    temp_coeff = 0;

	  std::string   temp_Str = "";

	  std::map<long long int, std::string> output;

	for	(
			std::pair<long long int, int> element :
			this->initial_states
		)
		{
			temp_ptle_num = this->get_Particle_at(entry, element.first);
			if ( temp_ptle_num > 0 )
			{
				temp_key1 = element.first;
					temp_key2 = temp_key1
							  - (long long int) std::pow(base, entry);

					temp_coeff = std::sqrt( temp_ptle_num );

				/**
					append the pair (temp_key, temp_coeff)
					into final_states[key=key(initial_states), value=(temp_key, temp_coeff)]
				*/
				temp_Str  = std::to_string(temp_key2)
						  + ","
						  + std::to_string(temp_coeff);
					output.insert((std::pair<long long int, std::string> (temp_key1, temp_Str)));

//					std::cout << " transforming: ";
//					this->print( this->to_ONum(temp_key1) );
//					std::cout << " => ";
//					this->print( this->to_ONum(temp_key2) );
//					std::cout << std::endl;

				temp_key1 = 0;
					temp_key2 = 0;
					temp_coeff = 0;
					temp_Str.clear();
			}

			temp_ptle_num = 0;
		} // for loop ends.

	std::cout << "@annihilate_at: " << std::endl;
//		this->print_Map(output);

	return output;
}

std::map<long long int, std::string> Basis_Representation::create_at (int entry, std::map<long long int, std::string> input_states)
{
//	Initial_States::iterator temp_ite
//    = this->initial_states.begin();

	          int          N = this->Particle_Num;
	          int       base = 1 + N;

	          int    temp_ptle_num = 0;
	long long int        temp_key1 = 0;
	long long int        temp_key2 = 0;
		   double       temp_coeff = 0;

	  std::string   out_Str, in_Str, first_token, second_token;

	  std::vector <std::string> vec_token;

	  std::map <long long int, std::string> output;

	for	(
			std::pair <long long int, std::string>
				element : input_states
		)
		{
			in_Str = element.second;
				vec_token   = this->split(in_Str, ',');
				first_token = vec_token[0];
				temp_ptle_num = this->get_Particle_at ( entry, std::stoi(first_token) );

			if (temp_ptle_num < N)
			{

				second_token = vec_token[1];
					temp_coeff = std::stod(second_token);
					temp_coeff = temp_coeff * std::sqrt(temp_ptle_num + 1);

				temp_key1 = element.first;
					temp_key2 = std::stoi(first_token) + (long long int) std::pow(base, entry);

				out_Str = std::to_string(temp_key2)
					+ "," + std::to_string(temp_coeff);
					output.insert(std::pair<long long int, std::string> (temp_key1, out_Str));

				temp_coeff = 0;
					temp_key1=0;
					temp_key2=0;
					second_token.clear();
					out_Str.clear();
			}

			first_token.clear();
				in_Str.clear();
				vec_token.clear();
				temp_ptle_num=0;
		} // for loop ends.

	return output;
}

void Basis_Representation::make_1PT_Operator ()
{
	// generate all Matrix Rep for tight-binding 1PT operators
	int L = this -> Lattice_Len;
	int OnePTCoupling = (-1) * this->t;


	for(int m = 0; m < L-1; m ++)
	{
		this->set_Final_States( this->hop(m,m+1) );
		this->add_TList(  this->to_Triplet_List( OnePTCoupling )  );

		this->final_states.clear();
	}

	for(int n = 1; n < L; n ++)
	// for_loop: print corresponding hermitian conjugates
	{
		this->set_Final_States( this->hop(n,n-1) );
		this->add_TList(  this->to_Triplet_List( OnePTCoupling )  );

		this->final_states.clear();
	}
}

std::map<long long int, std::string> Basis_Representation::hop(int i, int j)
{
	std::map<long long int, std::string> intermediate, output;
		intermediate = this->annihilate_at(i);
		output = this->create_at(j, intermediate);

	std::cout << "@hop print_Map" << std::endl;
//		this->print_Map(output);

	return output;
}

std::map<long long int, std::string> Basis_Representation::number_Operator(int entry)
{
	std::map<long long int, std::string> intermediate, output;
		intermediate = this->annihilate_at(entry);
		output = this->create_at(entry, intermediate);

	std::cout << "@number_Operator print_Map" << std::endl;
//		this->print_Map(output);

	return output;
}

std::map<long long int, std::string> Basis_Representation::number_Operator (int entry, std::map<long long int, std::string> in_states)
{
	std::cout << "@begin number_Operator[overload] in_states" << std::endl;
//	this->print_Map(in_states);

	std::map<long long int, std::string> intermediate, output;

	intermediate = this->annihilate_at(entry, in_states);
//	this->print_Map(intermediate);
	output = this->create_at(entry, intermediate);

	std::cout << "@number_Operator[overload] print_Map" << std::endl;
//		this->print_Map(output);

	return output;
}

std::map<long long int, std::string> Basis_Representation::add_constant (int constToAdd,
		std::map<long long int, std::string> statesToAdd)
{
	std::map<long long int, std::string> statesOut;

	std::vector <std::string> vec_token;
	std::string token1, tokenOutCoeff;
	long long int keyOut;

	int strToInt = 0;

	for	( std::pair <long long int, std::string> element : statesToAdd )
	{
		vec_token = this->split(element.second, ',');
			token1 = vec_token[0];
			keyOut = element.first;
			std::string strOut;

		strToInt = std::stoi(vec_token[1]);
			strToInt = constToAdd + strToInt;
			std::cout << "@add_constant key: "<< element.first
					  << " => n+1: "<< strToInt << std::endl;

		tokenOutCoeff = std::to_string(strToInt);

		strOut = std::to_string(keyOut) + "," + tokenOutCoeff;

		statesOut.insert(std::pair<long long int, std::string>
						(keyOut, strOut) );

		vec_token.clear();
		token1.clear();
		strOut.clear();
		tokenOutCoeff.clear();
		strToInt = 0;
		keyOut = 0;
	}

	return statesOut;
}

void Basis_Representation::make_2PT_Operator ()
{
	std::cout << "@begin make_2PT" << std::endl;

	std::map<long long int, std::string> intermediate, intermediate1, statesToStore;

	int L = this -> Lattice_Len;
	int TwoPTCoupling = (+1.0)/2.0 * this->J;

	for (int entry = 0; entry < L; entry ++)
	{
		intermediate  = this->number_Operator(entry);
		intermediate1 = this->add_constant(-1, intermediate);
		statesToStore = this->number_Operator(entry, intermediate1);
		this->set_Final_States( statesToStore );
		this->add_TList(  this->to_Triplet_List( TwoPTCoupling )  );

		this->final_states.clear();
		intermediate.clear();
		intermediate1.clear();
		statesToStore.clear();
	}

	std::cout << "@end make_2PT" << std::endl;
}

void Basis_Representation::set_Final_States(std::map<long long int, std::string> states_to_put)
{
	if (this->final_states.empty() != 1)
	{
		std::cout<< "Basis_Representation variable final_states already filled" << std::endl;
		return;
	}

	this -> final_states.insert( states_to_put.begin(), states_to_put.end() );
	std::cout << "@set_Final_States: " << std::endl;
//	this -> print_Map(this->final_states);
}

std::map<long long int, std::string> Basis_Representation::get_Final_States()
{
	return this->final_states;
}

void Basis_Representation::print_TList(std::vector< Eigen::Triplet<double> > TList_to_print)
{
//	for (
//			std::vector < Eigen::Triplet<double> >
//			element : TList_to_print
//		)
//	{
//		std::cout << element;
//	}
	for(std::size_t i=0; i<TList_to_print.size(); ++i)
	{
		std::cout <<TList_to_print[i].col() << " "
				  <<TList_to_print[i].row() << " "
				  <<TList_to_print[i].value()
				  <<std::endl;
	}
}

std::vector< Eigen::Triplet<double> > Basis_Representation::to_Triplet_List(int couplingToPut)
{
	std::vector< Eigen::Triplet<double> > output_TList;

	long long int initial_key = 0;
	long long int final_key = 0;
	double coeff_to_put = 0;
	double coupling = couplingToPut;
	int row_index = 0; // row runs thru. initial_states
	int col_index = 0; // column runs thru. final_states
	std::string Str_values, first_token, second_token;
	std::vector<std::string> vec_tokens;


	for	(
			std::pair<long long int, std::string> element :
			this->final_states
		)
	{
		Str_values = element.second;
			vec_tokens   = this->split(Str_values,',');
			first_token  = vec_tokens[0];
			second_token = vec_tokens[1];
			initial_key  = element.first;

		row_index = this->initial_states.at(initial_key);

		final_key = std::stoi(first_token);
			col_index = this->initial_states.at(final_key);

		coeff_to_put = coupling * std::stod(second_token);

		output_TList.push_back
			( Eigen::Triplet<double> (row_index, col_index, coeff_to_put) );
	}

	std::cout << "TList:" << std::endl;
		this -> print_TList( output_TList );
		std::cout << std::endl;

	return output_TList;
}

void Basis_Representation::set_TList (int couplingToPut)
{
	this->TList.clear();
//	long long int initial_state = 0;
//	std::vector<std::string> vec_tokens;
//		std::string Str_value, first_token, second_token;
//
//	Eigen::SparseMatrix<double> hamiltonian;
	std::vector<T> triplet_List;
		int MDim = this->final_states.size();

	triplet_List.reserve ( MDim );
		this->TList = this->to_Triplet_List( couplingToPut );
//		this->print_TList(triplet_List);
		std::cout << std::endl;

	this->final_states.clear();
//		std::cout << "@set_TList: final states empty?(1) "
//				  << this->final_states.empty()
//				  << std::endl;
	std::cout << "@end set_TList" << std::endl;
}

void Basis_Representation::add_TList (std::vector< Eigen::Triplet<double> > TListToAdd)
{
	for(std::size_t i=0; i<TListToAdd.size(); ++i)
	{
		this->TList.push_back( TListToAdd.at(i) );
	}

	std::cout << "@add_TList: " << std::endl;
//	this->print_TList(this->TList);
}

int Basis_Representation::get_Ini_States_size()
{
	return this->initial_states.size();
}

Eigen::SparseMatrix<double> Basis_Representation::make_Hamiltonian()
{
	this->make_2PT_Operator();
	this->make_1PT_Operator();
	this->make_Chem_Pot_Ops();

	int MDim = this->initial_states.size();
		Eigen::SparseMatrix<double> hamiltonian (MDim, MDim);
		hamiltonian.setFromTriplets(this->TList.begin(), this->TList.end());

		std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
//		std::cout << "Basis Vectors for Hamiltonian Marix" << std::endl;
		this->print_Ini_States();

//		std::cout << "@end Basis Vectors for Hamiltonian Matrix" << std::endl;

	return hamiltonian;
}

int Basis_Representation::getJ()
{
	return this->J;
}

int Basis_Representation::getT()
{
	return this->t;
}

void Basis_Representation::make_Chem_Pot_Ops()
{
	std::cout << "@begin make_Chem_Pot " << std::endl;

	int L = this -> Lattice_Len;
	auto ChemPotCoupling = (double) ( (-1) * this->mu );

	std::map<long long int, std::string> statesToStore;

	for (int entry = 0; entry < L; entry ++)
		{
			statesToStore  = this->number_Operator( entry );
			this->set_Final_States( statesToStore );
			this->add_TList(  this->to_Triplet_List( ChemPotCoupling )  );

			this->final_states.clear();
			statesToStore.clear();
		}

	std::cout << "@end make_Chem_Pot " << std::endl;
}

int Basis_Representation::coord_to_i(Basis_Representation::coordinate ctoi)
{
	int x = std::get<0>(ctoi);
	int y = std::get<1>(ctoi);

	int xLen = this->rowLen;

	int intOut;

	intOut = xLen * (y-1) + x;
	/**
	    intToOut is a bijection
	    from the coordinates of a square lattice.
	 */

	return intOut;
}

Basis_Representation::coordinate Basis_Representation::i_to_coord(int itoc)
{
	Basis_Representation::coordinate coordOut;

	int xOut, yOut;

	int xLen = this->rowLen,
		yLen = this->colLen;

	yOut = itoc / xLen;
	xOut = itoc % xLen;
	if (xOut == 0)
		xOut = xOut + xLen;
	if (yOut < yLen)
		yOut ++;


	coordOut = std::make_tuple(xOut,yOut);

	return coordOut;
}

int Basis_Representation::make_Coordinate(int x, int y)
{
	/**
	   make_Coordinate tests coord_to_i && i_to_coord
	 */

	Basis_Representation::coordinate coordOut, itoc;

	int ctoi;

	itoc = std::make_tuple(x,y);

	ctoi = this->coord_to_i(itoc);
	coordOut = this->i_to_coord(ctoi);

	std::cout << "@make_Coordinate coordinate: ("
			  << x << ", " << y << ")"
			  << " => " << ctoi
			  << " => ("
			  << std::get<0>(coordOut)
			  << ","
			  << std::get<1>(coordOut)
			  << ")"
			  << std::endl;

	return ctoi;
}








int  main ()
{
	int N = 1,    // total number of particles
		rowLen = 2,    // length/#-sites of the row.
		colLen = 2,    // of the column.
		D = 2,	  // lattice dim.

		sites = rowLen * colLen,
		// number of total sites
		// for a lattice of square geometry.

	    t  = 1,
		U  = 10,
		mu = 100;

	/**
		test variables
	 */

		int x = 1, y = 2;
		/**
		   coordinate domain [1,infty) x [1,infty)
		**/


		std::tuple<int, int> testCoord(x,y);

	int matDim = 0; // initialization: size of Sparse Matrix

	std::map<long long int, std::string>  states, states1;
	std::vector< Eigen::Triplet<double> > temp_TList;

	Eigen::SparseMatrix<double> SMat1 (matDim, matDim), SMat2 (matDim, matDim);
	Eigen::MatrixXd dense_mat;

	Basis_Representation me (N,sites,t,U,mu,rowLen,colLen);

	SMat1 = me.make_Hamiltonian();

//	dense_mat = Eigen::MatrixXd(SMat1);
//		std::cout << "@main hamiltonian matrix: " << std::endl
//				  << dense_mat << std::endl;
//
//	Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > eigen_sol( dense_mat );
//		std::cout << "eigenvalues: " << std::endl
//				  << eigen_sol.eigenvalues() << std::endl;
//
//
//
//	std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n";
//
//
//
//

	/**
		Aug 23
		Created CList
		Need to write i_to_coord, coord_to_i,
		Constructor includes make_CList(?)
	 */


	me.make_Coordinate(x,y);





















    return 0;
} // main.cpp ends here.


























