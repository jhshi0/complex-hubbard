
#ifndef QUANTUM_TOOLBOX_KINEMATICS
#define QUANTUM_TOOLBOX_KINEMATICS


/**
   Kinematics.h should be a header only file
   	   defined by two things:
 			Basis_Representation and Symmetry.
 				Basis_Representation is a class,
 					constructs the basis for the Fock Space, while
 		 				Symmetry is a public function of Basis_Representation
 		 					encodes its symmetries.
 		 						Together, they specify the kinematics for the many body system. */


namespace Quantum_Toolbox {

#include <bits/stdc++.h>
#include <string>


const int Lattice_Len = 2;

/**
	Basis_Representation is a class, which should do the following:
		prompt && store user defined parameters about the second quantized Hamiltonian, and
			use make_Vect_n to define a Basis(list) of all allowed states, and
				define a hash list[or Lin Table] from basis vector to basis vector, and
					output.

			Assumptions: LatDim==1[implicit] statistics s==+1 [implicit,bosons] spin
			*/

class Basis_Representation {



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
			Basis_Representation(int);
			Basis_Representation(int, int, int);



	private:

		/**
			- USER INPUT VARIABLES -
			Particle_Num			Total number of particles
			Lattice_Len				Total number of lattice sites
			t						coupling parameter for 1pt term
			J						coupling parameter for 2pt term

			- INTERNAL VARIABLES -
			basis_vect_n			a list of tuples (0/1 0/1)
			basis_Dim				size of vector of tuples basis_Vect_n
			lin_value				a list of decimal numbers encoding elements in basis_vect_n
			lin_key					a list of decimals encoding the keys for lin_value

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
			set_Int_Var				initialize all internal variables as zeros

			- GETTERS -
			get_Part_Num			returns Particle_Num
			get_Lattice_Len			returns Lattice_Len
			get_t					returns coupling t
			get_J					returns coupling J

			get_Basis_Dim			returns the size of the vector of tuples

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
				 	 	   t,
					 	   J;

			int Basis_Dim;

			//std::vector< std::tuple<int, int> > basis_vector_n;

			void make_Basis_Vect_n();
			void make_Vect_n(int, int);

			void set_Param(int, int, int);
			void set_Int_Var();

			int  get_Basis_Dim();

};

Basis_Representation::Basis_Representation()
{

	set_Param(0,0,0);
}

Basis_Representation::Basis_Representation(int A)
/**
	 @return constructor[number_of_particle, number_of_Lat_Sites]=> set the dimensions for Basis
 */
{
	set_Param(A,0,0);
	set_Int_Var();
}

Basis_Representation::Basis_Representation(int A, int B, int C)
{
	set_Param(A,B,C);
	set_Int_Var();
}

void Basis_Representation::set_Param(int A, int B, int C)
/**
   @ set to zero
 */
{
		Particle_Num = A,

				  t	 = B,
				  J  = C;
}

void Basis_Representation::set_Int_Var()
/**
   @param initialize basis_vect_n
 	 uses void make_Basis_Vect_n
 */
{
	make_Basis_Vect_n();
}



void Basis_Representation::make_Basis_Vect_n()
/**
 	 make_Vect_n
 	 	 @param uses Particle_Num && Lattice_Len
 	 	 	 assumes Tind Binding Model/nearest neighbor hopping
 	 	 	 	 	 writes down all allowed states in occupation basis.
 */
{
	int N = Particle_Num;
	int i = 0, j = N-i;

	for (int count = 0; count <= N; count++)
	{
		// append tuple (i,j) into basis_vector_n
		make_Vect_n(i,j);
		count++;
		i++;
	}
}


//void Basis_Representation::make_Vect_n(int n1, int n2)
/**
	uses pushback && make_tuple to define a tuple.
 */
/**
{
	// user input cannot exceed total Particle Number
	if (n1 > Particle_Num || n2 > Particle_Num)
	{
		std::cout << 'invalid user input at make_Vect_n' << std::endl;
		return;
	}

	basis_vector_n.push_back( std::make_tuple(n1,n2) );
	Basis_Dim++;
} */



















}// namespace Quantum_Toolbox ends here.

#endif /* KINEMATICS_H_ */

