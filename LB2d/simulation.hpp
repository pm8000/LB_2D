/**
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>

namespace lb {

/**
 *  @brief Simulation class implementing LB
 *
 *  This class holds a lattice as member (see @ref simulation::l) and
 *  carries out the simulation steps on top of it. The main methods of
 *  this class are @ref simulation::advect() and
 *  @ref simulation::collide().
 */
class simulation
{
public: // ctor

	/**
	 *  @brief Construct from domain size and flow parameters
	 *  @param[in] nx    extent in x direction
	 *  @param[in] ny    extent in y direction
	 *  @param[in] _Re   Reynolds number
	 *  @param[in] _Vmax mean flow velocity
	 */
	simulation(unsigned int nx, unsigned int ny, float_type _Re, float_type _Vmax)
	: l(nx, ny),
	  shift(velocity_set().size),
	  Re(_Re),
	  Vmax(_Vmax),
	  visc( /*fill in your code here*/ 0.001),
	  beta( 1/(2*visc/(this->cs*this->cs)+1)), //EDITED
	  time(0),
	  file_output(false), // set to true if you want to write files
	  output_freq(100),
	  output_index(0)
	{
		// define amount to shift populations for advection
		for (unsigned int i=0; i<velocity_set().size; ++i)
		{
			// **************************
			// * fill in your code here *
			int shift_x = velocity_set().c[0][i];
			int shift_y = velocity_set().c[1][i];
			// **************************
			shift[i] = shift_y * (nx+2) + 1 + shift_x; //EDITED
		}
	}

	/**
	 *  @brief Initialize the flow field
	 *
	 *  Initialization includes defining initial density, velocity and
	 *  populations. You can use Taylor-Green vortex flow conditions.
	 */


	void initialize() //EDITED
	{
		// **************************
		// * fill in your code here *
		// * the lines below are    *
		// * just examples          *
		// **************************
		const float_type pi(std::acos(-1.0));

		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
			    double lambda_x = 1; //wave length can be modified here
                double lambda_y = 1; //wave length can be modified here
                double K_x = 2*pi/(lambda_x*l.nx);
                double K_y = 2*pi/(lambda_y*l.ny);

			    //calculate analytical solution for t=0
				l.get_node(i,j).u()   = -Vmax*std::sqrt(K_y/K_x)*std::sin(K_y*j)*std::cos(K_x*i);
				l.get_node(i,j).v()   = Vmax*std::sqrt(K_x/K_y)*std::sin(K_x*i)*std::cos(K_y*j);
				l.get_node(i,j).rho() = 1.0 - Vmax*Vmax/(2*cs*cs*std::sqrt(K_x*K_x+K_y*K_y))*(K_y*K_y*std::cos(3*K_x*i) + K_x*K_x*std::cos(2*K_y*j));
			}
		}
	}

	/**
	 *  @brief advect the populations
	 *
	 *  Include periodic boundary conditions here also
	 */
	void advect()
	{
		// **************************
		// * fill in your code here *
		// **************************

	}

	/**  @brief apply wall boundary conditions */
	void wall_bc()
	{
		#pragma omp parallel for
		for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
		{
			// **************************
			// * fill in your code here *
			// **************************
		}
	}

	/** @brief collide the populations */
	void collide()
	{
		// **************************
		// * fill in your code here *
		// **************************

	}

	/** @brief LB step */
	void step()
	{
		advect();
		wall_bc();
		collide();

		// file io
		if ( file_output && ( ((time+1) % output_freq) == 0 || time == 0 ) )
		{
			write_fields();
			++output_index;
		}

		++time;
	}

public: // write to file

	/** write macroscopic variables to ascii file */
	void write_fields()
	{
		std::stringstream fns;
		fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
		l.write_fields(fns.str());
	}

public: // print

	/** print to output stream */
	friend std::ostream& operator<<(std::ostream& os, const simulation& sim)
	{
		os << "simulation parameters\n"
		   << "---------------------\n";
		os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
		os << "Re:     " << sim.Re << "\n";
		os << "Vmax:   " << sim.Vmax << "\n";
		os << "visc:   " << sim.visc << "\n";
		os << "beta:   " << sim.beta << "\n";
		return os;
	}

public: // members

	lattice l;                 ///< lattice
	std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
	const float_type Re;       ///< Reynolds number
	const float_type Vmax;     ///< mean flow velocity
	const float_type visc;     ///< viscosity
	const float_type beta;     ///< LB parameter beta
	const float_type cs = 1/std::sqrt(3); //speed of sound
	unsigned int time;         ///< simulation time
	bool file_output;          ///< flag whether to write files
	unsigned int output_freq;  ///< file output frequency
	unsigned int output_index; ///< index for file naming
};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED