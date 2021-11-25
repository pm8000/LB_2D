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
	  beta(1/(2*visc/(velocity_set().cs*velocity_set().cs)+1)), //EDITED
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
			shift[i] = shift_y * (nx+2) + shift_x; //EDITED
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
				velocity_set().equilibrate(l.get_node(i,j));
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

		//EDITED

		//advect populations towards top and rights
		//std::cout<<"before advection/BC"<<std::endl;
		//std::cout<<l<<std::endl;

		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
                for (unsigned int k=1; k<shift.size(); ++k)
                {
                    if (k<3 || k==5 || k==6)
                    {
                        continue;
                    }
                    l.get_node(l.index(i,j)+shift[k]).f(k)=l.get_node(i,j).f(k);
                }
			}
		}
		//advect populations towards bottom and left
		for (int j=static_cast<int>(l.ny)-1; j>-1; --j)
		{
			for (int i=static_cast<int>(l.nx)-1; i>-1; --i)
			{
                for (unsigned int k=1; k<shift.size(); ++k)
                {
                    if (k==3 || k==4 || k>=7)
                    {
                        continue;
                    }
                    l.get_node(l.index(i,j)+shift[k]).f(k)=l.get_node(i,j).f(k);
                }
			}
		}
	}

	/**  @brief apply wall boundary conditions */
	void wall_bc()
	{
		#pragma omp parallel for
			// **************************
			// * fill in your code here *
			// **************************
			//EDITED
			//bottom and top walls
			for (int i=0; i<static_cast<int>(l.nx);++i)
            {
                l.get_node(i,0).f(2)=l.get_node(i,-1).f(4);
                l.get_node(i,static_cast<int>(l.ny)-1).f(4)=l.get_node(i, static_cast<int>(l.ny)).f(2);
                l.get_node(i,0).f(5)=l.get_node(i-1,-1).f(7);
                l.get_node(i,0).f(6)=l.get_node(i+1,-1).f(8);
                l.get_node(i,static_cast<int>(l.ny)-1).f(7)=l.get_node(i+1, static_cast<int>(l.ny)).f(5);
                l.get_node(i,static_cast<int>(l.ny)-1).f(8)=l.get_node(i-1, static_cast<int>(l.nx)).f(6);

            }
            //left and right walls
            for(int j=0; j<static_cast<int>(l.ny);++j)
            {
                l.get_node(0,j).f(1)=l.get_node(-1,j).f(3);
                l.get_node(static_cast<int>(l.nx)-1,j).f(3)=l.get_node(static_cast<int>(l.nx), j).f(1);
                l.get_node(0,j).f(5)=l.get_node(-1,j-1).f(7);
                l.get_node(static_cast<int>(l.nx)-1,j).f(6)=l.get_node(static_cast<int>(l.nx), j-1).f(8);
                l.get_node(static_cast<int>(l.nx)-1,j).f(7)=l.get_node(static_cast<int>(l.nx), j+1).f(5);
                l.get_node(0,j).f(8)=l.get_node(-1,j+1).f(6);
            }
            //std::cout<<"after advection/BC"<<std::endl;
            //std::cout<<l<<std::endl;
	}

	/** @brief collide the populations */
	void collide()
	{
		// **************************
		// * fill in your code here *
		// **************************
		//EDITED
		for (int j=0; j<static_cast<int>(l.ny); ++j)
		{
			for (int i=0; i<static_cast<int>(l.nx); ++i)
			{
			    double rho_cur=0;
			    double u_cur=0;
			    double v_cur=0;
                for (unsigned long int k=0; k<shift.size(); ++k)
                {
                    rho_cur+=l.get_node(i,j).f(k);
                    u_cur+=l.get_node(i,j).f(k)*velocity_set().c[0][k];
                    v_cur+=l.get_node(i,j).f(k)*velocity_set().c[1][k];
                }
                l.get_node(i,j).rho()=rho_cur;
                l.get_node(i,j).u()=u_cur;
                l.get_node(i,j).v()=v_cur;
                float_type f[9];
                velocity_set().f_eq(f, rho_cur, u_cur, v_cur);
                for (unsigned long int k=0; k<shift.size(); ++k)
                {
                    l.get_node(i,j).f(k) += 2*beta*(f[k]-l.get_node(i,j).f(k));
                }

			}
		}
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
