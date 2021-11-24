
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>



int main(int argc, char *argv[])
{
	omp_set_num_threads(std::max(omp_get_max_threads(),omp_get_num_procs()));
	lb::simulation* sim = new lb::simulation(10,10,10000,0.05);
	sim->initialize();
	std::cout << *sim << std::endl;

	#ifdef USE_OPENGL_VISUALIZATION

        lb::visualization::initialize(sim,argc,argv);
		lb::visualization::get_instance().run();

	#else
        std::cout<<"else"<<std::endl;
		// Here are some hints for getting aquainted with the lattice class
		// ================================================================

		// how to print the lattice:
		// -------------------------

		std::cout << sim->l << std::endl;
		std::cout<<"first statement printed"<<std::endl;

		// how to access the lattice:
		// --------------------------

		// 1) access via node proxy
		sim->l.get_node(1,0).f(0) = 2;

		// 2) access data directly (make sure you know what you're doing)
		sim->l.f[0][sim->l.index(2,0)] = 3;

		// 3) using iterators to nodes
		(sim->l.begin() + sim->l.index(0,0))->f(0) = 1;

		//new
		std::cout<<"got here"<<std::endl;
        sim->step();

		std::cout << sim->l << std::endl;




		// use a loop like this to run the simulation

		/*for (unsigned int i=0; i<500; ++i)
		{
			sim->step();
		}*/

	#endif
    std::cout<<"end"<<std::endl;
	return 0;
}
