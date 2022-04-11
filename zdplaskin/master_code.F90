! Fortran script for running calculation on Ar + H2O plasma
! Runs calculation in a given X (position) for given simulation time and timestep

program main
	use options
	use ZDPlasKin
	implicit none
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!! initialize program !!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	

	! simulation run variables (updating on each step) - only init without values
	double precision:: t						! time (seconds) - original = time
	double precision:: right_sheath_thickness 	! thickness of the right sheath (cm) - original Sa
	double precision:: left_sheath_thickness 	! thickness of the left sheath (cm) - original Sb
	double precision:: bulk_start_position 		! position of begining of bulk plasma (after right sheath) (cm) - original Xa
	double precision:: bulk_end_position 		! position of the end of bulk plasma (before left sheath) (cm) - original Xb
	double precision:: Vab 						! ????
	double precision:: electron_density 		! electron density in bulk - original E
	! position dependent vectors
	double precision:: electric_field 			! electric field vector (position dependent) - original Eg
	double precision:: reduced_electric_field	! reduced electric field (Td) - original EN
	double precision:: omega_pe					! ????
	double precision:: electric_potential		! electric potential in the plasma - original Potential
	double precision:: all_neutral				! ????

	! physical constants
	double precision, parameter:: pi = 3.141593d+0 					! pi = 3.14... - original Pi
	double precision, parameter:: electron_charge = 1.60d-19 		! electron charge (C) - original ele
	double precision, parameter:: electron_mass = 9.109383d-31 	! electron mass (kg) - original m
	double precision, parameter:: permittivity_const = 8.85d-14	! permittivity of free space (F/cm) - original eps_0

	! simulation parameters - command line arguments
	double precision:: total_simulation_time	! total time (sec) to run simulation for
	double precision:: simulation_timestep		! time duration (sec) for single simulation step
	double precision:: x						! position (cm) to run the simulation for
	double precision:: initial_electron_density ! initial electron density - original Initial_elec

	! constants calculated from simulation parameters
	double precision:: rf_radial_frequency		   ! rf frequency in (rad / sec) - original Circular_frequency
	double precision:: rf_max_voltage			   ! maximal voltage of RF signal - original V0
	double precision:: rf_max_current 			   ! maximal current of RF signal - original I0
	double precision:: mean_sheath_thickness	   ! mean sheath thickness - original s0					
	double precision:: bulk_length				   ! bulk length (cm) - original d
	double precision:: max_sheath_thickness		   ! maximal thickness of the sheath (right + left) (cm) - original Sm
	double precision:: electric_field_distribution ! paramter for electric field distribution - original ea
	double precision:: Vp						   ! some constant parameter related to "distribution of the electric field"
	double precision:: plasma_dielectric_const 	   ! dielectric constant of the plasma - original eps_p
	integer:: j ! iteration parameter

	! parsing command line argumnet (run parameteres)
	call parse_arg(1, total_simulation_time)
	call parse_arg(2, simulation_timestep)
	call parse_arg(3, x)
	call parse_arg(4, initial_electron_density)	
	call set_parameters()

	! setting up additional variable values (must be done after declaration for some reason)
	rf_radial_frequency = 2 * pi * rf_frequency
	rf_max_voltage = (rf_power * rf_impedance)**0.50d0 * 2**0.50d0 
	rf_max_current = rf_max_voltage / rf_impedance
	mean_sheath_thickness = rf_max_current / (electron_charge * initial_electron_density * gap_area * rf_radial_frequency)
	bulk_length = gap_length - mean_sheath_thickness
	max_sheath_thickness = 2 * mean_sheath_thickness
	electric_field_distribution = initial_electron_density * electron_charge / permittivity_const
	Vp = 2 * electric_field_distribution * mean_sheath_thickness**2
  omega_pe = (electron_charge**2 *(initial_electron_density)/(permittivity_const * electron_mass) * 1.0d4)**0.50d0
	plasma_dielectric_const = permittivity_const * (1 - (omega_pe / rf_radial_frequency)**2)
	
	! ZDPlasKin initialization
	call ZDPlasKin_init()
	call ZDPlasKin_set_conditions(GAS_TEMPERATURE=gas_temperature, GAS_HEATING=.TRUE.)
	call ZDPlaskin_set_config(ATOL=1D-10,RTOL=1D-5)  

	! initialize concentration vector for t = 0
	do j=1,species_max
	   select case (species_name(j))          
	   case('AR') 
		call ZDPlasKin_set_density(species_name(j), (1 - water_fraction) * gas_density)
	   case('E') 
		  call ZDPlasKin_set_density(species_name(j), initial_electron_density)
	   case('H2O') 
		  call ZDPlasKin_set_density(species_name(j),gas_density * water_fraction)
	   end select
	end do

	! open the results (temp) file and write headers
	open(1,file=filename)
	write(1, "(*(G0, :, ','))") 'time', 'electric_field', 'reduced_electric_field', 'omega_pe', &
		'electric_potential','gas_temperature', (trim(species_name(j)), j=1,species_max)

	! initializing values for t = 0
	t = 0.0
	left_sheath_thickness = mean_sheath_thickness * (1 - cos(rf_radial_frequency * t))
	right_sheath_thickness = max_sheath_thickness - left_sheath_thickness
	bulk_start_position = left_sheath_thickness
	bulk_end_position = gap_length - right_sheath_thickness
	Vab = Vp * cos(rf_radial_frequency * t)
	
   call calculate_electric_potential(x, &
									   electric_field_distribution, &
									   permittivity_const, &
									   plasma_dielectric_const, &
									   bulk_start_position, &
									   bulk_end_position, &
									   gap_length, &
									   Vab, &
									   electric_field, &
									   electric_potential)
	reduced_electric_field = abs(electric_field / gas_density * 1.0d17)


	 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!! main program loop !!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! calculates the concentrations for each time T in the position X
	do while(t .lt. total_simulation_time)
		 ! defining location and charge parameters for time
		 left_sheath_thickness = mean_sheath_thickness * (1 - cos(rf_radial_frequency * t))
		 right_sheath_thickness = max_sheath_thickness - left_sheath_thickness
		 bulk_start_position = left_sheath_thickness
		 bulk_end_position = gap_length - right_sheath_thickness
		 Vab = Vp * cos(rf_radial_frequency * t)
		 print *, "t=", t, "out of", total_simulation_time
		 
		call calculate_electric_potential(x, &
											electric_field_distribution, &
											permittivity_const, &
											plasma_dielectric_const, &
											bulk_start_position, &
											bulk_end_position, &
											gap_length, &
											Vab, &
											electric_field, &
											electric_potential)  ! sets the electric field from the equations
		write(1, "(*(G0, :, ','))") t, electric_field, reduced_electric_field, omega_pe, &
						electric_potential, gas_temperature, (density(j), j = 1, species_max)
			
		! send reduced electric field data to ZDPlasKin
		call ZDPlasKin_set_conditions(REDUCED_FIELD = reduced_electric_field)
		! ZDPlasKin timestep progress
		call ZDPlasKin_timestep(t, simulation_timestep)
		! extract plasma properties
		! electron frequency
		electron_density = density(species_electrons)       
		omega_pe = ((electron_charge**2.0d0) * (electron_density) / (permittivity_const * electron_mass) * 1.0d4)**0.50d0
		! calculate gas_temperature
		call ZDPlasKin_get_conditions(GAS_TEMPERATURE = gas_temperature)

		! calculate reduced electric field (for the next iteration)
		call ZDPlasKin_get_density_total(ALL_NEUTRAL = all_neutral)
		reduced_electric_field = (electric_field / all_neutral) * 1.0d17   !!reduced electric field in Td
		! advance 
		t = t + simulation_timestep

		!writing file and screen
	end do
	! writing results for the last time
	write(1, "(*(G0, :, ','))") t, electric_field, reduced_electric_field, omega_pe, &
	electric_potential, gas_temperature, (density(j), j = 1, species_max)

	! closing results file
	close(1)
contains

subroutine calculate_electric_potential(x, &
										electric_field_distribution, &
										permittivity_const, &
										plasma_dielectric_const, &
										bulk_start_position, &
										bulk_end_position, &
										gap_length, &
										Vab, &
										output_electric_field, &
										output_potential)

	double precision,intent(in):: x, &
									electric_field_distribution, &
									permittivity_const, &
									plasma_dielectric_const, &
									bulk_start_position, &
									bulk_end_position, &
									gap_length, &
									Vab
	double precision,intent(out):: output_electric_field, output_potential
	if((x <= bulk_start_position .and. x >= 0) .and. bulk_start_position .ne. 0) then
	   output_electric_field = electric_field_distribution * (x - bulk_start_position)
	   output_potential = - electric_field_distribution * (x**2 / 2 - bulk_start_position * x) + Vab
	else if((x > bulk_start_position .and. x < bulk_end_position) .and. (bulk_start_position .ne. 0)) then
	   output_potential = electric_field_distribution * bulk_start_position**2 / 2 + Vab
	   output_electric_field = abs((permittivity_const / plasma_dielectric_const) * electric_field_distribution * bulk_start_position)
	else if((x >= bulk_start_position .and. x < bulk_end_position) .and. (bulk_start_position == 0)) then
		output_potential = electric_field_distribution * bulk_start_position**2 / 2 + Vab
		output_electric_field = abs((permittivity_const / plasma_dielectric_const) * electric_field_distribution * &
								(gap_length - bulk_end_position))
	else if(x >= bulk_end_position .and. x <= gap_length) then
	   output_potential = -electric_field_distribution * (x**2 / 2 - bulk_end_position * x + bulk_start_position**2 / 2) + &
	   						(electric_field_distribution / 2) * (gap_length - bulk_end_position)**2
	   output_electric_field = electric_field_distribution * (x - bulk_end_position)
	end if
end subroutine calculate_electric_potential
end program main