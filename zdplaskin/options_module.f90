module options
	implicit none
	public
	double precision:: rf_power					! power of the rf signal - original Power
	double precision:: rf_impedance				! impedance of the rf signal - original Impedance
	double precision:: rf_frequency 			! frequency of rf wave (Hz) - original Frequency
	double precision:: gas_pressure				! pressure of the gas in the chamber
	double precision:: gap_length				! ????
	double precision:: radius					! ????
	double precision:: gas_temperature			! ????
	double precision:: water_fraction			! fraction of water from the initial gas mixture
	character(len=99):: filename 				! path to save results
	double precision:: gap_area
	double precision:: gas_density
	! dummy variable to parse cmd arguments
	character(len=100), private:: arg
contains
subroutine set_parameters()
	implicit none
	call parse_arg(5, rf_power)
	call parse_arg(6, rf_impedance)
	call parse_arg(7, rf_frequency)
	call parse_arg(8, gas_pressure)
	call parse_arg(9, gap_length)
	call parse_arg(10, radius)
	call parse_arg(11, gas_temperature)
	call parse_arg(12, water_fraction)
	call get_command_argument(13, filename)
    gap_area    = 3.1415 * radius**2
    gas_density = 1.0d-6 * (101325.0d0 * gas_pressure / 760.0d0) / (1.38065d-23 * gas_temperature)
end subroutine set_parameters
subroutine parse_arg(arg_no, output)
	integer arg_no
	double precision:: output
	! passing command line arg to dummy variable
	call get_command_argument(arg_no, arg)
	! parsing variable to float
	read(arg, *) output
end subroutine parse_arg
end module options
