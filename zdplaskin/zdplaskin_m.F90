!
! ZDPLASKIN version 2.0a
! (c) 2008, Sergey Pancheshnyi (pancheshnyi@gmail.com)
!
! BOLSIG+
! (c) 2005, Gerjan Hagelaar (gerjan.hagelaar@laplace.univ-tlse.fr)
!
! http://www.zdplaskin.laplace.univ-tlse.fr/
! This software is provided "as is" without warranty and non-commercial use is freely
! granted provided proper reference is made in publications resulting from its use.
! Use of ZDPlasKin in commerical software requires a license.
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! Mon Jan 31 14:25:32 2022
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! ELEMENTS:  E  H  O AR 
!
!-----------------------------------------------------------------------------------------------------------------------------------
!
! MODULE ZDPlasKin
!
!-----------------------------------------------------------------------------------------------------------------------------------
module ZDPlasKin
  use dvode_f90_m, only : vode_opts
  implicit none
  public
!
! config
!
  integer, parameter :: species_max = 47, species_electrons = 1, species_length = 9, reactions_max = 143, reactions_length = 21
  double precision                          :: density(species_max)
  integer                                   :: species_charge(species_max)
  character(species_length)                 :: species_name(species_max)
  character(reactions_length)               :: reaction_sign(reactions_max)
  logical                                   :: lreaction_block(reactions_max)
!
! internal config
!
  double precision, parameter, private      :: vode_atol = 1.00D-10, vode_rtol = 1.00D-05, cfg_rtol = 1.00D-01
  integer, parameter, private               :: vode_neq = species_max + 1
  type (vode_opts), private                 :: vode_options
  integer, private                          :: vode_itask, vode_istate, ifile_unit = 5
  double precision, private                 :: stat_dens(species_max), stat_src(species_max), stat_rrt(reactions_max), stat_time, &
                                               dens_loc(vode_neq,0:3), rrt_loc(reactions_max), tsav = -huge(tsav), &
                                               mach_accur, mach_tiny
  double precision                          :: rrt(reactions_max), mrtm(species_max, reactions_max), ZDPlasKin_cfg(14)
  logical, private                          :: lZDPlasKin_init = .false., lprint, lstat_accum, ldensity_constant, &
                                               density_constant(species_max), lgas_heating
!
! qtplaskin config
!
  logical, private                          :: lqtplaskin, lqtplaskin_first = .true.
  double precision, parameter, private      :: qtplaskin_atol = 1.00D+00, qtplaskin_rtol = 1.00D-02
  character(32), allocatable                :: qtplaskin_user_names(:)
  double precision, allocatable             :: qtplaskin_user_data(:)
!
! physical constants
!
  double precision, parameter, private      :: eV_to_K = 1.16045052d4, q_elem = 1.60217662d-19, k_B = 1.38064852d-23
!
! bolsig+ config
!
  double precision, parameter, private      :: bolsig_rtol = 1.00D-03, bolsig_rtol_half = 3.16D-02, &
                                               bolsig_field_min = 1.00D-01, bolsig_field_max = 1.00D+03, &
                                               bolsig_eecol_frac_def = 1.00D-05
  double precision, private                 :: bolsig_eecol_frac
  integer, parameter, private               :: bolsig_species_max = 3, bolsig_species_length = 3, bolsig_rates_max = 28 
  character(*), parameter, private          :: bolsigfile = "bolsigdb.dat"
  integer                                   :: bolsig_pointer(bolsig_rates_max) = -1
  integer, private                          :: bolsig_species_index(bolsig_species_max) = -1, bolsig_collisions_max = 0 
  logical, private                          :: lbolsig_ignore_gas_temp, lbolsig_Maxwell_EEDF
  double precision, allocatable             :: bolsig_rates(:)
  character(bolsig_species_length), private :: bolsig_species(bolsig_species_max)
  interface
    subroutine ZDPlasKin_bolsig_Init(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_Init
    subroutine ZDPlasKin_bolsig_ReadCollisions(a)
      character(*), intent(in) :: a
    end subroutine ZDPlasKin_bolsig_ReadCollisions
    subroutine ZDPlasKin_bolsig_GetCollisions(i,j)
      integer, intent(out) :: i, j
    end subroutine ZDPlasKin_bolsig_GetCollisions
    subroutine ZDPlasKin_bolsig_GetSpeciesName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetSpeciesName
    subroutine ZDPlasKin_bolsig_GetReactionName(a,i)
      integer, intent(in) :: i
      character(*), intent(out) :: a
    end subroutine ZDPlasKin_bolsig_GetReactionName
    subroutine ZDPlasKin_bolsig_SolveBoltzmann(i,a,j,b)
      integer, intent(in) :: i, j
      double precision, intent(in)  :: a(i)
      double precision, intent(out) :: b(j)
    end subroutine ZDPlasKin_bolsig_SolveBoltzmann
    subroutine ZDPlasKin_bolsig_GetEEDF(i,a,b)
      integer, intent(in) :: i
      double precision, intent(out) :: a,b
    end subroutine ZDPlasKin_bolsig_GetEEDF
  end interface
!
! data section
!
  data species_charge(1:species_max) &
  /-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1,&
    1,-1, 1, 1,-1/
  data species_name(1:species_max) &
  /"E        ","H2O      ","H2O(J0)  ","H2O(J1)  ","H2O(J2)  ","H2O(J3)  ","H2O(V010)","H2O(V100)","HO2      ","H2O2     ",&
   "O3       ","H(N2)    ","H(N3)    ","H(N4)    ","H        ","H2       ","O2       ","OH       ","OH(A)    ","OH^-     ",&
   "O^-      ","H^-      ","H2O^+    ","H^+      ","OH^+     ","O^+      ","H2^+     ","H3O^+    ","O2^+     ","O(21S)   ",&
   "O(23P)   ","O(33S)   ","O(35S)   ","O(33P)   ","O(35P)   ","O        ","AR       ","AR^+     ","AR*      ","AR2^+    ",&
   "ARH^+    ","AR(W)^+  ","AR2(W)^+ ","E(W)     ","H2O(W)^+ ","OH(W)^+  ","H(W)^-   "/
  data reaction_sign(1:90) &
  /"bolsig:H2O->H2O*     ","bolsig:H2O->H+OH^-   ","bolsig:H2O->H2+O^-   ","bolsig:H2O->HO+H^-   ","bolsig:H2O->H(A)     ",&
   "bolsig:H2O->H(B)     ","bolsig:H2O->H(LA)    ","bolsig:H2O->O(A')    ","bolsig:H2O->O(B')    ","bolsig:H2O->O(1S)    ",&
   "bolsig:H2O->OH(X)    ","bolsig:H2O->OH(A)    ","bolsig:H2O->H2O(ROT) ","bolsig:H2O->H2O(ROT1)","bolsig:H2O->H2O(ROT2)",&
   "bolsig:H2O->H2O(ROT3)","bolsig:H2O->H2O(V010)","bolsig:H2O->H2O(V100)","bolsig:H2O->H2O^+    ","bolsig:H2O->H^+      ",&
   "bolsig:H2O->OH^+     ","bolsig:H2O->O^+      ","bolsig:H2O->H2^+     ","H^-+H2O=>OH^-+H2     ","H2O^++H2O=>H3O^++OH  ",&
   "H2O^++H2=>H3O^++H    ","H2O^++O2=>O2^++H2O   ","OH^++H2O=>H2O^++OH   ","OH^++H=>O^++H2       ","OH^++O=>O^++OH       ",&
   "H^++H2O=>H2O^++H     ","H2^++H2O=>H2O^++H2   ","H2^++H2O=>H3O^++H    ","O^++H2O=>H2O^++O     ","OH^++O2=>O2^++OH     ",&
   "O^-+O=>O2+E          ","O^-+O(21S)=>O2+E     ","O^-+O(23P)=>O2+E     ","O^-+O2=>O3+E         ","O^-+H2=>H2O+E        ",&
   "H^-+H=>H2+E          ","H^-+O2=>HO2+E        ","OH^-+H=>H2O+E        ","OH^-+O=>HO2+E        ","OH^-+O(21S)=>HO2+E   ",&
   "OH^-+O(23P)=>HO2+E   ","E+H2O^+=>H+OH        ","E+H2O^+=>H2+O        ","E+H2O^+=>H+H+O       ","E+H2^+=>H+H          ",&
   "E+OH^+=>O(21S)+H     ","E+H3O^+=>H+H2O       ","E+O2^+=>O+O          ","E+O^+=>O             ","E+E+H2O^+=>H2O+E     ",&
   "E+E+H3O^+=>H+H2O+E   ","H^-+H^+=>H+H(N2)     ","H^-+H^+=>H+H(N3)     ","H^-+H2^+=>H+H2       ","H^-+O^+=>H+O         ",&
   "H^-+OH^+=>H+OH       ","H^-+H2O^+=>H+H2O     ","O^-+O^+=>O+O         ","O^-+O2^+=>O+O2       ","O^-+H^+=>O+H         ",&
   "O^-+H2^+=>O+H2       ","O^-+OH^+=>O+OH       ","O^-+H2O^+=>O+H2O     ","OH^-+OH^+=>OH+OH(A)  ","OH^-+H^+=>OH+H       ",&
   "OH^-+H2^+=>OH+H2     ","OH^-+O^+=>OH+O       ","OH^-+H2O^+=>OH+H2O   ","H^-+H3O^+=>H2+H2O    ","O^-+H3O^+=>H+O+H2O   ",&
   "OH^-+H3O^+=>H2O+H2O  ","H^-+O2^+=>H+O2       ","OH^-+O2^+=>OH+O2     ","OH+OH=>O+H2O         ","OH+H2=>H2O+H         ",&
   "OH+O=>H+O2           ","OH+O=>HO2            ","OH+H=>O+H2           ","OH+H2O2=>HO2+H2O     ","OH+OH=>H2O2          ",&
   "OH+O3=>HO2+O2        ","H+HO2=>OH+OH         ","H+HO2=>H2O+O         ","OH+HO2=>H2O+O2       ","O+HO2=>OH+O2         "/
  data reaction_sign(91:143) &
  /"HO2+H=>H2+O2         ","HO2+HO2=>H2O2+O2     ","HO2+O3=>OH+O2+O2     ","H2O2+O=>HO2+OH       ","H2O2+H=>H2O+OH       ",&
   "H2O2+H=>HO2+H2       ","O+O3=>O2+O2          ","O+H2=>OH+H           ","H+O2+H2O=>HO2+H2O    ","H+O2=>HO2            ",&
   "H+O3=>OH+O2          ","O(21S)+H2O=>O+H2O    ","O(21S)+H2O=>OH+OH    ","O(23P)+H2O=>OH+OH    ","O(21S)=>O            ",&
   "O(33P)=>O(33S)       ","O(35P)=>O(35S)       ","O(33S)=>O(23P)       ","O(35S)=>O(23P)       ","OH(A)=>OH            ",&
   "H(N2)=>H             ","H(N3)=>H(N2)         ","H(N4)=>H(N2)         ","H2O(J0)=>H2O         ","H2O(J1)=>H2O         ",&
   "H2O(J2)=>H2O         ","H2O(J3)=>H2O         ","H2O(V010)=>H2O       ","H2O(V100)=>H2O       ","AR^++H=>AR+H^+       ",&
   "AR^++H^-=>AR+H       ","AR^++H2=>ARH^++H     ","AR^++OH^-=>AR+OH     ","AR^++OH^-=>AR+O+H    ","AR^++H2O=>AR+H2O^+   ",&
   "AR^++H2O=>ARH^++OH   ","H2O^++AR=>AR^++H2O   ","H3O^++AR=>ARH^++H2O  ","ARH^++O^-=>AR+O+H    ","bolsig:AR->AR        ",&
   "bolsig:AR->AR^+      ","bolsig:AR->AR*       ","bolsig:AR*->AR       ","bolsig:AR*->AR^+     ","AR2^++AR=>AR^++AR+AR ",&
   "AR*+AR*=>AR2^++E     ","AR*+AR+AR=>AR+AR+AR  ","AR^+=>AR(W)^+        ","AR2^+=>AR2(W)^+      ","H2O^+=>H2O(W)^+      ",&
   "OH^+=>OH(W)^+        ","H^-=>H(W)^-          ","E=>E(W)              "/
  data bolsig_species(1:bolsig_species_max) &
  /"AR ","H2O","AR*"/
contains
!-----------------------------------------------------------------------------------------------------------------------------------
!
! initialization
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_init()
  implicit none
  character(256) :: string
  integer :: i, j, k
  write(*,"(/,A)") "ZDPlasKin (version " // "2.0a" // ") INIT:"
  if( lZDPlasKin_init ) call ZDPlasKin_stop("   ERROR: the ZDPlasKin library has been initialized")
  write(string,*) species_max
  write(*,"(2x,A)")  "species        ... " // trim(adjustl(string))
  write(string,*) reactions_max
  write(*,"(2x,A)")  "reactions      ... " // trim(adjustl(string))
  if(species_max<=0 .or. reactions_max<=0) call ZDPlasKin_stop("   ERROR: wrong preconfig data")
  write(*,"(2x,A,$)") "BOLSIG+ loader ... " // trim(adjustl(bolsigfile)) // " : "
  call ZDPlasKin_bolsig_Init(bolsigfile)
  do i = 1, bolsig_species_max
    call ZDPlasKin_bolsig_ReadCollisions(trim(bolsig_species(i)))
    j = bolsig_collisions_max
    call ZDPlasKin_bolsig_GetCollisions(k,bolsig_collisions_max)
    if(bolsig_collisions_max <= j) then
      write(*,*)
      call ZDPlasKin_stop("ERROR: wrong file or missing data " // &
                        "(" // trim(adjustl(bolsigfile)) // ": <" // trim(bolsig_species(i)) // ">).")
    endif
  enddo
  if(bolsig_species_max /= k) then
    write(*,*)
    call ZDPlasKin_stop("ERROR: internal error in BOLSIG+ loader")
  endif
  write(string,*) bolsig_species_max
  write(*,"(A,$)") trim(adjustl(string)) // " species & "
  write(string,*) bolsig_collisions_max
  write(*,"(A)")   trim(adjustl(string)) // " collisions"
  write(*,"(2x,A,$)") "species  link  ... "
  j = 0
  do i = 1, bolsig_species_max
    j = j + 1
    k = 1
    do while(k<=species_max .and. bolsig_species_index(i)<=0)
      call ZDPlasKin_bolsig_GetSpeciesName(string,i)
      if(trim(species_name(k)) == trim(string)) then
        bolsig_species_index(i) = k
      else
        k = k + 1
      endif
    enddo
    if(bolsig_species_index(i) <= 0) call ZDPlasKin_stop("cannot find species link for <" // trim(string) // ">")
  enddo
  write(string,*) j
  write(*,"(A)") trim(adjustl(string))
  write(*,"(2x,A,$)") "process  link  ... "
  i = 1
  j = 1
  do while(i<=reactions_max .and. j<=bolsig_rates_max)
    if(reaction_sign(i)(1:7) == "bolsig:") then
      k = 1
      do while(k<=bolsig_collisions_max .and. bolsig_pointer(j)<=0)
        call ZDPlasKin_bolsig_GetReactionName(string,k)
        if(trim(string) == trim(reaction_sign(i)(8:))) then
          bolsig_pointer(j) = k
        else
          k = k + 1
        endif
      enddo
      if(bolsig_pointer(j) <= 0) call ZDPlasKin_stop("cannot find processes link for <" // trim(reaction_sign(i)) // ">")
      j = j + 1
    endif
    i = i + 1
  enddo
  if(j <= bolsig_rates_max) then
    call ZDPlasKin_stop("internal error")
  else
    write(string,*) bolsig_rates_max
    write(*,"(A)") trim(adjustl(string))
  endif
  i = 0
  do while((1.0d0+10.0d0**(i-1)) /= 1.0d0)
    i = i - 1
  enddo
  mach_accur = 10.0d0**i
  mach_tiny  = sqrt( tiny(mach_tiny) )
  lZDPlasKin_init = .true.
  call ZDPlasKin_reset()
  write(*,"(A,/)") "ZDPlasKin INIT DONE"
  return
end subroutine ZDPlasKin_init
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using implicit solver dvode_f90
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep(time,dtime)
  use dvode_f90_m, only : dvode_f90
  implicit none
  double precision, intent(in)    ::  time
  double precision, intent(inout) :: dtime
  double precision, save :: densav(vode_neq) = 0.0d0, cfgsav(3) = 0.0d0
  double precision :: tout
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(time < tsav) vode_istate = 1
  tsav = time
  if(dtime > 0.0d0) then
    vode_itask = 1
    tout = time + dtime
    if(dtime < mach_accur*abs(tout)) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: dtime parameter is too small (subroutine ZDPlasKin_timestep)")
  else
    vode_itask = 2
    tout = ( 1.0d0 + mach_accur ) * time + mach_tiny
  endif
  dens_loc(1:species_max,0) = density(:)
  dens_loc(1:species_max,1) = 0.5d0 * ( density(:) + abs( density(:) ) )
  if(any(dens_loc(1:species_max,1) /= densav(1:species_max))) vode_istate = 1
  densav(1:species_max) = dens_loc(1:species_max,1)
  if(vode_istate /= 1 .and. any( abs(cfgsav(:)-ZDPlasKin_cfg(1:3)) > cfg_rtol*abs(cfgsav(:)+ZDPlasKin_cfg(1:3)) )) vode_istate = 1
  cfgsav(:) = ZDPlasKin_cfg(1:3)
  if( lgas_heating ) then
    if(ZDPlasKin_cfg(1) /= densav(species_max+1)) vode_istate = 1
    densav(species_max+1) = ZDPlasKin_cfg(1)
  endif
  call dvode_f90(ZDPlasKin_fex,vode_neq,densav,tsav,tout,vode_itask,vode_istate,vode_options,j_fcn=ZDPlasKin_jex)
  if(vode_istate < 0) then
    write(*,"(A,1pd11.4)") "Tgas   =", ZDPlasKin_cfg(1)
    write(*,"(A,1pd11.4)") "    EN =", ZDPlasKin_cfg(3)
    write(*,"(A,1pd11.4)") "    Te =", ZDPlasKin_cfg(4)
    call ZDPlasKin_stop("ZDPlasKin ERROR: DVODE solver issued an error (subroutine ZDPlasKin_timestep)")
  endif
  if( lgas_heating ) ZDPlasKin_cfg(1) = densav(species_max+1)
  density(:) = dens_loc(1:species_max,0) - dens_loc(1:species_max,1) + densav(1:species_max)
  if(dtime <= 0.0d0) dtime = tsav - time
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep
!-----------------------------------------------------------------------------------------------------------------------------------
!
! timestep integration using explicit Euler method
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_timestep_explicit(time,dtime,rtol_loc,atol_loc,switch_implicit)
  implicit none
  double precision, intent(in) ::  time, rtol_loc, atol_loc
  double precision, intent(inout) :: dtime
  double precision, optional, intent(in) :: switch_implicit
  double precision :: time_loc, time_end, dtime_loc, dtime_max
  logical, save :: lwarn = .true.
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( lqtplaskin .and. lqtplaskin_first ) call ZDPlasKin_write_qtplaskin(time)
  if(rtol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: rtol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  if(atol_loc <= 0.0d0) call ZDPlasKin_stop("ZDPlasKin ERROR: atol_loc must be positive (subroutine ZDPlasKin_timestep_explicit)")
  tsav     = time
  time_loc = 0.0d0
  time_end = 0.5d0 * ( dtime + abs(dtime) ) + mach_tiny
  do while(time_loc < time_end)
    dens_loc(1:species_max,0) = density(:)
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    dens_loc(:,1) = 0.5d0 * ( dens_loc(:,0) + abs( dens_loc(:,0) ) )
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,1),dens_loc(:,2))
    where(dens_loc(:,2) >= 0.0d0)
      dens_loc(:,3) = + dens_loc(:,2) /    ( rtol_loc * dens_loc(:,1) + atol_loc ) 
    elsewhere
      dens_loc(:,3) = - dens_loc(:,2) / min( rtol_loc * dens_loc(:,1) + atol_loc , dens_loc(:,1) + mach_tiny )
    endwhere
    dtime_loc = 1.0d0 / ( maxval( dens_loc(:,3) ) + mach_tiny )
    if(dtime > 0.0d0) then
      dtime_max = dtime - time_loc
      dtime_loc = min( dtime_loc , dtime_max )
      if( present(switch_implicit) ) then
        if(dtime_loc*switch_implicit < dtime_max) then
          if(lprint .and. lwarn) then
            write(*,"(A,/,A,1pd9.2,A)") "ZDPlasKin INFO: low efficiency of Euler method (subroutine ZDPlasKin_timestep_explicit)", &
                        "                ZDPlasKin_timestep subroutine will be used in similar conditions (", switch_implicit, ")"
            lwarn = .false.
          endif
          time_loc = tsav
          density(:) = dens_loc(1:species_max,0)
          if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0)
          call ZDPlasKin_timestep(time_loc,dtime_max)
          return
        endif
      endif
    else
      dtime = dtime_loc
    endif
    time_loc = time_loc + dtime_loc
    tsav     = time     +  time_loc
    density(:) = dens_loc(1:species_max,0) + dtime_loc * dens_loc(1:species_max,2)
    if( lgas_heating ) ZDPlasKin_cfg(1) = dens_loc(species_max+1,0) + dtime_loc * dens_loc(species_max+1,2)
  enddo
  if( lstat_accum ) then
    call ZDPlasKin_get_rates(SOURCE_TERMS=dens_loc(1:species_max,0),REACTION_RATES=rrt_loc)
    stat_dens(:) = stat_dens(:) + dtime * density(:)
    stat_src(:)  = stat_src(:)  + dtime * dens_loc(1:species_max,0)
    stat_rrt(:)  = stat_rrt(:)  + dtime * rrt_loc(:)
    stat_time    = stat_time    + dtime
  endif
  if( lqtplaskin ) call ZDPlasKin_write_qtplaskin(time+dtime)
  return
end subroutine ZDPlasKin_timestep_explicit
!-----------------------------------------------------------------------------------------------------------------------------------
!
! update BOLSIG+ solution and get electron parameters
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_bolsig_rates(lbolsig_force)
  implicit none
  logical, optional, intent(in) :: lbolsig_force
  logical :: lforce
  integer :: i, j, k
  integer, save :: bolsig_points_max
  logical, save :: lfirst = .true., leecol = .true.
  double precision :: error, density_loc, cfg_loc(6+bolsig_species_max)
  double precision, save :: low_density_limit = bolsig_rtol, bolsig_mesh_a, bolsig_mesh_b
  double precision, save, allocatable :: bolsig_cfg(:,:), bolsig_reslt(:,:)
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    bolsig_mesh_a = 1.0d0 / log( 1.0d0 + bolsig_rtol )
    bolsig_mesh_b = bolsig_mesh_a * log( bolsig_field_min ) - 0.5d0
    bolsig_points_max = int( bolsig_mesh_a * log( bolsig_field_max ) - bolsig_mesh_b )
    allocate(bolsig_rates(bolsig_collisions_max), &
             bolsig_cfg(6+bolsig_species_max,0:bolsig_points_max), &
             bolsig_reslt(10+bolsig_collisions_max,0:bolsig_points_max),stat=i)
    if(i /= 0) call ZDPlasKin_stop("ZDPlasKin ERROR: memory allocation error (subroutine ZDPlasKin_bolsig_rates)")
    bolsig_cfg(:,:) = 0.0d0
    lfirst = .false.
  endif
  if( present(lbolsig_force) ) then
    lforce = lbolsig_force
  else
    lforce = .false.
  endif
  if(ZDPlasKin_cfg(1) <= 0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
  if(.not. lbolsig_Maxwell_EEDF ) then
    if(ZDPlasKin_cfg(2) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_bolsig_rates)")
    if(ZDPlasKin_cfg(3) < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_bolsig_rates)")
    ZDPlasKin_cfg(4) = 0.0d0
  else
    if(ZDPlasKin_cfg(4) <= 0.0d0) then
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELECTRON_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)")
    elseif(lprint .and. ZDPlasKin_cfg(4) < ZDPlasKin_cfg(1)) then
      write(*,"(A)") "ZDPlasKin INFO: ELECTRON_TEMPERATURE is below GAS_TEMPERATURE (subroutine ZDPlasKin_bolsig_rates)"
    endif
    ZDPlasKin_cfg(2:3) = 0.0d0
  endif
  density_loc = 0.5d0 * ( sum(density(bolsig_species_index(:))) + sum(abs(density(bolsig_species_index(:)))) )
  if(density_loc <= mach_tiny) then
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined densities configured for BOLSIG+ solver " // &
                                                                     "(subroutine ZDPlasKin_bolsig_rates)")
  elseif( lprint ) then
    call ZDPlasKin_get_density_total(ALL_NEUTRAL=error)
    error = abs( 1.0d0 - density_loc / error )
    if(error > low_density_limit) then
      write(*,"(A,1pd9.2)") "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds", error
      low_density_limit = sqrt( low_density_limit )
    endif
  endif
  cfg_loc(1:4) = ZDPlasKin_cfg(1:4)
  cfg_loc(2)   = cfg_loc(2) * 1.0d-6
  cfg_loc(6)   = 0.5d0 * ( density(species_electrons) + abs(density(species_electrons)) )
  cfg_loc(5)   = cfg_loc(6) * 1.0d6
  cfg_loc(6)   = cfg_loc(6) / density_loc
  if(cfg_loc(6) < bolsig_eecol_frac) then
    cfg_loc(6) = 0.0d0
  elseif(lprint .and. leecol) then
    write(*,"(A)") "ZDPlasKin INFO: set electron-electron collisions ON ..."
    leecol = .false.
  endif
  cfg_loc(7:) = 0.5d0 * ( density(bolsig_species_index(:)) + abs(density(bolsig_species_index(:))) ) / density_loc
  if(lbolsig_Maxwell_EEDF .or. bolsig_points_max==0) then
    lforce = .true.
    i = 0
  else
    error = min( max( cfg_loc(3) , bolsig_field_min ) , bolsig_field_max )
    i = int( bolsig_mesh_a * log( error ) - bolsig_mesh_b )
    i = max(0,min(i,bolsig_points_max))
  endif
  if( lforce ) then
    error = 2.0d0
  else
    if(lbolsig_ignore_gas_temp .and. bolsig_cfg(1,i)>0.0d0) then
      cfg_loc(1) = bolsig_cfg(1,i)
      error = 0.0d0
    else
      error = abs( ( cfg_loc(1) - bolsig_cfg(1,i) ) / ( 0.5d0 * ( cfg_loc(1) + bolsig_cfg(1,i) ) + mach_tiny) ) / bolsig_rtol_half
    endif
    if(error <= 1.0d0) then
      error = abs( ( cfg_loc(2) - bolsig_cfg(2,i) ) / ( 0.5d0 * ( cfg_loc(2) + bolsig_cfg(2,i) ) + mach_tiny) ) / bolsig_rtol_half
      if(error <= 1.0d0) then
        error = abs( ( cfg_loc(3) - bolsig_cfg(3,i) ) / ( 0.5d0 * ( cfg_loc(3) + bolsig_cfg(3,i) ) + mach_tiny) ) / bolsig_rtol
        if(error <= 1.0d0) then
          error = abs( ( max(cfg_loc(6),bolsig_eecol_frac) - max(bolsig_cfg(6,i),bolsig_eecol_frac) ) &
           / ( 0.5d0 * ( max(cfg_loc(6),bolsig_eecol_frac) + max(bolsig_cfg(6,i),bolsig_eecol_frac) ) + mach_tiny) ) &
           / bolsig_rtol_half
          if(error <= 1.0d0) error = maxval( abs( cfg_loc(7:) - bolsig_cfg(7:,i) ) ) &
                                     / ( 0.5d0 * maxval( cfg_loc(7:) + bolsig_cfg(7:,i) ) ) / bolsig_rtol
        endif
      endif
    endif
  endif
  if(error > 1.0d0) then
    j = 6 + bolsig_species_max
    k = 10 + bolsig_collisions_max
    bolsig_cfg(:,i) = cfg_loc(:)
    call ZDPlasKin_bolsig_SolveBoltzmann(j,bolsig_cfg(1:j,i),k,bolsig_reslt(1:k,i))
    if(.not. lbolsig_Maxwell_EEDF) then
      bolsig_reslt(2,i) = bolsig_reslt(2, i) * eV_to_K / 1.5d0
    else
      bolsig_reslt(2,i) = cfg_loc(4)
    endif
    bolsig_reslt(3, i) = bolsig_reslt(3, i) * 1.0d-2
    bolsig_reslt(4, i) = bolsig_reslt(4, i) * 1.0d-2 / density_loc
    bolsig_reslt(5, i) = bolsig_reslt(5, i) * 1.0d-2
    bolsig_reslt(6, i) = bolsig_reslt(6, i) * 1.0d-2
    bolsig_reslt(7:,i) = bolsig_reslt(7:,i) * 1.0d6
  endif
  ZDPlasKin_cfg(3:12) = bolsig_reslt(:10,i)
  bolsig_rates(:)     = bolsig_reslt(11:,i)
  return
end subroutine ZDPlasKin_bolsig_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get index of species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_species_index(str,i)
  implicit none
  character(*), intent(in ) :: str
  integer,      intent(out) :: i
  character(species_length) :: string
  integer :: j, istr
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  string = trim(adjustl(str))
  istr   = len_trim(string)
  do i = 1, istr
    j = iachar(string(i:i))
    if(j>=97 .and. j<=122) string(i:i) = achar(j-32)
  enddo
  i = 0
  j = 0
  do while(i==0 .and. j<species_max)
    j = j + 1
    if(string(1:istr) == trim(species_name(j))) i = j
  enddo
  if(i <= 0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: cannot identify species <"//trim(str)//"> (subroutine ZDPlasKin_get_species_index)")
  return
end subroutine ZDPlasKin_get_species_index
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get density for species
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in) :: string
  logical, optional, intent(in) :: LDENS_CONST
  double precision, optional, intent(in) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present(DENS) ) density(i) = DENS
  if( present(LDENS_CONST) ) then
    density_constant(i) = LDENS_CONST
    ldensity_constant   = any( density_constant(:) )
  endif
  return
end subroutine ZDPlasKin_set_density
subroutine ZDPlasKin_get_density(string,DENS,LDENS_CONST)
  implicit none
  character(*), intent(in ) :: string
  logical, optional, intent(out) :: LDENS_CONST
  double precision, optional, intent(out) :: DENS
  integer :: i
  call ZDPlasKin_get_species_index(string,i)
  if( present( DENS)       )  DENS       = density(i)
  if( present(LDENS_CONST) ) LDENS_CONST = density_constant(i)
  return
end subroutine ZDPlasKin_get_density
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get total densities
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_density_total(ALL_SPECIES,ALL_NEUTRAL,ALL_ION_POSITIVE,ALL_ION_NEGATIVE,ALL_CHARGE)
  double precision, optional, intent(out) :: ALL_SPECIES, ALL_NEUTRAL, ALL_ION_POSITIVE, ALL_ION_NEGATIVE, ALL_CHARGE
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ALL_SPECIES)      ) ALL_SPECIES      = sum(density(:))
  if( present(ALL_NEUTRAL)      ) ALL_NEUTRAL      = sum(density(:), mask = species_charge(:)==0)
  if( present(ALL_ION_POSITIVE) ) ALL_ION_POSITIVE = sum(density(:), mask = species_charge(:)>0)
  if( present(ALL_ION_NEGATIVE) ) ALL_ION_NEGATIVE = sum(density(:), mask = species_charge(:)<0) - density(species_electrons)
  if( present(ALL_CHARGE)       ) ALL_CHARGE       = sum(density(:) * dble(species_charge(:)))
  return
end subroutine ZDPlasKin_get_density_total
!-----------------------------------------------------------------------------------------------------------------------------------
!
! get species source terms & reaction rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_get_rates(SOURCE_TERMS,REACTION_RATES,SOURCE_TERMS_MATRIX,MEAN_DENSITY, &
                               MEAN_SOURCE_TERMS,MEAN_REACTION_RATES,MEAN_SOURCE_TERMS_MATRIX)
  double precision, optional, intent(out) :: SOURCE_TERMS(species_max), REACTION_RATES(reactions_max), &
                                             SOURCE_TERMS_MATRIX(species_max,reactions_max), MEAN_DENSITY(species_max), &
                                             MEAN_SOURCE_TERMS(species_max), MEAN_REACTION_RATES(reactions_max), &
                                             MEAN_SOURCE_TERMS_MATRIX(species_max,reactions_max)
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if(present(SOURCE_TERMS) .or. present(REACTION_RATES) .or. present(SOURCE_TERMS_MATRIX)) then
    dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
    if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
    call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
    if( present(SOURCE_TERMS)                 ) SOURCE_TERMS(:)   = dens_loc(1:species_max,1)
    if( present(REACTION_RATES)               ) REACTION_RATES(:) = rrt(:)
    if( present(SOURCE_TERMS_MATRIX)          ) call ZDPlasKin_reac_source_matrix(rrt(:),SOURCE_TERMS_MATRIX(:,:))
  endif
  if(present(MEAN_DENSITY)        .or. present(MEAN_SOURCE_TERMS) .or. &
     present(MEAN_REACTION_RATES) .or. present(MEAN_SOURCE_TERMS_MATRIX)) then
    if( lstat_accum ) then
      if(stat_time > 0.0d0) then
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = stat_dens(:) / stat_time
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = stat_src(:)  / stat_time
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = stat_rrt(:)  / stat_time
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) then
          call ZDPlasKin_reac_source_matrix(stat_rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
          MEAN_SOURCE_TERMS_MATRIX(:,:)  =  MEAN_SOURCE_TERMS_MATRIX(:,:) / stat_time
        endif
      else
        dens_loc(1:species_max,0) = 0.5d0 * ( density(:) + abs( density(:) ) )
        if( lgas_heating ) dens_loc(species_max+1,0) = ZDPlasKin_cfg(1)
        call ZDPlasKin_fex(vode_neq,tsav,dens_loc(:,0),dens_loc(:,1))
        if( present(MEAN_DENSITY)             ) MEAN_DENSITY        = density(:)
        if( present(MEAN_SOURCE_TERMS)        ) MEAN_SOURCE_TERMS   = dens_loc(1:species_max,1)
        if( present(MEAN_REACTION_RATES)      ) MEAN_REACTION_RATES = rrt(:)
        if( present(MEAN_SOURCE_TERMS_MATRIX) ) call ZDPlasKin_reac_source_matrix(rrt(:),MEAN_SOURCE_TERMS_MATRIX(:,:))
      endif
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_get_rates)")
    endif
  endif
  return
end subroutine ZDPlasKin_get_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set config
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_config(ATOL,RTOL,SILENCE_MODE,STAT_ACCUM,QTPLASKIN_SAVE,BOLSIG_EE_FRAC,BOLSIG_IGNORE_GAS_TEMPERATURE)
  use dvode_f90_m, only : set_intermediate_opts
  implicit none
  logical, optional, intent(in) :: SILENCE_MODE, STAT_ACCUM, QTPLASKIN_SAVE, BOLSIG_IGNORE_GAS_TEMPERATURE
  double precision, optional, intent(in) :: ATOL, RTOL, BOLSIG_EE_FRAC
  integer :: i
  logical, save :: lfirst = .true.
  integer, save :: bounded_components(vode_neq)
  double precision :: atol_loc, rtol_loc
  double precision, save :: atol_save = -1.0d0, rtol_save = -1.0d0
  if( lfirst ) then
    if(.not. lZDPlasKin_init) call ZDPlasKin_init()
    do i = 1, vode_neq
      bounded_components(i) = i
    enddo
    lfirst = .false.
  endif
  if( present(SILENCE_MODE) ) lprint = ( .not. SILENCE_MODE )
  if( present(BOLSIG_EE_FRAC) ) bolsig_eecol_frac = 0.5d0 * ( BOLSIG_EE_FRAC + abs(BOLSIG_EE_FRAC) )
  if( present(BOLSIG_IGNORE_GAS_TEMPERATURE) ) lbolsig_ignore_gas_temp = BOLSIG_IGNORE_GAS_TEMPERATURE
  if( present(STAT_ACCUM) ) then
    if( lprint ) then
      if(lstat_accum .neqv. STAT_ACCUM) then
        if( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set statistic acquisition OFF ..."
        endif
      elseif( STAT_ACCUM ) then
          write(*,"(A)") "ZDPlasKin INFO: reset statistic acquisition data ..."
      endif
    endif
    stat_dens(:) = 0.0d0
    stat_src(:)  = 0.0d0
    stat_rrt(:)  = 0.0d0
    stat_time    = 0.0d0
    lstat_accum  = STAT_ACCUM
  endif
  if( present(QTPLASKIN_SAVE) ) then
    if( lprint ) then
      if(lqtplaskin .neqv. QTPLASKIN_SAVE) then
        if( QTPLASKIN_SAVE ) then
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format ON ..."
        else
          write(*,"(A)") "ZDPlasKin INFO: set autosave in QTplaskin format OFF ..."
        endif
      endif
    endif
    lqtplaskin = QTPLASKIN_SAVE
  endif
  if( present(ATOL) ) then
    atol_loc = ATOL
  else
    atol_loc = atol_save
  endif
  if( present(RTOL) ) then
    rtol_loc = RTOL
  else
    rtol_loc = rtol_save
  endif
  if(min(atol_loc,rtol_loc)<0.0d0 .or. max(atol_loc,rtol_loc)<=0.0d0) &
    call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ATOL/RTOL (ZDPlasKin_set_config)")
  if(atol_loc/=atol_save .or. rtol_loc/=rtol_save) then
    atol_save = atol_loc
    rtol_save = rtol_loc
    if( lprint ) write(*,"(2(A,1pd9.2),A)") "ZDPlasKin INFO: set accuracy", atol_save, " (absolute) &", rtol_save, " (relative)"
    dens_loc(:,0) = 0.0d0
    dens_loc(:,1) = huge(dens_loc)
    vode_options  = set_intermediate_opts(abserr=atol_save,relerr=rtol_save, &
                                          dense_j=.true.,user_supplied_jacobian=.true., &
                                          constrained=bounded_components(:),clower=dens_loc(:,0),cupper=dens_loc(:,1))
    if(vode_istate /= 1) vode_istate = 3
  endif
  return
end subroutine ZDPlasKin_set_config
!-----------------------------------------------------------------------------------------------------------------------------------
!
! set/get conditions
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_set_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,GAS_HEATING,SPEC_HEAT_RATIO,HEAT_SOURCE,SOFT_RESET)
  implicit none
  double precision, optional, intent(in) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, ELEC_TEMPERATURE, &
                                            SPEC_HEAT_RATIO, HEAT_SOURCE
  logical,          optional, intent(in) :: GAS_HEATING, SOFT_RESET
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(GAS_TEMPERATURE) ) then
    if(GAS_TEMPERATURE <= 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined GAS_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(1) = GAS_TEMPERATURE
  endif
  if( present(REDUCED_FREQUENCY) ) then
    if(REDUCED_FREQUENCY < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FREQUENCY (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(2) = REDUCED_FREQUENCY
  endif
  if( present(REDUCED_FIELD) ) then
    if(REDUCED_FIELD < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined REDUCED_FIELD (subroutine ZDPlasKin_set_conditions)")
    ZDPlasKin_cfg(3) = REDUCED_FIELD
  endif
  if( present(SOFT_RESET) ) then
    if( SOFT_RESET ) vode_istate = 1
  endif
  if( present(ELEC_TEMPERATURE) ) then
    if(ELEC_TEMPERATURE < 0.0d0) &
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong or undefined ELEC_TEMPERATURE (subroutine ZDPlasKin_set_conditions)")
    if(ELEC_TEMPERATURE > 0.0d0) then
      lbolsig_Maxwell_EEDF = .true.
    else
      lbolsig_Maxwell_EEDF = .false.
    endif
    ZDPlasKin_cfg(4) = ELEC_TEMPERATURE
  endif
  if( present(GAS_HEATING) ) then
    if(lgas_heating .neqv. GAS_HEATING) then
      if( GAS_HEATING ) then
        if(present(SPEC_HEAT_RATIO) .or. ZDPlasKin_cfg(13)>0.0d0) then
          if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating ON ..."
        else
          ZDPlasKin_cfg(13) = 2.0d0/3.0d0
          if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set gas heating ON; specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
        endif
      else
        if( lprint ) write(*,"(A)") "ZDPlasKin INFO: set gas heating OFF ..."
      endif
      lgas_heating = GAS_HEATING
    endif
  endif
  if( present(SPEC_HEAT_RATIO) ) then
    if(SPEC_HEAT_RATIO > 1.0d0) then
      ZDPlasKin_cfg(13) = SPEC_HEAT_RATIO - 1.0d0
      if( lprint ) write(*,"(A,1pd9.2)") "ZDPlasKin INFO: set specific heat ratio =", ZDPlasKin_cfg(13) + 1.0d0
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: wrong value of SPEC_HEAT_RATIO (subroutine ZDPlasKin_set_conditions)")
    endif
  endif
  if( present(HEAT_SOURCE) ) then
    ZDPlasKin_cfg(14) = HEAT_SOURCE
    if( lprint ) write(*,"(A,1pd9.2,A)") "ZDPlasKin INFO: set heat source =", ZDPlasKin_cfg(14), " W/cm3"
  endif
end subroutine ZDPlasKin_set_conditions
subroutine ZDPlasKin_get_conditions(GAS_TEMPERATURE,REDUCED_FREQUENCY,REDUCED_FIELD, &
                                    ELEC_TEMPERATURE,ELEC_DRIFT_VELOCITY,ELEC_DIFF_COEFF,ELEC_MOBILITY_N, &
                                    ELEC_MU_EPS_N,ELEC_DIFF_EPS_N,ELEC_FREQUENCY_N, &
                                    ELEC_POWER_N,ELEC_POWER_ELASTIC_N,ELEC_POWER_INELASTIC_N,ELEC_EEDF)
  implicit none
  double precision, optional, intent(out) :: GAS_TEMPERATURE, REDUCED_FREQUENCY, REDUCED_FIELD, &
                                             ELEC_TEMPERATURE, ELEC_DRIFT_VELOCITY, ELEC_DIFF_COEFF, ELEC_MOBILITY_N, &
                                             ELEC_MU_EPS_N, ELEC_DIFF_EPS_N, ELEC_FREQUENCY_N, &
                                             ELEC_POWER_N, ELEC_POWER_ELASTIC_N, ELEC_POWER_INELASTIC_N
  double precision, optional, dimension(:,:), intent(out) :: ELEC_EEDF
  integer :: i
  double precision :: x,y
  if(.not. lZDPlasKin_init) call ZDPlasKin_init()
  if( present(ELEC_EEDF) ) then
    call ZDPlasKin_bolsig_rates(lbolsig_force=.true.)
  else
    call ZDPlasKin_bolsig_rates()
  endif
  if( present(GAS_TEMPERATURE)        ) GAS_TEMPERATURE        = ZDPlasKin_cfg(1)
  if( present(REDUCED_FREQUENCY)      ) REDUCED_FREQUENCY      = ZDPlasKin_cfg(2)
  if( present(REDUCED_FIELD)          ) REDUCED_FIELD          = ZDPlasKin_cfg(3)
  if( present(ELEC_TEMPERATURE)       ) ELEC_TEMPERATURE       = ZDPlasKin_cfg(4)
  if( present(ELEC_DRIFT_VELOCITY)    ) ELEC_DRIFT_VELOCITY    = ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
  if( present(ELEC_DIFF_COEFF)        ) ELEC_DIFF_COEFF        = ZDPlasKin_cfg(6)
  if( present(ELEC_MOBILITY_N)        ) ELEC_MOBILITY_N        = ZDPlasKin_cfg(5)
  if( present(ELEC_MU_EPS_N)          ) ELEC_MU_EPS_N          = ZDPlasKin_cfg(7)
  if( present(ELEC_DIFF_EPS_N)        ) ELEC_DIFF_EPS_N        = ZDPlasKin_cfg(8)
  if( present(ELEC_FREQUENCY_N)       ) ELEC_FREQUENCY_N       = ZDPlasKin_cfg(9)
  if( present(ELEC_POWER_N)           ) ELEC_POWER_N           = ZDPlasKin_cfg(10)
  if( present(ELEC_POWER_ELASTIC_N)   ) ELEC_POWER_ELASTIC_N   = ZDPlasKin_cfg(11)
  if( present(ELEC_POWER_INELASTIC_N) ) ELEC_POWER_INELASTIC_N = ZDPlasKin_cfg(12)
  if( present(ELEC_EEDF) ) then
    ELEC_EEDF = 0d0
  	 if( size(ELEC_EEDF,dim=1) < 2 ) then
      if(lprint) write(*,"(A)") &
  	     "ZDPlasKin WARNING: insufficient first dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
  	  else
  		y = 1.0d0
  		do i = 1, size(ELEC_EEDF,dim=2)
  		  call ZDPlasKin_bolsig_GetEEDF(i,x,y)
  		  if( x >= 0d0 .and. y > 0d0) then
  			ELEC_EEDF(1,i) = x
  			ELEC_EEDF(2,i) = y
  		  else
  			exit
  		  endif
  		enddo
  		if(lprint .and. y>0d0) write(*,"(A)") &
  		  "ZDPlasKin WARNING: insufficient second dimention of array ELEC_EEDF (subroutine ZDPlasKin_get_conditions)"
     endif
  endif
end subroutine ZDPlasKin_get_conditions
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reset
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reset()
  implicit none
  vode_istate         =  1
  density(:)          =  0.0d0
  ZDPlasKin_cfg(:)    =  0.0d0
  ldensity_constant   = .false.
  density_constant(:) = .false.
  lreaction_block(:)  = .false.
  lprint              = .true.
  lstat_accum         = .false.
  lqtplaskin          = .false.
  lgas_heating        = .false.
  bolsig_eecol_frac       = bolsig_eecol_frac_def
  lbolsig_ignore_gas_temp = .false.
  lbolsig_Maxwell_EEDF    = .false.
  write(*,"(A)") "ZDPlasKin INFO: reset data and configuration"
  call ZDPlasKin_set_config(ATOL=vode_atol,RTOL=vode_rtol)
  return
end subroutine ZDPlasKin_reset
!-----------------------------------------------------------------------------------------------------------------------------------
!
! stop
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_stop(string)
  implicit none
  character(*), intent(in) :: string
  if(string /= "") write(*,"(A)") trim(string)
  write(*,"(A,$)") "PRESS ENTER TO EXIT ... "
  read(*,*)
  stop
end subroutine ZDPlasKin_stop
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data to file
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_file(FILE_SPECIES,FILE_REACTIONS,FILE_SOURCE_MATRIX,FILE_UNIT)
  implicit none
  character(*), optional, intent(in) :: FILE_SPECIES, FILE_REACTIONS, FILE_SOURCE_MATRIX
  integer, optional, intent(in) :: FILE_UNIT
  logical :: lerror
  integer :: i
  if( present(FILE_UNIT) ) ifile_unit = FILE_UNIT
  if( present(FILE_SPECIES) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SPECIES)),action="write",err=100)
    do i = 1, species_max
      write(ifile_unit,111,err=100) i, species_name(i)
    enddo
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SPECIES)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
111 format(i2,1x,A9)
  endif
  if( present(FILE_REACTIONS) ) then
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_REACTIONS)),action="write",err=200)
    do i = 1, reactions_max
      write(ifile_unit,211,err=200) i, reaction_sign(i)
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_REACTIONS)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
211 format(i3,1x,A21)
  endif
  if( present(FILE_SOURCE_MATRIX) ) then
    if( lstat_accum ) then
      call ZDPlasKin_reac_source_matrix(stat_rrt(:),mrtm(:,:))
      if(stat_time > 0.0d0) mrtm(:,:) = mrtm(:,:) / stat_time
    else
      call ZDPlasKin_stop("ZDPlasKin ERROR: set statistics acquisition ON before (subroutine ZDPlasKin_write_file)")
    endif
    lerror = .true.
    open(ifile_unit,file=trim(adjustl(FILE_SOURCE_MATRIX)),action="write",err=300)
    write(ifile_unit,311,err=300) ( i, i = 1, species_max )
    write(ifile_unit,312,err=300) "N", "reaction", ( trim(species_name(i)), i = 1, species_max )
    do i = 1, reactions_max
      write(ifile_unit,313,err=300) i, reaction_sign(i), mrtm(:,i)
    enddo
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file <" &
                                    // trim(adjustl(FILE_SOURCE_MATRIX)) // "> (subroutine ZDPlasKin_write_file)")
    close(ifile_unit)
311 format(261x,47(1x,i9))
312 format(A3,1x,A21,1x,47(1x,A9))
313 format(i3,1x,A21,1x,47(1x,1pd9.2))
  endif
  return
end subroutine ZDPlasKin_write_file
!-----------------------------------------------------------------------------------------------------------------------------------
!
! save data in qtplaskin format
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_write_qtplaskin(time,LFORCE_WRITE)
  implicit none
  double precision, intent(in) :: time
  logical, optional, intent(in) :: LFORCE_WRITE
  integer, parameter :: idef_data = 5
  character(24), parameter :: qtplaskin_names(idef_data) = (/ "Reduced field [Td]      ", "Gas temperature [K]     ", &
                                  "Electron temperature [K]", "Current density [A/cm2] ", "Power density [W/cm3]   " /)
  double precision, save :: densav(0:species_max,2) = -huge(densav)
  double precision :: rtol, cond(idef_data)
  logical, save :: lfirst = .true.
  logical :: lerror
  integer, save :: iuser_data = 0
  integer :: i
  if( time < densav(0,1) ) lfirst = .true.
  if( lfirst ) then
    call ZDPlasKin_write_file(FILE_SPECIES="qt_species_list.txt",FILE_REACTIONS="qt_reactions_list.txt")
    if( allocated(qtplaskin_user_data) ) then
      iuser_data = size(qtplaskin_user_data)
      iuser_data = min(iuser_data,90)
      if( iuser_data > 0 ) then
        if( allocated(qtplaskin_user_names) ) then
          if( size(qtplaskin_user_names) /= iuser_data ) deallocate(qtplaskin_user_names)
        endif
        if( .not. allocated(qtplaskin_user_names) ) then
          allocate(qtplaskin_user_names(iuser_data))
          do i = 1, iuser_data
            write(qtplaskin_user_names(i),"(A,i2.2)") "user defined #", i
          enddo
        endif
      endif
    endif
    lerror = .true.
    open(ifile_unit,file="qt_conditions_list.txt",action="write",err=100)
    do i = 1, idef_data
      write(ifile_unit,"(i3,1x,A)",err=100) i, trim(adjustl(qtplaskin_names(i)))
    enddo
    if( iuser_data > 0 ) then
      do i = 1, iuser_data
        write(ifile_unit,"(i3,1x,A)",err=100) (i+idef_data), trim(adjustl(qtplaskin_user_names(i)))
      enddo
    endif
    lerror = .false.
100 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions_list.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    rrt(:) = 1.0d0
    call ZDPlasKin_reac_source_matrix(rrt(:),mrtm(:,:))
    open(ifile_unit,file="qt_matrix.txt",action="write",err=200)
    do i = 1, species_max
      write(ifile_unit,"(143(i3))",err=200) int(mrtm(i,:))
    enddo
    lerror = .false.
200 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_matrix.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_densities.txt",action="write",err=300)
    write(ifile_unit,"(1x,A14,47(111x,i2.2))",err=300) "Time_s", ( i, i = 1, species_max )
    lerror = .false.
300 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_densities.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_conditions.txt",action="write",err=400)
    write(ifile_unit,"(1x,A12,$)",err=400) "Time_s"
    do i = 1, idef_data + iuser_data
      write(ifile_unit,"(11x,i2.2,$)",err=400) i
    enddo
    write(ifile_unit,*,err=400)
    lerror = .false.
400 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_conditions.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
    open(ifile_unit,file="qt_rates.txt",action="write",err=500)
    write(ifile_unit,"(1x,A12,143(101x,i3.3))",err=500) "Time_s", ( i, i = 1, reactions_max )
    lerror = .false.
500 if( lerror ) call ZDPlasKin_stop("ZDPlasKin ERROR: cannot write to file " // &
                                     "<qt_rates.txt> (subroutine writer_save_qtplaskin)")
    close(ifile_unit)
  endif
  if( present(LFORCE_WRITE) ) then
    if( LFORCE_WRITE ) lfirst = .true.
  endif
  rtol = 10.0d0 ** ( floor( log10( abs(densav(0,1)) + tiny(rtol) ) ) - 6 )
  if( ( time - densav(0,1) ) >= rtol .or. lfirst ) then
    densav(0,2) = time
    densav(1:species_max,2) = density(:)
    where( densav(:,2) < 1.0d-99 ) densav(:,2) = 0.0d0
    if( time > 2.0d0 * densav(0,1) ) then
      rtol = huge(rtol)
    else
      rtol = maxval( abs(densav(1:,1)-densav(1:,2)) / ( abs(densav(1:,1)+densav(1:,2))/2.0d0 + qtplaskin_atol ) )
    endif
    if( rtol > qtplaskin_rtol .or. lfirst ) then
      open(ifile_unit,file="qt_densities.txt",access="append")
      write(ifile_unit,"(1pe15.6,47(1pe13.4))") densav(0,2), densav(1:,2)
      close(ifile_unit)
      open(ifile_unit,file="qt_conditions.txt",access="append")
      cond(1) = ZDPlasKin_cfg(3)
      cond(2) = ZDPlasKin_cfg(1)
      cond(3) = ZDPlasKin_cfg(4)
      cond(4) = q_elem * density(species_electrons) * ZDPlasKin_cfg(5) * ZDPlasKin_cfg(3) * 1.0d-17
      call ZDPlasKin_get_density_total(ALL_NEUTRAL=cond(5))
      cond(5) = cond(4) * cond(5) * ZDPlasKin_cfg(3) * 1.0d-17
      where( abs(cond(:)) < 1.0d-99 ) cond(:) = 0.0d0
      write(ifile_unit,"(6(1pe13.4),$)") densav(0,2), cond(:)
      if( iuser_data > 0 ) then
        where( abs(qtplaskin_user_data(1:iuser_data)) < 1.0d-99 ) qtplaskin_user_data(1:iuser_data) = 0.0d0
        write(ifile_unit,"(90(1pe13.4))") qtplaskin_user_data(1:iuser_data)
      else
        write(ifile_unit,*)
      endif
      close(ifile_unit)
      call ZDPlasKin_get_rates(REACTION_RATES=rrt_loc)
      where( abs(rrt_loc(:)) < 1.0d-99 ) rrt_loc(:) = 0.0d0
      open(ifile_unit,file="qt_rates.txt",access="append")
      write(ifile_unit,"(144(1pe13.4))") densav(0,2), rrt_loc(:)
      close(ifile_unit)
      densav(:,1) = densav(:,2)
    endif
  endif
  lfirst = .false.
  lqtplaskin_first = .false.
  return
end subroutine ZDPlasKin_write_qtplaskin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction sensitivity acquisition
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_source_matrix(reac_rate_local,reac_source_local)
  implicit none
  double precision, intent(in)  :: reac_rate_local(reactions_max)
  double precision, intent(out) :: reac_source_local(species_max,reactions_max)
  reac_source_local(:,:) = 0.0d0
  reac_source_local(01,002) = - reac_rate_local(002) 
  reac_source_local(02,002) = - reac_rate_local(002) 
  reac_source_local(15,002) = + reac_rate_local(002) 
  reac_source_local(20,002) = + reac_rate_local(002) 
  reac_source_local(01,003) = - reac_rate_local(003) 
  reac_source_local(02,003) = - reac_rate_local(003) 
  reac_source_local(16,003) = + reac_rate_local(003) 
  reac_source_local(21,003) = + reac_rate_local(003) 
  reac_source_local(01,004) = - reac_rate_local(004) 
  reac_source_local(02,004) = - reac_rate_local(004) 
  reac_source_local(18,004) = + reac_rate_local(004) 
  reac_source_local(22,004) = + reac_rate_local(004) 
  reac_source_local(02,005) = - reac_rate_local(005) 
  reac_source_local(12,005) = + reac_rate_local(005) 
  reac_source_local(18,005) = + reac_rate_local(005) 
  reac_source_local(02,006) = - reac_rate_local(006) 
  reac_source_local(13,006) = + reac_rate_local(006) 
  reac_source_local(18,006) = + reac_rate_local(006) 
  reac_source_local(02,007) = - reac_rate_local(007) 
  reac_source_local(14,007) = + reac_rate_local(007) 
  reac_source_local(18,007) = + reac_rate_local(007) 
  reac_source_local(02,008) = - reac_rate_local(008) 
  reac_source_local(16,008) = + reac_rate_local(008) 
  reac_source_local(34,008) = + reac_rate_local(008) 
  reac_source_local(02,009) = - reac_rate_local(009) 
  reac_source_local(16,009) = + reac_rate_local(009) 
  reac_source_local(35,009) = + reac_rate_local(009) 
  reac_source_local(02,010) = - reac_rate_local(010) 
  reac_source_local(16,010) = + reac_rate_local(010) 
  reac_source_local(30,010) = + reac_rate_local(010) 
  reac_source_local(02,011) = - reac_rate_local(011) 
  reac_source_local(15,011) = + reac_rate_local(011) 
  reac_source_local(18,011) = + reac_rate_local(011) 
  reac_source_local(02,012) = - reac_rate_local(012) 
  reac_source_local(15,012) = + reac_rate_local(012) 
  reac_source_local(19,012) = + reac_rate_local(012) 
  reac_source_local(02,013) = - reac_rate_local(013) 
  reac_source_local(03,013) = + reac_rate_local(013) 
  reac_source_local(02,014) = - reac_rate_local(014) 
  reac_source_local(04,014) = + reac_rate_local(014) 
  reac_source_local(02,015) = - reac_rate_local(015) 
  reac_source_local(05,015) = + reac_rate_local(015) 
  reac_source_local(02,016) = - reac_rate_local(016) 
  reac_source_local(06,016) = + reac_rate_local(016) 
  reac_source_local(02,017) = - reac_rate_local(017) 
  reac_source_local(07,017) = + reac_rate_local(017) 
  reac_source_local(02,018) = - reac_rate_local(018) 
  reac_source_local(08,018) = + reac_rate_local(018) 
  reac_source_local(01,019) = + reac_rate_local(019) 
  reac_source_local(02,019) = - reac_rate_local(019) 
  reac_source_local(23,019) = + reac_rate_local(019) 
  reac_source_local(01,020) = + reac_rate_local(020) 
  reac_source_local(02,020) = - reac_rate_local(020) 
  reac_source_local(18,020) = + reac_rate_local(020) 
  reac_source_local(24,020) = + reac_rate_local(020) 
  reac_source_local(01,021) = + reac_rate_local(021) 
  reac_source_local(02,021) = - reac_rate_local(021) 
  reac_source_local(15,021) = + reac_rate_local(021) 
  reac_source_local(25,021) = + reac_rate_local(021) 
  reac_source_local(01,022) = + reac_rate_local(022) 
  reac_source_local(02,022) = - reac_rate_local(022) 
  reac_source_local(16,022) = + reac_rate_local(022) 
  reac_source_local(26,022) = + reac_rate_local(022) 
  reac_source_local(01,023) = + reac_rate_local(023) 
  reac_source_local(02,023) = - reac_rate_local(023) 
  reac_source_local(27,023) = + reac_rate_local(023) 
  reac_source_local(36,023) = + reac_rate_local(023) 
  reac_source_local(02,024) = - reac_rate_local(024) 
  reac_source_local(16,024) = + reac_rate_local(024) 
  reac_source_local(20,024) = + reac_rate_local(024) 
  reac_source_local(22,024) = - reac_rate_local(024) 
  reac_source_local(02,025) = - reac_rate_local(025) 
  reac_source_local(18,025) = + reac_rate_local(025) 
  reac_source_local(23,025) = - reac_rate_local(025) 
  reac_source_local(28,025) = + reac_rate_local(025) 
  reac_source_local(15,026) = + reac_rate_local(026) 
  reac_source_local(16,026) = - reac_rate_local(026) 
  reac_source_local(23,026) = - reac_rate_local(026) 
  reac_source_local(28,026) = + reac_rate_local(026) 
  reac_source_local(02,027) = + reac_rate_local(027) 
  reac_source_local(17,027) = - reac_rate_local(027) 
  reac_source_local(23,027) = - reac_rate_local(027) 
  reac_source_local(29,027) = + reac_rate_local(027) 
  reac_source_local(02,028) = - reac_rate_local(028) 
  reac_source_local(18,028) = + reac_rate_local(028) 
  reac_source_local(23,028) = + reac_rate_local(028) 
  reac_source_local(25,028) = - reac_rate_local(028) 
  reac_source_local(15,029) = - reac_rate_local(029) 
  reac_source_local(16,029) = + reac_rate_local(029) 
  reac_source_local(25,029) = - reac_rate_local(029) 
  reac_source_local(26,029) = + reac_rate_local(029) 
  reac_source_local(18,030) = + reac_rate_local(030) 
  reac_source_local(25,030) = - reac_rate_local(030) 
  reac_source_local(26,030) = + reac_rate_local(030) 
  reac_source_local(36,030) = - reac_rate_local(030) 
  reac_source_local(02,031) = - reac_rate_local(031) 
  reac_source_local(15,031) = + reac_rate_local(031) 
  reac_source_local(23,031) = + reac_rate_local(031) 
  reac_source_local(24,031) = - reac_rate_local(031) 
  reac_source_local(02,032) = - reac_rate_local(032) 
  reac_source_local(16,032) = + reac_rate_local(032) 
  reac_source_local(23,032) = + reac_rate_local(032) 
  reac_source_local(27,032) = - reac_rate_local(032) 
  reac_source_local(02,033) = - reac_rate_local(033) 
  reac_source_local(15,033) = + reac_rate_local(033) 
  reac_source_local(27,033) = - reac_rate_local(033) 
  reac_source_local(28,033) = + reac_rate_local(033) 
  reac_source_local(02,034) = - reac_rate_local(034) 
  reac_source_local(23,034) = + reac_rate_local(034) 
  reac_source_local(26,034) = - reac_rate_local(034) 
  reac_source_local(36,034) = + reac_rate_local(034) 
  reac_source_local(17,035) = - reac_rate_local(035) 
  reac_source_local(18,035) = + reac_rate_local(035) 
  reac_source_local(25,035) = - reac_rate_local(035) 
  reac_source_local(29,035) = + reac_rate_local(035) 
  reac_source_local(01,036) = + reac_rate_local(036) 
  reac_source_local(17,036) = + reac_rate_local(036) 
  reac_source_local(21,036) = - reac_rate_local(036) 
  reac_source_local(36,036) = - reac_rate_local(036) 
  reac_source_local(01,037) = + reac_rate_local(037) 
  reac_source_local(17,037) = + reac_rate_local(037) 
  reac_source_local(21,037) = - reac_rate_local(037) 
  reac_source_local(30,037) = - reac_rate_local(037) 
  reac_source_local(01,038) = + reac_rate_local(038) 
  reac_source_local(17,038) = + reac_rate_local(038) 
  reac_source_local(21,038) = - reac_rate_local(038) 
  reac_source_local(31,038) = - reac_rate_local(038) 
  reac_source_local(01,039) = + reac_rate_local(039) 
  reac_source_local(11,039) = + reac_rate_local(039) 
  reac_source_local(17,039) = - reac_rate_local(039) 
  reac_source_local(21,039) = - reac_rate_local(039) 
  reac_source_local(01,040) = + reac_rate_local(040) 
  reac_source_local(02,040) = + reac_rate_local(040) 
  reac_source_local(16,040) = - reac_rate_local(040) 
  reac_source_local(21,040) = - reac_rate_local(040) 
  reac_source_local(01,041) = + reac_rate_local(041) 
  reac_source_local(15,041) = - reac_rate_local(041) 
  reac_source_local(16,041) = + reac_rate_local(041) 
  reac_source_local(22,041) = - reac_rate_local(041) 
  reac_source_local(01,042) = + reac_rate_local(042) 
  reac_source_local(09,042) = + reac_rate_local(042) 
  reac_source_local(17,042) = - reac_rate_local(042) 
  reac_source_local(22,042) = - reac_rate_local(042) 
  reac_source_local(01,043) = + reac_rate_local(043) 
  reac_source_local(02,043) = + reac_rate_local(043) 
  reac_source_local(15,043) = - reac_rate_local(043) 
  reac_source_local(20,043) = - reac_rate_local(043) 
  reac_source_local(01,044) = + reac_rate_local(044) 
  reac_source_local(09,044) = + reac_rate_local(044) 
  reac_source_local(20,044) = - reac_rate_local(044) 
  reac_source_local(36,044) = - reac_rate_local(044) 
  reac_source_local(01,045) = + reac_rate_local(045) 
  reac_source_local(09,045) = + reac_rate_local(045) 
  reac_source_local(20,045) = - reac_rate_local(045) 
  reac_source_local(30,045) = - reac_rate_local(045) 
  reac_source_local(01,046) = + reac_rate_local(046) 
  reac_source_local(09,046) = + reac_rate_local(046) 
  reac_source_local(20,046) = - reac_rate_local(046) 
  reac_source_local(31,046) = - reac_rate_local(046) 
  reac_source_local(01,047) = - reac_rate_local(047) 
  reac_source_local(15,047) = + reac_rate_local(047) 
  reac_source_local(18,047) = + reac_rate_local(047) 
  reac_source_local(23,047) = - reac_rate_local(047) 
  reac_source_local(01,048) = - reac_rate_local(048) 
  reac_source_local(16,048) = + reac_rate_local(048) 
  reac_source_local(23,048) = - reac_rate_local(048) 
  reac_source_local(36,048) = + reac_rate_local(048) 
  reac_source_local(01,049) = - reac_rate_local(049) 
  reac_source_local(15,049) = + reac_rate_local(049) * 2.d0
  reac_source_local(23,049) = - reac_rate_local(049) 
  reac_source_local(36,049) = + reac_rate_local(049) 
  reac_source_local(01,050) = - reac_rate_local(050) 
  reac_source_local(15,050) = + reac_rate_local(050) * 2.d0
  reac_source_local(27,050) = - reac_rate_local(050) 
  reac_source_local(01,051) = - reac_rate_local(051) 
  reac_source_local(15,051) = + reac_rate_local(051) 
  reac_source_local(25,051) = - reac_rate_local(051) 
  reac_source_local(30,051) = + reac_rate_local(051) 
  reac_source_local(01,052) = - reac_rate_local(052) 
  reac_source_local(02,052) = + reac_rate_local(052) 
  reac_source_local(15,052) = + reac_rate_local(052) 
  reac_source_local(28,052) = - reac_rate_local(052) 
  reac_source_local(01,053) = - reac_rate_local(053) 
  reac_source_local(29,053) = - reac_rate_local(053) 
  reac_source_local(36,053) = + reac_rate_local(053) * 2.d0
  reac_source_local(01,054) = - reac_rate_local(054) 
  reac_source_local(26,054) = - reac_rate_local(054) 
  reac_source_local(36,054) = + reac_rate_local(054) 
  reac_source_local(01,055) = - reac_rate_local(055) 
  reac_source_local(02,055) = + reac_rate_local(055) 
  reac_source_local(23,055) = - reac_rate_local(055) 
  reac_source_local(01,056) = - reac_rate_local(056) 
  reac_source_local(02,056) = + reac_rate_local(056) 
  reac_source_local(15,056) = + reac_rate_local(056) 
  reac_source_local(28,056) = - reac_rate_local(056) 
  reac_source_local(12,057) = + reac_rate_local(057) 
  reac_source_local(15,057) = + reac_rate_local(057) 
  reac_source_local(22,057) = - reac_rate_local(057) 
  reac_source_local(24,057) = - reac_rate_local(057) 
  reac_source_local(13,058) = + reac_rate_local(058) 
  reac_source_local(15,058) = + reac_rate_local(058) 
  reac_source_local(22,058) = - reac_rate_local(058) 
  reac_source_local(24,058) = - reac_rate_local(058) 
  reac_source_local(15,059) = + reac_rate_local(059) 
  reac_source_local(16,059) = + reac_rate_local(059) 
  reac_source_local(22,059) = - reac_rate_local(059) 
  reac_source_local(27,059) = - reac_rate_local(059) 
  reac_source_local(15,060) = + reac_rate_local(060) 
  reac_source_local(22,060) = - reac_rate_local(060) 
  reac_source_local(26,060) = - reac_rate_local(060) 
  reac_source_local(36,060) = + reac_rate_local(060) 
  reac_source_local(15,061) = + reac_rate_local(061) 
  reac_source_local(18,061) = + reac_rate_local(061) 
  reac_source_local(22,061) = - reac_rate_local(061) 
  reac_source_local(25,061) = - reac_rate_local(061) 
  reac_source_local(02,062) = + reac_rate_local(062) 
  reac_source_local(15,062) = + reac_rate_local(062) 
  reac_source_local(22,062) = - reac_rate_local(062) 
  reac_source_local(23,062) = - reac_rate_local(062) 
  reac_source_local(21,063) = - reac_rate_local(063) 
  reac_source_local(26,063) = - reac_rate_local(063) 
  reac_source_local(36,063) = + reac_rate_local(063) * 2.d0
  reac_source_local(17,064) = + reac_rate_local(064) 
  reac_source_local(21,064) = - reac_rate_local(064) 
  reac_source_local(29,064) = - reac_rate_local(064) 
  reac_source_local(36,064) = + reac_rate_local(064) 
  reac_source_local(15,065) = + reac_rate_local(065) 
  reac_source_local(21,065) = - reac_rate_local(065) 
  reac_source_local(24,065) = - reac_rate_local(065) 
  reac_source_local(36,065) = + reac_rate_local(065) 
  reac_source_local(16,066) = + reac_rate_local(066) 
  reac_source_local(21,066) = - reac_rate_local(066) 
  reac_source_local(27,066) = - reac_rate_local(066) 
  reac_source_local(36,066) = + reac_rate_local(066) 
  reac_source_local(18,067) = + reac_rate_local(067) 
  reac_source_local(21,067) = - reac_rate_local(067) 
  reac_source_local(25,067) = - reac_rate_local(067) 
  reac_source_local(36,067) = + reac_rate_local(067) 
  reac_source_local(02,068) = + reac_rate_local(068) 
  reac_source_local(21,068) = - reac_rate_local(068) 
  reac_source_local(23,068) = - reac_rate_local(068) 
  reac_source_local(36,068) = + reac_rate_local(068) 
  reac_source_local(18,069) = + reac_rate_local(069) 
  reac_source_local(19,069) = + reac_rate_local(069) 
  reac_source_local(20,069) = - reac_rate_local(069) 
  reac_source_local(25,069) = - reac_rate_local(069) 
  reac_source_local(15,070) = + reac_rate_local(070) 
  reac_source_local(18,070) = + reac_rate_local(070) 
  reac_source_local(20,070) = - reac_rate_local(070) 
  reac_source_local(24,070) = - reac_rate_local(070) 
  reac_source_local(16,071) = + reac_rate_local(071) 
  reac_source_local(18,071) = + reac_rate_local(071) 
  reac_source_local(20,071) = - reac_rate_local(071) 
  reac_source_local(27,071) = - reac_rate_local(071) 
  reac_source_local(18,072) = + reac_rate_local(072) 
  reac_source_local(20,072) = - reac_rate_local(072) 
  reac_source_local(26,072) = - reac_rate_local(072) 
  reac_source_local(36,072) = + reac_rate_local(072) 
  reac_source_local(02,073) = + reac_rate_local(073) 
  reac_source_local(18,073) = + reac_rate_local(073) 
  reac_source_local(20,073) = - reac_rate_local(073) 
  reac_source_local(23,073) = - reac_rate_local(073) 
  reac_source_local(02,074) = + reac_rate_local(074) 
  reac_source_local(16,074) = + reac_rate_local(074) 
  reac_source_local(22,074) = - reac_rate_local(074) 
  reac_source_local(28,074) = - reac_rate_local(074) 
  reac_source_local(02,075) = + reac_rate_local(075) 
  reac_source_local(15,075) = + reac_rate_local(075) 
  reac_source_local(21,075) = - reac_rate_local(075) 
  reac_source_local(28,075) = - reac_rate_local(075) 
  reac_source_local(36,075) = + reac_rate_local(075) 
  reac_source_local(02,076) = + reac_rate_local(076) * 2.d0
  reac_source_local(20,076) = - reac_rate_local(076) 
  reac_source_local(28,076) = - reac_rate_local(076) 
  reac_source_local(15,077) = + reac_rate_local(077) 
  reac_source_local(17,077) = + reac_rate_local(077) 
  reac_source_local(22,077) = - reac_rate_local(077) 
  reac_source_local(29,077) = - reac_rate_local(077) 
  reac_source_local(17,078) = + reac_rate_local(078) 
  reac_source_local(18,078) = + reac_rate_local(078) 
  reac_source_local(20,078) = - reac_rate_local(078) 
  reac_source_local(29,078) = - reac_rate_local(078) 
  reac_source_local(02,079) = + reac_rate_local(079) 
  reac_source_local(18,079) = - reac_rate_local(079) * 2.d0
  reac_source_local(36,079) = + reac_rate_local(079) 
  reac_source_local(02,080) = + reac_rate_local(080) 
  reac_source_local(15,080) = + reac_rate_local(080) 
  reac_source_local(16,080) = - reac_rate_local(080) 
  reac_source_local(18,080) = - reac_rate_local(080) 
  reac_source_local(15,081) = + reac_rate_local(081) 
  reac_source_local(17,081) = + reac_rate_local(081) 
  reac_source_local(18,081) = - reac_rate_local(081) 
  reac_source_local(36,081) = - reac_rate_local(081) 
  reac_source_local(09,082) = + reac_rate_local(082) 
  reac_source_local(18,082) = - reac_rate_local(082) 
  reac_source_local(36,082) = - reac_rate_local(082) 
  reac_source_local(15,083) = - reac_rate_local(083) 
  reac_source_local(16,083) = + reac_rate_local(083) 
  reac_source_local(18,083) = - reac_rate_local(083) 
  reac_source_local(36,083) = + reac_rate_local(083) 
  reac_source_local(02,084) = + reac_rate_local(084) 
  reac_source_local(09,084) = + reac_rate_local(084) 
  reac_source_local(10,084) = - reac_rate_local(084) 
  reac_source_local(18,084) = - reac_rate_local(084) 
  reac_source_local(10,085) = + reac_rate_local(085) 
  reac_source_local(18,085) = - reac_rate_local(085) * 2.d0
  reac_source_local(09,086) = + reac_rate_local(086) 
  reac_source_local(11,086) = - reac_rate_local(086) 
  reac_source_local(17,086) = + reac_rate_local(086) 
  reac_source_local(18,086) = - reac_rate_local(086) 
  reac_source_local(09,087) = - reac_rate_local(087) 
  reac_source_local(15,087) = - reac_rate_local(087) 
  reac_source_local(18,087) = + reac_rate_local(087) * 2.d0
  reac_source_local(02,088) = + reac_rate_local(088) 
  reac_source_local(09,088) = - reac_rate_local(088) 
  reac_source_local(15,088) = - reac_rate_local(088) 
  reac_source_local(36,088) = + reac_rate_local(088) 
  reac_source_local(02,089) = + reac_rate_local(089) 
  reac_source_local(09,089) = - reac_rate_local(089) 
  reac_source_local(17,089) = + reac_rate_local(089) 
  reac_source_local(18,089) = - reac_rate_local(089) 
  reac_source_local(09,090) = - reac_rate_local(090) 
  reac_source_local(17,090) = + reac_rate_local(090) 
  reac_source_local(18,090) = + reac_rate_local(090) 
  reac_source_local(36,090) = - reac_rate_local(090) 
  reac_source_local(09,091) = - reac_rate_local(091) 
  reac_source_local(15,091) = - reac_rate_local(091) 
  reac_source_local(16,091) = + reac_rate_local(091) 
  reac_source_local(17,091) = + reac_rate_local(091) 
  reac_source_local(09,092) = - reac_rate_local(092) * 2.d0
  reac_source_local(10,092) = + reac_rate_local(092) 
  reac_source_local(17,092) = + reac_rate_local(092) 
  reac_source_local(09,093) = - reac_rate_local(093) 
  reac_source_local(11,093) = - reac_rate_local(093) 
  reac_source_local(17,093) = + reac_rate_local(093) * 2.d0
  reac_source_local(18,093) = + reac_rate_local(093) 
  reac_source_local(09,094) = + reac_rate_local(094) 
  reac_source_local(10,094) = - reac_rate_local(094) 
  reac_source_local(18,094) = + reac_rate_local(094) 
  reac_source_local(36,094) = - reac_rate_local(094) 
  reac_source_local(02,095) = + reac_rate_local(095) 
  reac_source_local(10,095) = - reac_rate_local(095) 
  reac_source_local(15,095) = - reac_rate_local(095) 
  reac_source_local(18,095) = + reac_rate_local(095) 
  reac_source_local(09,096) = + reac_rate_local(096) 
  reac_source_local(10,096) = - reac_rate_local(096) 
  reac_source_local(15,096) = - reac_rate_local(096) 
  reac_source_local(16,096) = + reac_rate_local(096) 
  reac_source_local(11,097) = - reac_rate_local(097) 
  reac_source_local(17,097) = + reac_rate_local(097) * 2.d0
  reac_source_local(36,097) = - reac_rate_local(097) 
  reac_source_local(15,098) = + reac_rate_local(098) 
  reac_source_local(16,098) = - reac_rate_local(098) 
  reac_source_local(18,098) = + reac_rate_local(098) 
  reac_source_local(36,098) = - reac_rate_local(098) 
  reac_source_local(09,099) = + reac_rate_local(099) 
  reac_source_local(15,099) = - reac_rate_local(099) 
  reac_source_local(17,099) = - reac_rate_local(099) 
  reac_source_local(09,100) = + reac_rate_local(100) 
  reac_source_local(15,100) = - reac_rate_local(100) 
  reac_source_local(17,100) = - reac_rate_local(100) 
  reac_source_local(11,101) = - reac_rate_local(101) 
  reac_source_local(15,101) = - reac_rate_local(101) 
  reac_source_local(17,101) = + reac_rate_local(101) 
  reac_source_local(18,101) = + reac_rate_local(101) 
  reac_source_local(30,102) = - reac_rate_local(102) 
  reac_source_local(36,102) = + reac_rate_local(102) 
  reac_source_local(02,103) = - reac_rate_local(103) 
  reac_source_local(18,103) = + reac_rate_local(103) * 2.d0
  reac_source_local(30,103) = - reac_rate_local(103) 
  reac_source_local(02,104) = - reac_rate_local(104) 
  reac_source_local(18,104) = + reac_rate_local(104) * 2.d0
  reac_source_local(31,104) = - reac_rate_local(104) 
  reac_source_local(30,105) = - reac_rate_local(105) 
  reac_source_local(36,105) = + reac_rate_local(105) 
  reac_source_local(32,106) = + reac_rate_local(106) 
  reac_source_local(34,106) = - reac_rate_local(106) 
  reac_source_local(33,107) = + reac_rate_local(107) 
  reac_source_local(35,107) = - reac_rate_local(107) 
  reac_source_local(31,108) = + reac_rate_local(108) 
  reac_source_local(32,108) = - reac_rate_local(108) 
  reac_source_local(31,109) = + reac_rate_local(109) 
  reac_source_local(33,109) = - reac_rate_local(109) 
  reac_source_local(18,110) = + reac_rate_local(110) 
  reac_source_local(19,110) = - reac_rate_local(110) 
  reac_source_local(12,111) = - reac_rate_local(111) 
  reac_source_local(15,111) = + reac_rate_local(111) 
  reac_source_local(12,112) = + reac_rate_local(112) 
  reac_source_local(13,112) = - reac_rate_local(112) 
  reac_source_local(12,113) = + reac_rate_local(113) 
  reac_source_local(14,113) = - reac_rate_local(113) 
  reac_source_local(02,114) = + reac_rate_local(114) 
  reac_source_local(03,114) = - reac_rate_local(114) 
  reac_source_local(02,115) = + reac_rate_local(115) 
  reac_source_local(04,115) = - reac_rate_local(115) 
  reac_source_local(02,116) = + reac_rate_local(116) 
  reac_source_local(05,116) = - reac_rate_local(116) 
  reac_source_local(02,117) = + reac_rate_local(117) 
  reac_source_local(06,117) = - reac_rate_local(117) 
  reac_source_local(02,118) = + reac_rate_local(118) 
  reac_source_local(07,118) = - reac_rate_local(118) 
  reac_source_local(02,119) = + reac_rate_local(119) 
  reac_source_local(08,119) = - reac_rate_local(119) 
  reac_source_local(15,120) = - reac_rate_local(120) 
  reac_source_local(24,120) = + reac_rate_local(120) 
  reac_source_local(37,120) = + reac_rate_local(120) 
  reac_source_local(38,120) = - reac_rate_local(120) 
  reac_source_local(15,121) = + reac_rate_local(121) 
  reac_source_local(22,121) = - reac_rate_local(121) 
  reac_source_local(37,121) = + reac_rate_local(121) 
  reac_source_local(38,121) = - reac_rate_local(121) 
  reac_source_local(15,122) = + reac_rate_local(122) 
  reac_source_local(16,122) = - reac_rate_local(122) 
  reac_source_local(38,122) = - reac_rate_local(122) 
  reac_source_local(41,122) = + reac_rate_local(122) 
  reac_source_local(18,123) = + reac_rate_local(123) 
  reac_source_local(20,123) = - reac_rate_local(123) 
  reac_source_local(37,123) = + reac_rate_local(123) 
  reac_source_local(38,123) = - reac_rate_local(123) 
  reac_source_local(15,124) = + reac_rate_local(124) 
  reac_source_local(20,124) = - reac_rate_local(124) 
  reac_source_local(36,124) = + reac_rate_local(124) 
  reac_source_local(37,124) = + reac_rate_local(124) 
  reac_source_local(38,124) = - reac_rate_local(124) 
  reac_source_local(02,125) = - reac_rate_local(125) 
  reac_source_local(23,125) = + reac_rate_local(125) 
  reac_source_local(37,125) = + reac_rate_local(125) 
  reac_source_local(38,125) = - reac_rate_local(125) 
  reac_source_local(02,126) = - reac_rate_local(126) 
  reac_source_local(18,126) = + reac_rate_local(126) 
  reac_source_local(38,126) = - reac_rate_local(126) 
  reac_source_local(41,126) = + reac_rate_local(126) 
  reac_source_local(02,127) = + reac_rate_local(127) 
  reac_source_local(23,127) = - reac_rate_local(127) 
  reac_source_local(37,127) = - reac_rate_local(127) 
  reac_source_local(38,127) = + reac_rate_local(127) 
  reac_source_local(02,128) = + reac_rate_local(128) 
  reac_source_local(28,128) = - reac_rate_local(128) 
  reac_source_local(37,128) = - reac_rate_local(128) 
  reac_source_local(41,128) = + reac_rate_local(128) 
  reac_source_local(15,129) = + reac_rate_local(129) 
  reac_source_local(21,129) = - reac_rate_local(129) 
  reac_source_local(36,129) = + reac_rate_local(129) 
  reac_source_local(37,129) = + reac_rate_local(129) 
  reac_source_local(41,129) = - reac_rate_local(129) 
  reac_source_local(01,131) = + reac_rate_local(131) 
  reac_source_local(37,131) = - reac_rate_local(131) 
  reac_source_local(38,131) = + reac_rate_local(131) 
  reac_source_local(37,132) = - reac_rate_local(132) 
  reac_source_local(39,132) = + reac_rate_local(132) 
  reac_source_local(37,133) = + reac_rate_local(133) 
  reac_source_local(39,133) = - reac_rate_local(133) 
  reac_source_local(01,134) = + reac_rate_local(134) 
  reac_source_local(38,134) = + reac_rate_local(134) 
  reac_source_local(39,134) = - reac_rate_local(134) 
  reac_source_local(37,135) = + reac_rate_local(135) 
  reac_source_local(38,135) = + reac_rate_local(135) 
  reac_source_local(40,135) = - reac_rate_local(135) 
  reac_source_local(01,136) = + reac_rate_local(136) 
  reac_source_local(39,136) = - reac_rate_local(136) * 2.d0
  reac_source_local(40,136) = + reac_rate_local(136) 
  reac_source_local(37,137) = + reac_rate_local(137) 
  reac_source_local(39,137) = - reac_rate_local(137) 
  reac_source_local(38,138) = - reac_rate_local(138) 
  reac_source_local(42,138) = + reac_rate_local(138) 
  reac_source_local(40,139) = - reac_rate_local(139) 
  reac_source_local(43,139) = + reac_rate_local(139) 
  reac_source_local(23,140) = - reac_rate_local(140) 
  reac_source_local(45,140) = + reac_rate_local(140) 
  reac_source_local(25,141) = - reac_rate_local(141) 
  reac_source_local(46,141) = + reac_rate_local(141) 
  reac_source_local(22,142) = - reac_rate_local(142) 
  reac_source_local(47,142) = + reac_rate_local(142) 
  reac_source_local(01,143) = - reac_rate_local(143) 
  reac_source_local(44,143) = + reac_rate_local(143) 
  return
end subroutine ZDPlasKin_reac_source_matrix
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction source terms
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_fex(neq,t,y,ydot)
  implicit none
  integer,          intent(in)  :: neq
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: ydot(neq)
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(48)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  rrt(001) = rrt(001) * density(01) * density(02) 
  rrt(002) = rrt(002) * density(01) * density(02) 
  rrt(003) = rrt(003) * density(01) * density(02) 
  rrt(004) = rrt(004) * density(01) * density(02) 
  rrt(005) = rrt(005) * density(01) * density(02) 
  rrt(006) = rrt(006) * density(01) * density(02) 
  rrt(007) = rrt(007) * density(01) * density(02) 
  rrt(008) = rrt(008) * density(01) * density(02) 
  rrt(009) = rrt(009) * density(01) * density(02) 
  rrt(010) = rrt(010) * density(01) * density(02) 
  rrt(011) = rrt(011) * density(01) * density(02) 
  rrt(012) = rrt(012) * density(01) * density(02) 
  rrt(013) = rrt(013) * density(01) * density(02) 
  rrt(014) = rrt(014) * density(01) * density(02) 
  rrt(015) = rrt(015) * density(01) * density(02) 
  rrt(016) = rrt(016) * density(01) * density(02) 
  rrt(017) = rrt(017) * density(01) * density(02) 
  rrt(018) = rrt(018) * density(01) * density(02) 
  rrt(019) = rrt(019) * density(01) * density(02) 
  rrt(020) = rrt(020) * density(01) * density(02) 
  rrt(021) = rrt(021) * density(01) * density(02) 
  rrt(022) = rrt(022) * density(01) * density(02) 
  rrt(023) = rrt(023) * density(01) * density(02) 
  rrt(024) = rrt(024) * density(02) * density(22) 
  rrt(025) = rrt(025) * density(02) * density(23) 
  rrt(026) = rrt(026) * density(16) * density(23) 
  rrt(027) = rrt(027) * density(17) * density(23) 
  rrt(028) = rrt(028) * density(02) * density(25) 
  rrt(029) = rrt(029) * density(15) * density(25) 
  rrt(030) = rrt(030) * density(25) * density(36) 
  rrt(031) = rrt(031) * density(02) * density(24) 
  rrt(032) = rrt(032) * density(02) * density(27) 
  rrt(033) = rrt(033) * density(02) * density(27) 
  rrt(034) = rrt(034) * density(02) * density(26) 
  rrt(035) = rrt(035) * density(17) * density(25) 
  rrt(036) = rrt(036) * density(21) * density(36) 
  rrt(037) = rrt(037) * density(21) * density(30) 
  rrt(038) = rrt(038) * density(21) * density(31) 
  rrt(039) = rrt(039) * density(17) * density(21) 
  rrt(040) = rrt(040) * density(16) * density(21) 
  rrt(041) = rrt(041) * density(15) * density(22) 
  rrt(042) = rrt(042) * density(17) * density(22) 
  rrt(043) = rrt(043) * density(15) * density(20) 
  rrt(044) = rrt(044) * density(20) * density(36) 
  rrt(045) = rrt(045) * density(20) * density(30) 
  rrt(046) = rrt(046) * density(20) * density(31) 
  rrt(047) = rrt(047) * density(01) * density(23) 
  rrt(048) = rrt(048) * density(01) * density(23) 
  rrt(049) = rrt(049) * density(01) * density(23) 
  rrt(050) = rrt(050) * density(01) * density(27) 
  rrt(051) = rrt(051) * density(01) * density(25) 
  rrt(052) = rrt(052) * density(01) * density(28) 
  rrt(053) = rrt(053) * density(01) * density(29) 
  rrt(054) = rrt(054) * density(01) * density(26) 
  rrt(055) = rrt(055) * density(01)**2 * density(23) 
  rrt(056) = rrt(056) * density(01)**2 * density(28) 
  rrt(057) = rrt(057) * density(22) * density(24) 
  rrt(058) = rrt(058) * density(22) * density(24) 
  rrt(059) = rrt(059) * density(22) * density(27) 
  rrt(060) = rrt(060) * density(22) * density(26) 
  rrt(061) = rrt(061) * density(22) * density(25) 
  rrt(062) = rrt(062) * density(22) * density(23) 
  rrt(063) = rrt(063) * density(21) * density(26) 
  rrt(064) = rrt(064) * density(21) * density(29) 
  rrt(065) = rrt(065) * density(21) * density(24) 
  rrt(066) = rrt(066) * density(21) * density(27) 
  rrt(067) = rrt(067) * density(21) * density(25) 
  rrt(068) = rrt(068) * density(21) * density(23) 
  rrt(069) = rrt(069) * density(20) * density(25) 
  rrt(070) = rrt(070) * density(20) * density(24) 
  rrt(071) = rrt(071) * density(20) * density(27) 
  rrt(072) = rrt(072) * density(20) * density(26) 
  rrt(073) = rrt(073) * density(20) * density(23) 
  rrt(074) = rrt(074) * density(22) * density(28) 
  rrt(075) = rrt(075) * density(21) * density(28) 
  rrt(076) = rrt(076) * density(20) * density(28) 
  rrt(077) = rrt(077) * density(22) * density(29) 
  rrt(078) = rrt(078) * density(20) * density(29) 
  rrt(079) = rrt(079) * density(18)**2 
  rrt(080) = rrt(080) * density(16) * density(18) 
  rrt(081) = rrt(081) * density(18) * density(36) 
  rrt(082) = rrt(082) * density(18) * density(36) 
  rrt(083) = rrt(083) * density(15) * density(18) 
  rrt(084) = rrt(084) * density(10) * density(18) 
  rrt(085) = rrt(085) * density(18)**2 
  rrt(086) = rrt(086) * density(11) * density(18) 
  rrt(087) = rrt(087) * density(09) * density(15) 
  rrt(088) = rrt(088) * density(09) * density(15) 
  rrt(089) = rrt(089) * density(09) * density(18) 
  rrt(090) = rrt(090) * density(09) * density(36) 
  rrt(091) = rrt(091) * density(09) * density(15) 
  rrt(092) = rrt(092) * density(09)**2 
  rrt(093) = rrt(093) * density(09) * density(11) 
  rrt(094) = rrt(094) * density(10) * density(36) 
  rrt(095) = rrt(095) * density(10) * density(15) 
  rrt(096) = rrt(096) * density(10) * density(15) 
  rrt(097) = rrt(097) * density(11) * density(36) 
  rrt(098) = rrt(098) * density(16) * density(36) 
  rrt(099) = rrt(099) * density(02) * density(15) * density(17) 
  rrt(100) = rrt(100) * density(15) * density(17) 
  rrt(101) = rrt(101) * density(11) * density(15) 
  rrt(102) = rrt(102) * density(02) * density(30) 
  rrt(103) = rrt(103) * density(02) * density(30) 
  rrt(104) = rrt(104) * density(02) * density(31) 
  rrt(105) = rrt(105) * density(30) 
  rrt(106) = rrt(106) * density(34) 
  rrt(107) = rrt(107) * density(35) 
  rrt(108) = rrt(108) * density(32) 
  rrt(109) = rrt(109) * density(33) 
  rrt(110) = rrt(110) * density(19) 
  rrt(111) = rrt(111) * density(12) 
  rrt(112) = rrt(112) * density(13) 
  rrt(113) = rrt(113) * density(14) 
  rrt(114) = rrt(114) * density(03) 
  rrt(115) = rrt(115) * density(04) 
  rrt(116) = rrt(116) * density(05) 
  rrt(117) = rrt(117) * density(06) 
  rrt(118) = rrt(118) * density(07) 
  rrt(119) = rrt(119) * density(08) 
  rrt(120) = rrt(120) * density(15) * density(38) 
  rrt(121) = rrt(121) * density(22) * density(38) 
  rrt(122) = rrt(122) * density(16) * density(38) 
  rrt(123) = rrt(123) * density(20) * density(38) 
  rrt(124) = rrt(124) * density(20) * density(38) 
  rrt(125) = rrt(125) * density(02) * density(38) 
  rrt(126) = rrt(126) * density(02) * density(38) 
  rrt(127) = rrt(127) * density(23) * density(37) 
  rrt(128) = rrt(128) * density(28) * density(37) 
  rrt(129) = rrt(129) * density(21) * density(41) 
  rrt(130) = rrt(130) * density(01) * density(37) 
  rrt(131) = rrt(131) * density(01) * density(37) 
  rrt(132) = rrt(132) * density(01) * density(37) 
  rrt(133) = rrt(133) * density(01) * density(39) 
  rrt(134) = rrt(134) * density(01) * density(39) 
  rrt(135) = rrt(135) * density(37) * density(40) 
  rrt(136) = rrt(136) * density(39)**2 
  rrt(137) = rrt(137) * density(37)**2 * density(39) 
  rrt(138) = rrt(138) * density(38) 
  rrt(139) = rrt(139) * density(40) 
  rrt(140) = rrt(140) * density(23) 
  rrt(141) = rrt(141) * density(25) 
  rrt(142) = rrt(142) * density(22) 
  rrt(143) = rrt(143) * density(01) 
  ydot(01) = -rrt(002)-rrt(003)-rrt(004)+rrt(019)+rrt(020)+rrt(021)+rrt(022)+rrt(023)+rrt(036)+rrt(037)+rrt(038)+rrt(039)+rrt(040)&
             +rrt(041)+rrt(042)+rrt(043)+rrt(044)+rrt(045)+rrt(046)-rrt(047)-rrt(048)-rrt(049)-rrt(050)-rrt(051)-rrt(052)-rrt(053)&
             -rrt(054)-rrt(055)-rrt(056)+rrt(131)+rrt(134)+rrt(136)-rrt(143) 
  ydot(02) = -rrt(002)-rrt(003)-rrt(004)-rrt(005)-rrt(006)-rrt(007)-rrt(008)-rrt(009)-rrt(010)-rrt(011)-rrt(012)-rrt(013)-rrt(014)&
             -rrt(015)-rrt(016)-rrt(017)-rrt(018)-rrt(019)-rrt(020)-rrt(021)-rrt(022)-rrt(023)-rrt(024)-rrt(025)+rrt(027)-rrt(028)&
             -rrt(031)-rrt(032)-rrt(033)-rrt(034)+rrt(040)+rrt(043)+rrt(052)+rrt(055)+rrt(056)+rrt(062)+rrt(068)+rrt(073)+rrt(074)&
             +rrt(075)+  2.d0 * rrt(076)+rrt(079)+rrt(080)+rrt(084)+rrt(088)+rrt(089)+rrt(095)-rrt(103)-rrt(104)+rrt(114)+rrt(115)&
             +rrt(116)+rrt(117)+rrt(118)+rrt(119)-rrt(125)-rrt(126)+rrt(127)+rrt(128) 
  ydot(03) = +rrt(013)-rrt(114) 
  ydot(04) = +rrt(014)-rrt(115) 
  ydot(05) = +rrt(015)-rrt(116) 
  ydot(06) = +rrt(016)-rrt(117) 
  ydot(07) = +rrt(017)-rrt(118) 
  ydot(08) = +rrt(018)-rrt(119) 
  ydot(09) = +rrt(042)+rrt(044)+rrt(045)+rrt(046)+rrt(082)+rrt(084)+rrt(086)-rrt(087)-rrt(088)-rrt(089)-rrt(090)-rrt(091)&
             -  2.d0 * rrt(092)-rrt(093)+rrt(094)+rrt(096)+rrt(099)+rrt(100) 
  ydot(10) = -rrt(084)+rrt(085)+rrt(092)-rrt(094)-rrt(095)-rrt(096) 
  ydot(11) = +rrt(039)-rrt(086)-rrt(093)-rrt(097)-rrt(101) 
  ydot(12) = +rrt(005)+rrt(057)-rrt(111)+rrt(112)+rrt(113) 
  ydot(13) = +rrt(006)+rrt(058)-rrt(112) 
  ydot(14) = +rrt(007)-rrt(113) 
  ydot(15) = +rrt(002)+rrt(011)+rrt(012)+rrt(021)+rrt(026)-rrt(029)+rrt(031)+rrt(033)-rrt(041)-rrt(043)+rrt(047)+  2.d0 * rrt(049)&
             +  2.d0 * rrt(050)+rrt(051)+rrt(052)+rrt(056)+rrt(057)+rrt(058)+rrt(059)+rrt(060)+rrt(061)+rrt(062)+rrt(065)+rrt(070)&
             +rrt(075)+rrt(077)+rrt(080)+rrt(081)-rrt(083)-rrt(087)-rrt(088)-rrt(091)-rrt(095)-rrt(096)+rrt(098)-rrt(099)-rrt(100)&
             -rrt(101)+rrt(111)-rrt(120)+rrt(121)+rrt(122)+rrt(124)+rrt(129) 
  ydot(16) = +rrt(003)+rrt(008)+rrt(009)+rrt(010)+rrt(022)+rrt(024)-rrt(026)+rrt(029)+rrt(032)-rrt(040)+rrt(041)+rrt(048)+rrt(059)&
             +rrt(066)+rrt(071)+rrt(074)-rrt(080)+rrt(083)+rrt(091)+rrt(096)-rrt(098)-rrt(122) 
  ydot(17) = -rrt(027)-rrt(035)+rrt(036)+rrt(037)+rrt(038)-rrt(039)-rrt(042)+rrt(064)+rrt(077)+rrt(078)+rrt(081)+rrt(086)+rrt(089)&
             +rrt(090)+rrt(091)+rrt(092)+  2.d0 * rrt(093)+  2.d0 * rrt(097)-rrt(099)-rrt(100)+rrt(101) 
  ydot(18) = +rrt(004)+rrt(005)+rrt(006)+rrt(007)+rrt(011)+rrt(020)+rrt(025)+rrt(028)+rrt(030)+rrt(035)+rrt(047)+rrt(061)+rrt(067)&
             +rrt(069)+rrt(070)+rrt(071)+rrt(072)+rrt(073)+rrt(078)-  2.d0 * rrt(079)-rrt(080)-rrt(081)-rrt(082)-rrt(083)-rrt(084)&
             -  2.d0 * rrt(085)-rrt(086)+  2.d0 * rrt(087)-rrt(089)+rrt(090)+rrt(093)+rrt(094)+rrt(095)+rrt(098)+rrt(101)&
             +  2.d0 * rrt(103)+  2.d0 * rrt(104)+rrt(110)+rrt(123)+rrt(126) 
  ydot(19) = +rrt(012)+rrt(069)-rrt(110) 
  ydot(20) = +rrt(002)+rrt(024)-rrt(043)-rrt(044)-rrt(045)-rrt(046)-rrt(069)-rrt(070)-rrt(071)-rrt(072)-rrt(073)-rrt(076)-rrt(078)&
             -rrt(123)-rrt(124) 
  ydot(21) = +rrt(003)-rrt(036)-rrt(037)-rrt(038)-rrt(039)-rrt(040)-rrt(063)-rrt(064)-rrt(065)-rrt(066)-rrt(067)-rrt(068)-rrt(075)&
             -rrt(129) 
  ydot(22) = +rrt(004)-rrt(024)-rrt(041)-rrt(042)-rrt(057)-rrt(058)-rrt(059)-rrt(060)-rrt(061)-rrt(062)-rrt(074)-rrt(077)-rrt(121)&
             -rrt(142) 
  ydot(23) = +rrt(019)-rrt(025)-rrt(026)-rrt(027)+rrt(028)+rrt(031)+rrt(032)+rrt(034)-rrt(047)-rrt(048)-rrt(049)-rrt(055)-rrt(062)&
             -rrt(068)-rrt(073)+rrt(125)-rrt(127)-rrt(140) 
  ydot(24) = +rrt(020)-rrt(031)-rrt(057)-rrt(058)-rrt(065)-rrt(070)+rrt(120) 
  ydot(25) = +rrt(021)-rrt(028)-rrt(029)-rrt(030)-rrt(035)-rrt(051)-rrt(061)-rrt(067)-rrt(069)-rrt(141) 
  ydot(26) = +rrt(022)+rrt(029)+rrt(030)-rrt(034)-rrt(054)-rrt(060)-rrt(063)-rrt(072) 
  ydot(27) = +rrt(023)-rrt(032)-rrt(033)-rrt(050)-rrt(059)-rrt(066)-rrt(071) 
  ydot(28) = +rrt(025)+rrt(026)+rrt(033)-rrt(052)-rrt(056)-rrt(074)-rrt(075)-rrt(076)-rrt(128) 
  ydot(29) = +rrt(027)+rrt(035)-rrt(053)-rrt(064)-rrt(077)-rrt(078) 
  ydot(30) = +rrt(010)-rrt(037)-rrt(045)+rrt(051)-rrt(102)-rrt(103)-rrt(105) 
  ydot(31) = -rrt(038)-rrt(046)-rrt(104)+rrt(108)+rrt(109) 
  ydot(32) = +rrt(106)-rrt(108) 
  ydot(33) = +rrt(107)-rrt(109) 
  ydot(34) = +rrt(008)-rrt(106) 
  ydot(35) = +rrt(009)-rrt(107) 
  ydot(36) = +rrt(023)-rrt(030)+rrt(034)-rrt(036)-rrt(044)+rrt(048)+rrt(049)+  2.d0 * rrt(053)+rrt(054)+rrt(060)+  2.d0 * rrt(063)&
             +rrt(064)+rrt(065)+rrt(066)+rrt(067)+rrt(068)+rrt(072)+rrt(075)+rrt(079)-rrt(081)-rrt(082)+rrt(083)+rrt(088)-rrt(090)&
             -rrt(094)-rrt(097)-rrt(098)+rrt(102)+rrt(105)+rrt(124)+rrt(129) 
  ydot(37) = +rrt(120)+rrt(121)+rrt(123)+rrt(124)+rrt(125)-rrt(127)-rrt(128)+rrt(129)-rrt(131)-rrt(132)+rrt(133)+rrt(135)+rrt(137) 
  ydot(38) = -rrt(120)-rrt(121)-rrt(122)-rrt(123)-rrt(124)-rrt(125)-rrt(126)+rrt(127)+rrt(131)+rrt(134)+rrt(135)-rrt(138) 
  ydot(39) = +rrt(132)-rrt(133)-rrt(134)-  2.d0 * rrt(136)-rrt(137) 
  ydot(40) = -rrt(135)+rrt(136)-rrt(139) 
  ydot(41) = +rrt(122)+rrt(126)+rrt(128)-rrt(129) 
  ydot(42) = +rrt(138) 
  ydot(43) = +rrt(139) 
  ydot(44) = +rrt(143) 
  ydot(45) = +rrt(140) 
  ydot(46) = +rrt(141) 
  ydot(47) = +rrt(142) 
  if( ldensity_constant ) where( density_constant(:) ) ydot(1:species_max) = 0.0d0
  ydot(48) = 0.0d0
  if( lgas_heating ) then
    ydot(48) = ( ZDPlasKin_cfg(14)/k_B + ydot(48) ) / ( sum(density(1:species_max)) - density(species_electrons) ) &
            + eV_to_K * ZDPlasKin_cfg(11) * density(species_electrons)
    ydot(48) = ydot(48) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_fex
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction jacobian
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_jex(neq,t,y,ml,mu,pd,nrpd)
  implicit none
  integer,          intent(in)  :: neq, ml, mu, nrpd
  double precision, intent(in)  :: t, y(neq)
  double precision, intent(out) :: pd(nrpd,neq)
  integer                       :: i
  if( lgas_heating ) ZDPlasKin_cfg(1) = y(48)
  density(:) = y(1:species_max)
  call ZDPlasKin_reac_rates(t)
  pd(01,01) = pd(01,01) - rrt(002) * density(02) 
  pd(01,02) = pd(01,02) - rrt(002) * density(01) 
  pd(02,01) = pd(02,01) - rrt(002) * density(02) 
  pd(02,02) = pd(02,02) - rrt(002) * density(01) 
  pd(15,01) = pd(15,01) + rrt(002) * density(02) 
  pd(15,02) = pd(15,02) + rrt(002) * density(01) 
  pd(20,01) = pd(20,01) + rrt(002) * density(02) 
  pd(20,02) = pd(20,02) + rrt(002) * density(01) 
  pd(01,01) = pd(01,01) - rrt(003) * density(02) 
  pd(01,02) = pd(01,02) - rrt(003) * density(01) 
  pd(02,01) = pd(02,01) - rrt(003) * density(02) 
  pd(02,02) = pd(02,02) - rrt(003) * density(01) 
  pd(16,01) = pd(16,01) + rrt(003) * density(02) 
  pd(16,02) = pd(16,02) + rrt(003) * density(01) 
  pd(21,01) = pd(21,01) + rrt(003) * density(02) 
  pd(21,02) = pd(21,02) + rrt(003) * density(01) 
  pd(01,01) = pd(01,01) - rrt(004) * density(02) 
  pd(01,02) = pd(01,02) - rrt(004) * density(01) 
  pd(02,01) = pd(02,01) - rrt(004) * density(02) 
  pd(02,02) = pd(02,02) - rrt(004) * density(01) 
  pd(18,01) = pd(18,01) + rrt(004) * density(02) 
  pd(18,02) = pd(18,02) + rrt(004) * density(01) 
  pd(22,01) = pd(22,01) + rrt(004) * density(02) 
  pd(22,02) = pd(22,02) + rrt(004) * density(01) 
  pd(02,01) = pd(02,01) - rrt(005) * density(02) 
  pd(02,02) = pd(02,02) - rrt(005) * density(01) 
  pd(12,01) = pd(12,01) + rrt(005) * density(02) 
  pd(12,02) = pd(12,02) + rrt(005) * density(01) 
  pd(18,01) = pd(18,01) + rrt(005) * density(02) 
  pd(18,02) = pd(18,02) + rrt(005) * density(01) 
  pd(02,01) = pd(02,01) - rrt(006) * density(02) 
  pd(02,02) = pd(02,02) - rrt(006) * density(01) 
  pd(13,01) = pd(13,01) + rrt(006) * density(02) 
  pd(13,02) = pd(13,02) + rrt(006) * density(01) 
  pd(18,01) = pd(18,01) + rrt(006) * density(02) 
  pd(18,02) = pd(18,02) + rrt(006) * density(01) 
  pd(02,01) = pd(02,01) - rrt(007) * density(02) 
  pd(02,02) = pd(02,02) - rrt(007) * density(01) 
  pd(14,01) = pd(14,01) + rrt(007) * density(02) 
  pd(14,02) = pd(14,02) + rrt(007) * density(01) 
  pd(18,01) = pd(18,01) + rrt(007) * density(02) 
  pd(18,02) = pd(18,02) + rrt(007) * density(01) 
  pd(02,01) = pd(02,01) - rrt(008) * density(02) 
  pd(02,02) = pd(02,02) - rrt(008) * density(01) 
  pd(16,01) = pd(16,01) + rrt(008) * density(02) 
  pd(16,02) = pd(16,02) + rrt(008) * density(01) 
  pd(34,01) = pd(34,01) + rrt(008) * density(02) 
  pd(34,02) = pd(34,02) + rrt(008) * density(01) 
  pd(02,01) = pd(02,01) - rrt(009) * density(02) 
  pd(02,02) = pd(02,02) - rrt(009) * density(01) 
  pd(16,01) = pd(16,01) + rrt(009) * density(02) 
  pd(16,02) = pd(16,02) + rrt(009) * density(01) 
  pd(35,01) = pd(35,01) + rrt(009) * density(02) 
  pd(35,02) = pd(35,02) + rrt(009) * density(01) 
  pd(02,01) = pd(02,01) - rrt(010) * density(02) 
  pd(02,02) = pd(02,02) - rrt(010) * density(01) 
  pd(16,01) = pd(16,01) + rrt(010) * density(02) 
  pd(16,02) = pd(16,02) + rrt(010) * density(01) 
  pd(30,01) = pd(30,01) + rrt(010) * density(02) 
  pd(30,02) = pd(30,02) + rrt(010) * density(01) 
  pd(02,01) = pd(02,01) - rrt(011) * density(02) 
  pd(02,02) = pd(02,02) - rrt(011) * density(01) 
  pd(15,01) = pd(15,01) + rrt(011) * density(02) 
  pd(15,02) = pd(15,02) + rrt(011) * density(01) 
  pd(18,01) = pd(18,01) + rrt(011) * density(02) 
  pd(18,02) = pd(18,02) + rrt(011) * density(01) 
  pd(02,01) = pd(02,01) - rrt(012) * density(02) 
  pd(02,02) = pd(02,02) - rrt(012) * density(01) 
  pd(15,01) = pd(15,01) + rrt(012) * density(02) 
  pd(15,02) = pd(15,02) + rrt(012) * density(01) 
  pd(19,01) = pd(19,01) + rrt(012) * density(02) 
  pd(19,02) = pd(19,02) + rrt(012) * density(01) 
  pd(02,01) = pd(02,01) - rrt(013) * density(02) 
  pd(02,02) = pd(02,02) - rrt(013) * density(01) 
  pd(03,01) = pd(03,01) + rrt(013) * density(02) 
  pd(03,02) = pd(03,02) + rrt(013) * density(01) 
  pd(02,01) = pd(02,01) - rrt(014) * density(02) 
  pd(02,02) = pd(02,02) - rrt(014) * density(01) 
  pd(04,01) = pd(04,01) + rrt(014) * density(02) 
  pd(04,02) = pd(04,02) + rrt(014) * density(01) 
  pd(02,01) = pd(02,01) - rrt(015) * density(02) 
  pd(02,02) = pd(02,02) - rrt(015) * density(01) 
  pd(05,01) = pd(05,01) + rrt(015) * density(02) 
  pd(05,02) = pd(05,02) + rrt(015) * density(01) 
  pd(02,01) = pd(02,01) - rrt(016) * density(02) 
  pd(02,02) = pd(02,02) - rrt(016) * density(01) 
  pd(06,01) = pd(06,01) + rrt(016) * density(02) 
  pd(06,02) = pd(06,02) + rrt(016) * density(01) 
  pd(02,01) = pd(02,01) - rrt(017) * density(02) 
  pd(02,02) = pd(02,02) - rrt(017) * density(01) 
  pd(07,01) = pd(07,01) + rrt(017) * density(02) 
  pd(07,02) = pd(07,02) + rrt(017) * density(01) 
  pd(02,01) = pd(02,01) - rrt(018) * density(02) 
  pd(02,02) = pd(02,02) - rrt(018) * density(01) 
  pd(08,01) = pd(08,01) + rrt(018) * density(02) 
  pd(08,02) = pd(08,02) + rrt(018) * density(01) 
  pd(01,01) = pd(01,01) + rrt(019) * density(02) 
  pd(01,02) = pd(01,02) + rrt(019) * density(01) 
  pd(02,01) = pd(02,01) - rrt(019) * density(02) 
  pd(02,02) = pd(02,02) - rrt(019) * density(01) 
  pd(23,01) = pd(23,01) + rrt(019) * density(02) 
  pd(23,02) = pd(23,02) + rrt(019) * density(01) 
  pd(01,01) = pd(01,01) + rrt(020) * density(02) 
  pd(01,02) = pd(01,02) + rrt(020) * density(01) 
  pd(02,01) = pd(02,01) - rrt(020) * density(02) 
  pd(02,02) = pd(02,02) - rrt(020) * density(01) 
  pd(18,01) = pd(18,01) + rrt(020) * density(02) 
  pd(18,02) = pd(18,02) + rrt(020) * density(01) 
  pd(24,01) = pd(24,01) + rrt(020) * density(02) 
  pd(24,02) = pd(24,02) + rrt(020) * density(01) 
  pd(01,01) = pd(01,01) + rrt(021) * density(02) 
  pd(01,02) = pd(01,02) + rrt(021) * density(01) 
  pd(02,01) = pd(02,01) - rrt(021) * density(02) 
  pd(02,02) = pd(02,02) - rrt(021) * density(01) 
  pd(15,01) = pd(15,01) + rrt(021) * density(02) 
  pd(15,02) = pd(15,02) + rrt(021) * density(01) 
  pd(25,01) = pd(25,01) + rrt(021) * density(02) 
  pd(25,02) = pd(25,02) + rrt(021) * density(01) 
  pd(01,01) = pd(01,01) + rrt(022) * density(02) 
  pd(01,02) = pd(01,02) + rrt(022) * density(01) 
  pd(02,01) = pd(02,01) - rrt(022) * density(02) 
  pd(02,02) = pd(02,02) - rrt(022) * density(01) 
  pd(16,01) = pd(16,01) + rrt(022) * density(02) 
  pd(16,02) = pd(16,02) + rrt(022) * density(01) 
  pd(26,01) = pd(26,01) + rrt(022) * density(02) 
  pd(26,02) = pd(26,02) + rrt(022) * density(01) 
  pd(01,01) = pd(01,01) + rrt(023) * density(02) 
  pd(01,02) = pd(01,02) + rrt(023) * density(01) 
  pd(02,01) = pd(02,01) - rrt(023) * density(02) 
  pd(02,02) = pd(02,02) - rrt(023) * density(01) 
  pd(27,01) = pd(27,01) + rrt(023) * density(02) 
  pd(27,02) = pd(27,02) + rrt(023) * density(01) 
  pd(36,01) = pd(36,01) + rrt(023) * density(02) 
  pd(36,02) = pd(36,02) + rrt(023) * density(01) 
  pd(02,02) = pd(02,02) - rrt(024) * density(22) 
  pd(02,22) = pd(02,22) - rrt(024) * density(02) 
  pd(16,02) = pd(16,02) + rrt(024) * density(22) 
  pd(16,22) = pd(16,22) + rrt(024) * density(02) 
  pd(20,02) = pd(20,02) + rrt(024) * density(22) 
  pd(20,22) = pd(20,22) + rrt(024) * density(02) 
  pd(22,02) = pd(22,02) - rrt(024) * density(22) 
  pd(22,22) = pd(22,22) - rrt(024) * density(02) 
  pd(02,02) = pd(02,02) - rrt(025) * density(23) 
  pd(02,23) = pd(02,23) - rrt(025) * density(02) 
  pd(18,02) = pd(18,02) + rrt(025) * density(23) 
  pd(18,23) = pd(18,23) + rrt(025) * density(02) 
  pd(23,02) = pd(23,02) - rrt(025) * density(23) 
  pd(23,23) = pd(23,23) - rrt(025) * density(02) 
  pd(28,02) = pd(28,02) + rrt(025) * density(23) 
  pd(28,23) = pd(28,23) + rrt(025) * density(02) 
  pd(15,16) = pd(15,16) + rrt(026) * density(23) 
  pd(15,23) = pd(15,23) + rrt(026) * density(16) 
  pd(16,16) = pd(16,16) - rrt(026) * density(23) 
  pd(16,23) = pd(16,23) - rrt(026) * density(16) 
  pd(23,16) = pd(23,16) - rrt(026) * density(23) 
  pd(23,23) = pd(23,23) - rrt(026) * density(16) 
  pd(28,16) = pd(28,16) + rrt(026) * density(23) 
  pd(28,23) = pd(28,23) + rrt(026) * density(16) 
  pd(02,17) = pd(02,17) + rrt(027) * density(23) 
  pd(02,23) = pd(02,23) + rrt(027) * density(17) 
  pd(17,17) = pd(17,17) - rrt(027) * density(23) 
  pd(17,23) = pd(17,23) - rrt(027) * density(17) 
  pd(23,17) = pd(23,17) - rrt(027) * density(23) 
  pd(23,23) = pd(23,23) - rrt(027) * density(17) 
  pd(29,17) = pd(29,17) + rrt(027) * density(23) 
  pd(29,23) = pd(29,23) + rrt(027) * density(17) 
  pd(02,02) = pd(02,02) - rrt(028) * density(25) 
  pd(02,25) = pd(02,25) - rrt(028) * density(02) 
  pd(18,02) = pd(18,02) + rrt(028) * density(25) 
  pd(18,25) = pd(18,25) + rrt(028) * density(02) 
  pd(23,02) = pd(23,02) + rrt(028) * density(25) 
  pd(23,25) = pd(23,25) + rrt(028) * density(02) 
  pd(25,02) = pd(25,02) - rrt(028) * density(25) 
  pd(25,25) = pd(25,25) - rrt(028) * density(02) 
  pd(15,15) = pd(15,15) - rrt(029) * density(25) 
  pd(15,25) = pd(15,25) - rrt(029) * density(15) 
  pd(16,15) = pd(16,15) + rrt(029) * density(25) 
  pd(16,25) = pd(16,25) + rrt(029) * density(15) 
  pd(25,15) = pd(25,15) - rrt(029) * density(25) 
  pd(25,25) = pd(25,25) - rrt(029) * density(15) 
  pd(26,15) = pd(26,15) + rrt(029) * density(25) 
  pd(26,25) = pd(26,25) + rrt(029) * density(15) 
  pd(18,25) = pd(18,25) + rrt(030) * density(36) 
  pd(18,36) = pd(18,36) + rrt(030) * density(25) 
  pd(25,25) = pd(25,25) - rrt(030) * density(36) 
  pd(25,36) = pd(25,36) - rrt(030) * density(25) 
  pd(26,25) = pd(26,25) + rrt(030) * density(36) 
  pd(26,36) = pd(26,36) + rrt(030) * density(25) 
  pd(36,25) = pd(36,25) - rrt(030) * density(36) 
  pd(36,36) = pd(36,36) - rrt(030) * density(25) 
  pd(02,02) = pd(02,02) - rrt(031) * density(24) 
  pd(02,24) = pd(02,24) - rrt(031) * density(02) 
  pd(15,02) = pd(15,02) + rrt(031) * density(24) 
  pd(15,24) = pd(15,24) + rrt(031) * density(02) 
  pd(23,02) = pd(23,02) + rrt(031) * density(24) 
  pd(23,24) = pd(23,24) + rrt(031) * density(02) 
  pd(24,02) = pd(24,02) - rrt(031) * density(24) 
  pd(24,24) = pd(24,24) - rrt(031) * density(02) 
  pd(02,02) = pd(02,02) - rrt(032) * density(27) 
  pd(02,27) = pd(02,27) - rrt(032) * density(02) 
  pd(16,02) = pd(16,02) + rrt(032) * density(27) 
  pd(16,27) = pd(16,27) + rrt(032) * density(02) 
  pd(23,02) = pd(23,02) + rrt(032) * density(27) 
  pd(23,27) = pd(23,27) + rrt(032) * density(02) 
  pd(27,02) = pd(27,02) - rrt(032) * density(27) 
  pd(27,27) = pd(27,27) - rrt(032) * density(02) 
  pd(02,02) = pd(02,02) - rrt(033) * density(27) 
  pd(02,27) = pd(02,27) - rrt(033) * density(02) 
  pd(15,02) = pd(15,02) + rrt(033) * density(27) 
  pd(15,27) = pd(15,27) + rrt(033) * density(02) 
  pd(27,02) = pd(27,02) - rrt(033) * density(27) 
  pd(27,27) = pd(27,27) - rrt(033) * density(02) 
  pd(28,02) = pd(28,02) + rrt(033) * density(27) 
  pd(28,27) = pd(28,27) + rrt(033) * density(02) 
  pd(02,02) = pd(02,02) - rrt(034) * density(26) 
  pd(02,26) = pd(02,26) - rrt(034) * density(02) 
  pd(23,02) = pd(23,02) + rrt(034) * density(26) 
  pd(23,26) = pd(23,26) + rrt(034) * density(02) 
  pd(26,02) = pd(26,02) - rrt(034) * density(26) 
  pd(26,26) = pd(26,26) - rrt(034) * density(02) 
  pd(36,02) = pd(36,02) + rrt(034) * density(26) 
  pd(36,26) = pd(36,26) + rrt(034) * density(02) 
  pd(17,17) = pd(17,17) - rrt(035) * density(25) 
  pd(17,25) = pd(17,25) - rrt(035) * density(17) 
  pd(18,17) = pd(18,17) + rrt(035) * density(25) 
  pd(18,25) = pd(18,25) + rrt(035) * density(17) 
  pd(25,17) = pd(25,17) - rrt(035) * density(25) 
  pd(25,25) = pd(25,25) - rrt(035) * density(17) 
  pd(29,17) = pd(29,17) + rrt(035) * density(25) 
  pd(29,25) = pd(29,25) + rrt(035) * density(17) 
  pd(01,21) = pd(01,21) + rrt(036) * density(36) 
  pd(01,36) = pd(01,36) + rrt(036) * density(21) 
  pd(17,21) = pd(17,21) + rrt(036) * density(36) 
  pd(17,36) = pd(17,36) + rrt(036) * density(21) 
  pd(21,21) = pd(21,21) - rrt(036) * density(36) 
  pd(21,36) = pd(21,36) - rrt(036) * density(21) 
  pd(36,21) = pd(36,21) - rrt(036) * density(36) 
  pd(36,36) = pd(36,36) - rrt(036) * density(21) 
  pd(01,21) = pd(01,21) + rrt(037) * density(30) 
  pd(01,30) = pd(01,30) + rrt(037) * density(21) 
  pd(17,21) = pd(17,21) + rrt(037) * density(30) 
  pd(17,30) = pd(17,30) + rrt(037) * density(21) 
  pd(21,21) = pd(21,21) - rrt(037) * density(30) 
  pd(21,30) = pd(21,30) - rrt(037) * density(21) 
  pd(30,21) = pd(30,21) - rrt(037) * density(30) 
  pd(30,30) = pd(30,30) - rrt(037) * density(21) 
  pd(01,21) = pd(01,21) + rrt(038) * density(31) 
  pd(01,31) = pd(01,31) + rrt(038) * density(21) 
  pd(17,21) = pd(17,21) + rrt(038) * density(31) 
  pd(17,31) = pd(17,31) + rrt(038) * density(21) 
  pd(21,21) = pd(21,21) - rrt(038) * density(31) 
  pd(21,31) = pd(21,31) - rrt(038) * density(21) 
  pd(31,21) = pd(31,21) - rrt(038) * density(31) 
  pd(31,31) = pd(31,31) - rrt(038) * density(21) 
  pd(01,17) = pd(01,17) + rrt(039) * density(21) 
  pd(01,21) = pd(01,21) + rrt(039) * density(17) 
  pd(11,17) = pd(11,17) + rrt(039) * density(21) 
  pd(11,21) = pd(11,21) + rrt(039) * density(17) 
  pd(17,17) = pd(17,17) - rrt(039) * density(21) 
  pd(17,21) = pd(17,21) - rrt(039) * density(17) 
  pd(21,17) = pd(21,17) - rrt(039) * density(21) 
  pd(21,21) = pd(21,21) - rrt(039) * density(17) 
  pd(01,16) = pd(01,16) + rrt(040) * density(21) 
  pd(01,21) = pd(01,21) + rrt(040) * density(16) 
  pd(02,16) = pd(02,16) + rrt(040) * density(21) 
  pd(02,21) = pd(02,21) + rrt(040) * density(16) 
  pd(16,16) = pd(16,16) - rrt(040) * density(21) 
  pd(16,21) = pd(16,21) - rrt(040) * density(16) 
  pd(21,16) = pd(21,16) - rrt(040) * density(21) 
  pd(21,21) = pd(21,21) - rrt(040) * density(16) 
  pd(01,15) = pd(01,15) + rrt(041) * density(22) 
  pd(01,22) = pd(01,22) + rrt(041) * density(15) 
  pd(15,15) = pd(15,15) - rrt(041) * density(22) 
  pd(15,22) = pd(15,22) - rrt(041) * density(15) 
  pd(16,15) = pd(16,15) + rrt(041) * density(22) 
  pd(16,22) = pd(16,22) + rrt(041) * density(15) 
  pd(22,15) = pd(22,15) - rrt(041) * density(22) 
  pd(22,22) = pd(22,22) - rrt(041) * density(15) 
  pd(01,17) = pd(01,17) + rrt(042) * density(22) 
  pd(01,22) = pd(01,22) + rrt(042) * density(17) 
  pd(09,17) = pd(09,17) + rrt(042) * density(22) 
  pd(09,22) = pd(09,22) + rrt(042) * density(17) 
  pd(17,17) = pd(17,17) - rrt(042) * density(22) 
  pd(17,22) = pd(17,22) - rrt(042) * density(17) 
  pd(22,17) = pd(22,17) - rrt(042) * density(22) 
  pd(22,22) = pd(22,22) - rrt(042) * density(17) 
  pd(01,15) = pd(01,15) + rrt(043) * density(20) 
  pd(01,20) = pd(01,20) + rrt(043) * density(15) 
  pd(02,15) = pd(02,15) + rrt(043) * density(20) 
  pd(02,20) = pd(02,20) + rrt(043) * density(15) 
  pd(15,15) = pd(15,15) - rrt(043) * density(20) 
  pd(15,20) = pd(15,20) - rrt(043) * density(15) 
  pd(20,15) = pd(20,15) - rrt(043) * density(20) 
  pd(20,20) = pd(20,20) - rrt(043) * density(15) 
  pd(01,20) = pd(01,20) + rrt(044) * density(36) 
  pd(01,36) = pd(01,36) + rrt(044) * density(20) 
  pd(09,20) = pd(09,20) + rrt(044) * density(36) 
  pd(09,36) = pd(09,36) + rrt(044) * density(20) 
  pd(20,20) = pd(20,20) - rrt(044) * density(36) 
  pd(20,36) = pd(20,36) - rrt(044) * density(20) 
  pd(36,20) = pd(36,20) - rrt(044) * density(36) 
  pd(36,36) = pd(36,36) - rrt(044) * density(20) 
  pd(01,20) = pd(01,20) + rrt(045) * density(30) 
  pd(01,30) = pd(01,30) + rrt(045) * density(20) 
  pd(09,20) = pd(09,20) + rrt(045) * density(30) 
  pd(09,30) = pd(09,30) + rrt(045) * density(20) 
  pd(20,20) = pd(20,20) - rrt(045) * density(30) 
  pd(20,30) = pd(20,30) - rrt(045) * density(20) 
  pd(30,20) = pd(30,20) - rrt(045) * density(30) 
  pd(30,30) = pd(30,30) - rrt(045) * density(20) 
  pd(01,20) = pd(01,20) + rrt(046) * density(31) 
  pd(01,31) = pd(01,31) + rrt(046) * density(20) 
  pd(09,20) = pd(09,20) + rrt(046) * density(31) 
  pd(09,31) = pd(09,31) + rrt(046) * density(20) 
  pd(20,20) = pd(20,20) - rrt(046) * density(31) 
  pd(20,31) = pd(20,31) - rrt(046) * density(20) 
  pd(31,20) = pd(31,20) - rrt(046) * density(31) 
  pd(31,31) = pd(31,31) - rrt(046) * density(20) 
  pd(01,01) = pd(01,01) - rrt(047) * density(23) 
  pd(01,23) = pd(01,23) - rrt(047) * density(01) 
  pd(15,01) = pd(15,01) + rrt(047) * density(23) 
  pd(15,23) = pd(15,23) + rrt(047) * density(01) 
  pd(18,01) = pd(18,01) + rrt(047) * density(23) 
  pd(18,23) = pd(18,23) + rrt(047) * density(01) 
  pd(23,01) = pd(23,01) - rrt(047) * density(23) 
  pd(23,23) = pd(23,23) - rrt(047) * density(01) 
  pd(01,01) = pd(01,01) - rrt(048) * density(23) 
  pd(01,23) = pd(01,23) - rrt(048) * density(01) 
  pd(16,01) = pd(16,01) + rrt(048) * density(23) 
  pd(16,23) = pd(16,23) + rrt(048) * density(01) 
  pd(23,01) = pd(23,01) - rrt(048) * density(23) 
  pd(23,23) = pd(23,23) - rrt(048) * density(01) 
  pd(36,01) = pd(36,01) + rrt(048) * density(23) 
  pd(36,23) = pd(36,23) + rrt(048) * density(01) 
  pd(01,01) = pd(01,01) - rrt(049) * density(23) 
  pd(01,23) = pd(01,23) - rrt(049) * density(01) 
  pd(15,01) = pd(15,01) + rrt(049) * density(23) * 2.0d0
  pd(15,23) = pd(15,23) + rrt(049) * density(01) * 2.0d0
  pd(23,01) = pd(23,01) - rrt(049) * density(23) 
  pd(23,23) = pd(23,23) - rrt(049) * density(01) 
  pd(36,01) = pd(36,01) + rrt(049) * density(23) 
  pd(36,23) = pd(36,23) + rrt(049) * density(01) 
  pd(01,01) = pd(01,01) - rrt(050) * density(27) 
  pd(01,27) = pd(01,27) - rrt(050) * density(01) 
  pd(15,01) = pd(15,01) + rrt(050) * density(27) * 2.0d0
  pd(15,27) = pd(15,27) + rrt(050) * density(01) * 2.0d0
  pd(27,01) = pd(27,01) - rrt(050) * density(27) 
  pd(27,27) = pd(27,27) - rrt(050) * density(01) 
  pd(01,01) = pd(01,01) - rrt(051) * density(25) 
  pd(01,25) = pd(01,25) - rrt(051) * density(01) 
  pd(15,01) = pd(15,01) + rrt(051) * density(25) 
  pd(15,25) = pd(15,25) + rrt(051) * density(01) 
  pd(25,01) = pd(25,01) - rrt(051) * density(25) 
  pd(25,25) = pd(25,25) - rrt(051) * density(01) 
  pd(30,01) = pd(30,01) + rrt(051) * density(25) 
  pd(30,25) = pd(30,25) + rrt(051) * density(01) 
  pd(01,01) = pd(01,01) - rrt(052) * density(28) 
  pd(01,28) = pd(01,28) - rrt(052) * density(01) 
  pd(02,01) = pd(02,01) + rrt(052) * density(28) 
  pd(02,28) = pd(02,28) + rrt(052) * density(01) 
  pd(15,01) = pd(15,01) + rrt(052) * density(28) 
  pd(15,28) = pd(15,28) + rrt(052) * density(01) 
  pd(28,01) = pd(28,01) - rrt(052) * density(28) 
  pd(28,28) = pd(28,28) - rrt(052) * density(01) 
  pd(01,01) = pd(01,01) - rrt(053) * density(29) 
  pd(01,29) = pd(01,29) - rrt(053) * density(01) 
  pd(29,01) = pd(29,01) - rrt(053) * density(29) 
  pd(29,29) = pd(29,29) - rrt(053) * density(01) 
  pd(36,01) = pd(36,01) + rrt(053) * density(29) * 2.0d0
  pd(36,29) = pd(36,29) + rrt(053) * density(01) * 2.0d0
  pd(01,01) = pd(01,01) - rrt(054) * density(26) 
  pd(01,26) = pd(01,26) - rrt(054) * density(01) 
  pd(26,01) = pd(26,01) - rrt(054) * density(26) 
  pd(26,26) = pd(26,26) - rrt(054) * density(01) 
  pd(36,01) = pd(36,01) + rrt(054) * density(26) 
  pd(36,26) = pd(36,26) + rrt(054) * density(01) 
  pd(01,01) = pd(01,01) - rrt(055) * density(01) * density(23) * 2.0d0
  pd(01,23) = pd(01,23) - rrt(055) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(055) * density(01) * density(23) * 2.0d0
  pd(02,23) = pd(02,23) + rrt(055) * density(01)**2 
  pd(23,01) = pd(23,01) - rrt(055) * density(01) * density(23) * 2.0d0
  pd(23,23) = pd(23,23) - rrt(055) * density(01)**2 
  pd(01,01) = pd(01,01) - rrt(056) * density(01) * density(28) * 2.0d0
  pd(01,28) = pd(01,28) - rrt(056) * density(01)**2 
  pd(02,01) = pd(02,01) + rrt(056) * density(01) * density(28) * 2.0d0
  pd(02,28) = pd(02,28) + rrt(056) * density(01)**2 
  pd(15,01) = pd(15,01) + rrt(056) * density(01) * density(28) * 2.0d0
  pd(15,28) = pd(15,28) + rrt(056) * density(01)**2 
  pd(28,01) = pd(28,01) - rrt(056) * density(01) * density(28) * 2.0d0
  pd(28,28) = pd(28,28) - rrt(056) * density(01)**2 
  pd(12,22) = pd(12,22) + rrt(057) * density(24) 
  pd(12,24) = pd(12,24) + rrt(057) * density(22) 
  pd(15,22) = pd(15,22) + rrt(057) * density(24) 
  pd(15,24) = pd(15,24) + rrt(057) * density(22) 
  pd(22,22) = pd(22,22) - rrt(057) * density(24) 
  pd(22,24) = pd(22,24) - rrt(057) * density(22) 
  pd(24,22) = pd(24,22) - rrt(057) * density(24) 
  pd(24,24) = pd(24,24) - rrt(057) * density(22) 
  pd(13,22) = pd(13,22) + rrt(058) * density(24) 
  pd(13,24) = pd(13,24) + rrt(058) * density(22) 
  pd(15,22) = pd(15,22) + rrt(058) * density(24) 
  pd(15,24) = pd(15,24) + rrt(058) * density(22) 
  pd(22,22) = pd(22,22) - rrt(058) * density(24) 
  pd(22,24) = pd(22,24) - rrt(058) * density(22) 
  pd(24,22) = pd(24,22) - rrt(058) * density(24) 
  pd(24,24) = pd(24,24) - rrt(058) * density(22) 
  pd(15,22) = pd(15,22) + rrt(059) * density(27) 
  pd(15,27) = pd(15,27) + rrt(059) * density(22) 
  pd(16,22) = pd(16,22) + rrt(059) * density(27) 
  pd(16,27) = pd(16,27) + rrt(059) * density(22) 
  pd(22,22) = pd(22,22) - rrt(059) * density(27) 
  pd(22,27) = pd(22,27) - rrt(059) * density(22) 
  pd(27,22) = pd(27,22) - rrt(059) * density(27) 
  pd(27,27) = pd(27,27) - rrt(059) * density(22) 
  pd(15,22) = pd(15,22) + rrt(060) * density(26) 
  pd(15,26) = pd(15,26) + rrt(060) * density(22) 
  pd(22,22) = pd(22,22) - rrt(060) * density(26) 
  pd(22,26) = pd(22,26) - rrt(060) * density(22) 
  pd(26,22) = pd(26,22) - rrt(060) * density(26) 
  pd(26,26) = pd(26,26) - rrt(060) * density(22) 
  pd(36,22) = pd(36,22) + rrt(060) * density(26) 
  pd(36,26) = pd(36,26) + rrt(060) * density(22) 
  pd(15,22) = pd(15,22) + rrt(061) * density(25) 
  pd(15,25) = pd(15,25) + rrt(061) * density(22) 
  pd(18,22) = pd(18,22) + rrt(061) * density(25) 
  pd(18,25) = pd(18,25) + rrt(061) * density(22) 
  pd(22,22) = pd(22,22) - rrt(061) * density(25) 
  pd(22,25) = pd(22,25) - rrt(061) * density(22) 
  pd(25,22) = pd(25,22) - rrt(061) * density(25) 
  pd(25,25) = pd(25,25) - rrt(061) * density(22) 
  pd(02,22) = pd(02,22) + rrt(062) * density(23) 
  pd(02,23) = pd(02,23) + rrt(062) * density(22) 
  pd(15,22) = pd(15,22) + rrt(062) * density(23) 
  pd(15,23) = pd(15,23) + rrt(062) * density(22) 
  pd(22,22) = pd(22,22) - rrt(062) * density(23) 
  pd(22,23) = pd(22,23) - rrt(062) * density(22) 
  pd(23,22) = pd(23,22) - rrt(062) * density(23) 
  pd(23,23) = pd(23,23) - rrt(062) * density(22) 
  pd(21,21) = pd(21,21) - rrt(063) * density(26) 
  pd(21,26) = pd(21,26) - rrt(063) * density(21) 
  pd(26,21) = pd(26,21) - rrt(063) * density(26) 
  pd(26,26) = pd(26,26) - rrt(063) * density(21) 
  pd(36,21) = pd(36,21) + rrt(063) * density(26) * 2.0d0
  pd(36,26) = pd(36,26) + rrt(063) * density(21) * 2.0d0
  pd(17,21) = pd(17,21) + rrt(064) * density(29) 
  pd(17,29) = pd(17,29) + rrt(064) * density(21) 
  pd(21,21) = pd(21,21) - rrt(064) * density(29) 
  pd(21,29) = pd(21,29) - rrt(064) * density(21) 
  pd(29,21) = pd(29,21) - rrt(064) * density(29) 
  pd(29,29) = pd(29,29) - rrt(064) * density(21) 
  pd(36,21) = pd(36,21) + rrt(064) * density(29) 
  pd(36,29) = pd(36,29) + rrt(064) * density(21) 
  pd(15,21) = pd(15,21) + rrt(065) * density(24) 
  pd(15,24) = pd(15,24) + rrt(065) * density(21) 
  pd(21,21) = pd(21,21) - rrt(065) * density(24) 
  pd(21,24) = pd(21,24) - rrt(065) * density(21) 
  pd(24,21) = pd(24,21) - rrt(065) * density(24) 
  pd(24,24) = pd(24,24) - rrt(065) * density(21) 
  pd(36,21) = pd(36,21) + rrt(065) * density(24) 
  pd(36,24) = pd(36,24) + rrt(065) * density(21) 
  pd(16,21) = pd(16,21) + rrt(066) * density(27) 
  pd(16,27) = pd(16,27) + rrt(066) * density(21) 
  pd(21,21) = pd(21,21) - rrt(066) * density(27) 
  pd(21,27) = pd(21,27) - rrt(066) * density(21) 
  pd(27,21) = pd(27,21) - rrt(066) * density(27) 
  pd(27,27) = pd(27,27) - rrt(066) * density(21) 
  pd(36,21) = pd(36,21) + rrt(066) * density(27) 
  pd(36,27) = pd(36,27) + rrt(066) * density(21) 
  pd(18,21) = pd(18,21) + rrt(067) * density(25) 
  pd(18,25) = pd(18,25) + rrt(067) * density(21) 
  pd(21,21) = pd(21,21) - rrt(067) * density(25) 
  pd(21,25) = pd(21,25) - rrt(067) * density(21) 
  pd(25,21) = pd(25,21) - rrt(067) * density(25) 
  pd(25,25) = pd(25,25) - rrt(067) * density(21) 
  pd(36,21) = pd(36,21) + rrt(067) * density(25) 
  pd(36,25) = pd(36,25) + rrt(067) * density(21) 
  pd(02,21) = pd(02,21) + rrt(068) * density(23) 
  pd(02,23) = pd(02,23) + rrt(068) * density(21) 
  pd(21,21) = pd(21,21) - rrt(068) * density(23) 
  pd(21,23) = pd(21,23) - rrt(068) * density(21) 
  pd(23,21) = pd(23,21) - rrt(068) * density(23) 
  pd(23,23) = pd(23,23) - rrt(068) * density(21) 
  pd(36,21) = pd(36,21) + rrt(068) * density(23) 
  pd(36,23) = pd(36,23) + rrt(068) * density(21) 
  pd(18,20) = pd(18,20) + rrt(069) * density(25) 
  pd(18,25) = pd(18,25) + rrt(069) * density(20) 
  pd(19,20) = pd(19,20) + rrt(069) * density(25) 
  pd(19,25) = pd(19,25) + rrt(069) * density(20) 
  pd(20,20) = pd(20,20) - rrt(069) * density(25) 
  pd(20,25) = pd(20,25) - rrt(069) * density(20) 
  pd(25,20) = pd(25,20) - rrt(069) * density(25) 
  pd(25,25) = pd(25,25) - rrt(069) * density(20) 
  pd(15,20) = pd(15,20) + rrt(070) * density(24) 
  pd(15,24) = pd(15,24) + rrt(070) * density(20) 
  pd(18,20) = pd(18,20) + rrt(070) * density(24) 
  pd(18,24) = pd(18,24) + rrt(070) * density(20) 
  pd(20,20) = pd(20,20) - rrt(070) * density(24) 
  pd(20,24) = pd(20,24) - rrt(070) * density(20) 
  pd(24,20) = pd(24,20) - rrt(070) * density(24) 
  pd(24,24) = pd(24,24) - rrt(070) * density(20) 
  pd(16,20) = pd(16,20) + rrt(071) * density(27) 
  pd(16,27) = pd(16,27) + rrt(071) * density(20) 
  pd(18,20) = pd(18,20) + rrt(071) * density(27) 
  pd(18,27) = pd(18,27) + rrt(071) * density(20) 
  pd(20,20) = pd(20,20) - rrt(071) * density(27) 
  pd(20,27) = pd(20,27) - rrt(071) * density(20) 
  pd(27,20) = pd(27,20) - rrt(071) * density(27) 
  pd(27,27) = pd(27,27) - rrt(071) * density(20) 
  pd(18,20) = pd(18,20) + rrt(072) * density(26) 
  pd(18,26) = pd(18,26) + rrt(072) * density(20) 
  pd(20,20) = pd(20,20) - rrt(072) * density(26) 
  pd(20,26) = pd(20,26) - rrt(072) * density(20) 
  pd(26,20) = pd(26,20) - rrt(072) * density(26) 
  pd(26,26) = pd(26,26) - rrt(072) * density(20) 
  pd(36,20) = pd(36,20) + rrt(072) * density(26) 
  pd(36,26) = pd(36,26) + rrt(072) * density(20) 
  pd(02,20) = pd(02,20) + rrt(073) * density(23) 
  pd(02,23) = pd(02,23) + rrt(073) * density(20) 
  pd(18,20) = pd(18,20) + rrt(073) * density(23) 
  pd(18,23) = pd(18,23) + rrt(073) * density(20) 
  pd(20,20) = pd(20,20) - rrt(073) * density(23) 
  pd(20,23) = pd(20,23) - rrt(073) * density(20) 
  pd(23,20) = pd(23,20) - rrt(073) * density(23) 
  pd(23,23) = pd(23,23) - rrt(073) * density(20) 
  pd(02,22) = pd(02,22) + rrt(074) * density(28) 
  pd(02,28) = pd(02,28) + rrt(074) * density(22) 
  pd(16,22) = pd(16,22) + rrt(074) * density(28) 
  pd(16,28) = pd(16,28) + rrt(074) * density(22) 
  pd(22,22) = pd(22,22) - rrt(074) * density(28) 
  pd(22,28) = pd(22,28) - rrt(074) * density(22) 
  pd(28,22) = pd(28,22) - rrt(074) * density(28) 
  pd(28,28) = pd(28,28) - rrt(074) * density(22) 
  pd(02,21) = pd(02,21) + rrt(075) * density(28) 
  pd(02,28) = pd(02,28) + rrt(075) * density(21) 
  pd(15,21) = pd(15,21) + rrt(075) * density(28) 
  pd(15,28) = pd(15,28) + rrt(075) * density(21) 
  pd(21,21) = pd(21,21) - rrt(075) * density(28) 
  pd(21,28) = pd(21,28) - rrt(075) * density(21) 
  pd(28,21) = pd(28,21) - rrt(075) * density(28) 
  pd(28,28) = pd(28,28) - rrt(075) * density(21) 
  pd(36,21) = pd(36,21) + rrt(075) * density(28) 
  pd(36,28) = pd(36,28) + rrt(075) * density(21) 
  pd(02,20) = pd(02,20) + rrt(076) * density(28) * 2.0d0
  pd(02,28) = pd(02,28) + rrt(076) * density(20) * 2.0d0
  pd(20,20) = pd(20,20) - rrt(076) * density(28) 
  pd(20,28) = pd(20,28) - rrt(076) * density(20) 
  pd(28,20) = pd(28,20) - rrt(076) * density(28) 
  pd(28,28) = pd(28,28) - rrt(076) * density(20) 
  pd(15,22) = pd(15,22) + rrt(077) * density(29) 
  pd(15,29) = pd(15,29) + rrt(077) * density(22) 
  pd(17,22) = pd(17,22) + rrt(077) * density(29) 
  pd(17,29) = pd(17,29) + rrt(077) * density(22) 
  pd(22,22) = pd(22,22) - rrt(077) * density(29) 
  pd(22,29) = pd(22,29) - rrt(077) * density(22) 
  pd(29,22) = pd(29,22) - rrt(077) * density(29) 
  pd(29,29) = pd(29,29) - rrt(077) * density(22) 
  pd(17,20) = pd(17,20) + rrt(078) * density(29) 
  pd(17,29) = pd(17,29) + rrt(078) * density(20) 
  pd(18,20) = pd(18,20) + rrt(078) * density(29) 
  pd(18,29) = pd(18,29) + rrt(078) * density(20) 
  pd(20,20) = pd(20,20) - rrt(078) * density(29) 
  pd(20,29) = pd(20,29) - rrt(078) * density(20) 
  pd(29,20) = pd(29,20) - rrt(078) * density(29) 
  pd(29,29) = pd(29,29) - rrt(078) * density(20) 
  pd(02,18) = pd(02,18) + rrt(079) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(079) * density(18) * 4.0d0
  pd(36,18) = pd(36,18) + rrt(079) * density(18) * 2.0d0
  pd(02,16) = pd(02,16) + rrt(080) * density(18) 
  pd(02,18) = pd(02,18) + rrt(080) * density(16) 
  pd(15,16) = pd(15,16) + rrt(080) * density(18) 
  pd(15,18) = pd(15,18) + rrt(080) * density(16) 
  pd(16,16) = pd(16,16) - rrt(080) * density(18) 
  pd(16,18) = pd(16,18) - rrt(080) * density(16) 
  pd(18,16) = pd(18,16) - rrt(080) * density(18) 
  pd(18,18) = pd(18,18) - rrt(080) * density(16) 
  pd(15,18) = pd(15,18) + rrt(081) * density(36) 
  pd(15,36) = pd(15,36) + rrt(081) * density(18) 
  pd(17,18) = pd(17,18) + rrt(081) * density(36) 
  pd(17,36) = pd(17,36) + rrt(081) * density(18) 
  pd(18,18) = pd(18,18) - rrt(081) * density(36) 
  pd(18,36) = pd(18,36) - rrt(081) * density(18) 
  pd(36,18) = pd(36,18) - rrt(081) * density(36) 
  pd(36,36) = pd(36,36) - rrt(081) * density(18) 
  pd(09,18) = pd(09,18) + rrt(082) * density(36) 
  pd(09,36) = pd(09,36) + rrt(082) * density(18) 
  pd(18,18) = pd(18,18) - rrt(082) * density(36) 
  pd(18,36) = pd(18,36) - rrt(082) * density(18) 
  pd(36,18) = pd(36,18) - rrt(082) * density(36) 
  pd(36,36) = pd(36,36) - rrt(082) * density(18) 
  pd(15,15) = pd(15,15) - rrt(083) * density(18) 
  pd(15,18) = pd(15,18) - rrt(083) * density(15) 
  pd(16,15) = pd(16,15) + rrt(083) * density(18) 
  pd(16,18) = pd(16,18) + rrt(083) * density(15) 
  pd(18,15) = pd(18,15) - rrt(083) * density(18) 
  pd(18,18) = pd(18,18) - rrt(083) * density(15) 
  pd(36,15) = pd(36,15) + rrt(083) * density(18) 
  pd(36,18) = pd(36,18) + rrt(083) * density(15) 
  pd(02,10) = pd(02,10) + rrt(084) * density(18) 
  pd(02,18) = pd(02,18) + rrt(084) * density(10) 
  pd(09,10) = pd(09,10) + rrt(084) * density(18) 
  pd(09,18) = pd(09,18) + rrt(084) * density(10) 
  pd(10,10) = pd(10,10) - rrt(084) * density(18) 
  pd(10,18) = pd(10,18) - rrt(084) * density(10) 
  pd(18,10) = pd(18,10) - rrt(084) * density(18) 
  pd(18,18) = pd(18,18) - rrt(084) * density(10) 
  pd(10,18) = pd(10,18) + rrt(085) * density(18) * 2.0d0
  pd(18,18) = pd(18,18) - rrt(085) * density(18) * 4.0d0
  pd(09,11) = pd(09,11) + rrt(086) * density(18) 
  pd(09,18) = pd(09,18) + rrt(086) * density(11) 
  pd(11,11) = pd(11,11) - rrt(086) * density(18) 
  pd(11,18) = pd(11,18) - rrt(086) * density(11) 
  pd(17,11) = pd(17,11) + rrt(086) * density(18) 
  pd(17,18) = pd(17,18) + rrt(086) * density(11) 
  pd(18,11) = pd(18,11) - rrt(086) * density(18) 
  pd(18,18) = pd(18,18) - rrt(086) * density(11) 
  pd(09,09) = pd(09,09) - rrt(087) * density(15) 
  pd(09,15) = pd(09,15) - rrt(087) * density(09) 
  pd(15,09) = pd(15,09) - rrt(087) * density(15) 
  pd(15,15) = pd(15,15) - rrt(087) * density(09) 
  pd(18,09) = pd(18,09) + rrt(087) * density(15) * 2.0d0
  pd(18,15) = pd(18,15) + rrt(087) * density(09) * 2.0d0
  pd(02,09) = pd(02,09) + rrt(088) * density(15) 
  pd(02,15) = pd(02,15) + rrt(088) * density(09) 
  pd(09,09) = pd(09,09) - rrt(088) * density(15) 
  pd(09,15) = pd(09,15) - rrt(088) * density(09) 
  pd(15,09) = pd(15,09) - rrt(088) * density(15) 
  pd(15,15) = pd(15,15) - rrt(088) * density(09) 
  pd(36,09) = pd(36,09) + rrt(088) * density(15) 
  pd(36,15) = pd(36,15) + rrt(088) * density(09) 
  pd(02,09) = pd(02,09) + rrt(089) * density(18) 
  pd(02,18) = pd(02,18) + rrt(089) * density(09) 
  pd(09,09) = pd(09,09) - rrt(089) * density(18) 
  pd(09,18) = pd(09,18) - rrt(089) * density(09) 
  pd(17,09) = pd(17,09) + rrt(089) * density(18) 
  pd(17,18) = pd(17,18) + rrt(089) * density(09) 
  pd(18,09) = pd(18,09) - rrt(089) * density(18) 
  pd(18,18) = pd(18,18) - rrt(089) * density(09) 
  pd(09,09) = pd(09,09) - rrt(090) * density(36) 
  pd(09,36) = pd(09,36) - rrt(090) * density(09) 
  pd(17,09) = pd(17,09) + rrt(090) * density(36) 
  pd(17,36) = pd(17,36) + rrt(090) * density(09) 
  pd(18,09) = pd(18,09) + rrt(090) * density(36) 
  pd(18,36) = pd(18,36) + rrt(090) * density(09) 
  pd(36,09) = pd(36,09) - rrt(090) * density(36) 
  pd(36,36) = pd(36,36) - rrt(090) * density(09) 
  pd(09,09) = pd(09,09) - rrt(091) * density(15) 
  pd(09,15) = pd(09,15) - rrt(091) * density(09) 
  pd(15,09) = pd(15,09) - rrt(091) * density(15) 
  pd(15,15) = pd(15,15) - rrt(091) * density(09) 
  pd(16,09) = pd(16,09) + rrt(091) * density(15) 
  pd(16,15) = pd(16,15) + rrt(091) * density(09) 
  pd(17,09) = pd(17,09) + rrt(091) * density(15) 
  pd(17,15) = pd(17,15) + rrt(091) * density(09) 
  pd(09,09) = pd(09,09) - rrt(092) * density(09) * 4.0d0
  pd(10,09) = pd(10,09) + rrt(092) * density(09) * 2.0d0
  pd(17,09) = pd(17,09) + rrt(092) * density(09) * 2.0d0
  pd(09,09) = pd(09,09) - rrt(093) * density(11) 
  pd(09,11) = pd(09,11) - rrt(093) * density(09) 
  pd(11,09) = pd(11,09) - rrt(093) * density(11) 
  pd(11,11) = pd(11,11) - rrt(093) * density(09) 
  pd(17,09) = pd(17,09) + rrt(093) * density(11) * 2.0d0
  pd(17,11) = pd(17,11) + rrt(093) * density(09) * 2.0d0
  pd(18,09) = pd(18,09) + rrt(093) * density(11) 
  pd(18,11) = pd(18,11) + rrt(093) * density(09) 
  pd(09,10) = pd(09,10) + rrt(094) * density(36) 
  pd(09,36) = pd(09,36) + rrt(094) * density(10) 
  pd(10,10) = pd(10,10) - rrt(094) * density(36) 
  pd(10,36) = pd(10,36) - rrt(094) * density(10) 
  pd(18,10) = pd(18,10) + rrt(094) * density(36) 
  pd(18,36) = pd(18,36) + rrt(094) * density(10) 
  pd(36,10) = pd(36,10) - rrt(094) * density(36) 
  pd(36,36) = pd(36,36) - rrt(094) * density(10) 
  pd(02,10) = pd(02,10) + rrt(095) * density(15) 
  pd(02,15) = pd(02,15) + rrt(095) * density(10) 
  pd(10,10) = pd(10,10) - rrt(095) * density(15) 
  pd(10,15) = pd(10,15) - rrt(095) * density(10) 
  pd(15,10) = pd(15,10) - rrt(095) * density(15) 
  pd(15,15) = pd(15,15) - rrt(095) * density(10) 
  pd(18,10) = pd(18,10) + rrt(095) * density(15) 
  pd(18,15) = pd(18,15) + rrt(095) * density(10) 
  pd(09,10) = pd(09,10) + rrt(096) * density(15) 
  pd(09,15) = pd(09,15) + rrt(096) * density(10) 
  pd(10,10) = pd(10,10) - rrt(096) * density(15) 
  pd(10,15) = pd(10,15) - rrt(096) * density(10) 
  pd(15,10) = pd(15,10) - rrt(096) * density(15) 
  pd(15,15) = pd(15,15) - rrt(096) * density(10) 
  pd(16,10) = pd(16,10) + rrt(096) * density(15) 
  pd(16,15) = pd(16,15) + rrt(096) * density(10) 
  pd(11,11) = pd(11,11) - rrt(097) * density(36) 
  pd(11,36) = pd(11,36) - rrt(097) * density(11) 
  pd(17,11) = pd(17,11) + rrt(097) * density(36) * 2.0d0
  pd(17,36) = pd(17,36) + rrt(097) * density(11) * 2.0d0
  pd(36,11) = pd(36,11) - rrt(097) * density(36) 
  pd(36,36) = pd(36,36) - rrt(097) * density(11) 
  pd(15,16) = pd(15,16) + rrt(098) * density(36) 
  pd(15,36) = pd(15,36) + rrt(098) * density(16) 
  pd(16,16) = pd(16,16) - rrt(098) * density(36) 
  pd(16,36) = pd(16,36) - rrt(098) * density(16) 
  pd(18,16) = pd(18,16) + rrt(098) * density(36) 
  pd(18,36) = pd(18,36) + rrt(098) * density(16) 
  pd(36,16) = pd(36,16) - rrt(098) * density(36) 
  pd(36,36) = pd(36,36) - rrt(098) * density(16) 
  pd(09,02) = pd(09,02) + rrt(099) * density(15) * density(17) 
  pd(09,15) = pd(09,15) + rrt(099) * density(02) * density(17) 
  pd(09,17) = pd(09,17) + rrt(099) * density(02) * density(15) 
  pd(15,02) = pd(15,02) - rrt(099) * density(15) * density(17) 
  pd(15,15) = pd(15,15) - rrt(099) * density(02) * density(17) 
  pd(15,17) = pd(15,17) - rrt(099) * density(02) * density(15) 
  pd(17,02) = pd(17,02) - rrt(099) * density(15) * density(17) 
  pd(17,15) = pd(17,15) - rrt(099) * density(02) * density(17) 
  pd(17,17) = pd(17,17) - rrt(099) * density(02) * density(15) 
  pd(09,15) = pd(09,15) + rrt(100) * density(17) 
  pd(09,17) = pd(09,17) + rrt(100) * density(15) 
  pd(15,15) = pd(15,15) - rrt(100) * density(17) 
  pd(15,17) = pd(15,17) - rrt(100) * density(15) 
  pd(17,15) = pd(17,15) - rrt(100) * density(17) 
  pd(17,17) = pd(17,17) - rrt(100) * density(15) 
  pd(11,11) = pd(11,11) - rrt(101) * density(15) 
  pd(11,15) = pd(11,15) - rrt(101) * density(11) 
  pd(15,11) = pd(15,11) - rrt(101) * density(15) 
  pd(15,15) = pd(15,15) - rrt(101) * density(11) 
  pd(17,11) = pd(17,11) + rrt(101) * density(15) 
  pd(17,15) = pd(17,15) + rrt(101) * density(11) 
  pd(18,11) = pd(18,11) + rrt(101) * density(15) 
  pd(18,15) = pd(18,15) + rrt(101) * density(11) 
  pd(30,02) = pd(30,02) - rrt(102) * density(30) 
  pd(30,30) = pd(30,30) - rrt(102) * density(02) 
  pd(36,02) = pd(36,02) + rrt(102) * density(30) 
  pd(36,30) = pd(36,30) + rrt(102) * density(02) 
  pd(02,02) = pd(02,02) - rrt(103) * density(30) 
  pd(02,30) = pd(02,30) - rrt(103) * density(02) 
  pd(18,02) = pd(18,02) + rrt(103) * density(30) * 2.0d0
  pd(18,30) = pd(18,30) + rrt(103) * density(02) * 2.0d0
  pd(30,02) = pd(30,02) - rrt(103) * density(30) 
  pd(30,30) = pd(30,30) - rrt(103) * density(02) 
  pd(02,02) = pd(02,02) - rrt(104) * density(31) 
  pd(02,31) = pd(02,31) - rrt(104) * density(02) 
  pd(18,02) = pd(18,02) + rrt(104) * density(31) * 2.0d0
  pd(18,31) = pd(18,31) + rrt(104) * density(02) * 2.0d0
  pd(31,02) = pd(31,02) - rrt(104) * density(31) 
  pd(31,31) = pd(31,31) - rrt(104) * density(02) 
  pd(30,30) = pd(30,30) - rrt(105) 
  pd(36,30) = pd(36,30) + rrt(105) 
  pd(32,34) = pd(32,34) + rrt(106) 
  pd(34,34) = pd(34,34) - rrt(106) 
  pd(33,35) = pd(33,35) + rrt(107) 
  pd(35,35) = pd(35,35) - rrt(107) 
  pd(31,32) = pd(31,32) + rrt(108) 
  pd(32,32) = pd(32,32) - rrt(108) 
  pd(31,33) = pd(31,33) + rrt(109) 
  pd(33,33) = pd(33,33) - rrt(109) 
  pd(18,19) = pd(18,19) + rrt(110) 
  pd(19,19) = pd(19,19) - rrt(110) 
  pd(12,12) = pd(12,12) - rrt(111) 
  pd(15,12) = pd(15,12) + rrt(111) 
  pd(12,13) = pd(12,13) + rrt(112) 
  pd(13,13) = pd(13,13) - rrt(112) 
  pd(12,14) = pd(12,14) + rrt(113) 
  pd(14,14) = pd(14,14) - rrt(113) 
  pd(02,03) = pd(02,03) + rrt(114) 
  pd(03,03) = pd(03,03) - rrt(114) 
  pd(02,04) = pd(02,04) + rrt(115) 
  pd(04,04) = pd(04,04) - rrt(115) 
  pd(02,05) = pd(02,05) + rrt(116) 
  pd(05,05) = pd(05,05) - rrt(116) 
  pd(02,06) = pd(02,06) + rrt(117) 
  pd(06,06) = pd(06,06) - rrt(117) 
  pd(02,07) = pd(02,07) + rrt(118) 
  pd(07,07) = pd(07,07) - rrt(118) 
  pd(02,08) = pd(02,08) + rrt(119) 
  pd(08,08) = pd(08,08) - rrt(119) 
  pd(15,15) = pd(15,15) - rrt(120) * density(38) 
  pd(15,38) = pd(15,38) - rrt(120) * density(15) 
  pd(24,15) = pd(24,15) + rrt(120) * density(38) 
  pd(24,38) = pd(24,38) + rrt(120) * density(15) 
  pd(37,15) = pd(37,15) + rrt(120) * density(38) 
  pd(37,38) = pd(37,38) + rrt(120) * density(15) 
  pd(38,15) = pd(38,15) - rrt(120) * density(38) 
  pd(38,38) = pd(38,38) - rrt(120) * density(15) 
  pd(15,22) = pd(15,22) + rrt(121) * density(38) 
  pd(15,38) = pd(15,38) + rrt(121) * density(22) 
  pd(22,22) = pd(22,22) - rrt(121) * density(38) 
  pd(22,38) = pd(22,38) - rrt(121) * density(22) 
  pd(37,22) = pd(37,22) + rrt(121) * density(38) 
  pd(37,38) = pd(37,38) + rrt(121) * density(22) 
  pd(38,22) = pd(38,22) - rrt(121) * density(38) 
  pd(38,38) = pd(38,38) - rrt(121) * density(22) 
  pd(15,16) = pd(15,16) + rrt(122) * density(38) 
  pd(15,38) = pd(15,38) + rrt(122) * density(16) 
  pd(16,16) = pd(16,16) - rrt(122) * density(38) 
  pd(16,38) = pd(16,38) - rrt(122) * density(16) 
  pd(38,16) = pd(38,16) - rrt(122) * density(38) 
  pd(38,38) = pd(38,38) - rrt(122) * density(16) 
  pd(41,16) = pd(41,16) + rrt(122) * density(38) 
  pd(41,38) = pd(41,38) + rrt(122) * density(16) 
  pd(18,20) = pd(18,20) + rrt(123) * density(38) 
  pd(18,38) = pd(18,38) + rrt(123) * density(20) 
  pd(20,20) = pd(20,20) - rrt(123) * density(38) 
  pd(20,38) = pd(20,38) - rrt(123) * density(20) 
  pd(37,20) = pd(37,20) + rrt(123) * density(38) 
  pd(37,38) = pd(37,38) + rrt(123) * density(20) 
  pd(38,20) = pd(38,20) - rrt(123) * density(38) 
  pd(38,38) = pd(38,38) - rrt(123) * density(20) 
  pd(15,20) = pd(15,20) + rrt(124) * density(38) 
  pd(15,38) = pd(15,38) + rrt(124) * density(20) 
  pd(20,20) = pd(20,20) - rrt(124) * density(38) 
  pd(20,38) = pd(20,38) - rrt(124) * density(20) 
  pd(36,20) = pd(36,20) + rrt(124) * density(38) 
  pd(36,38) = pd(36,38) + rrt(124) * density(20) 
  pd(37,20) = pd(37,20) + rrt(124) * density(38) 
  pd(37,38) = pd(37,38) + rrt(124) * density(20) 
  pd(38,20) = pd(38,20) - rrt(124) * density(38) 
  pd(38,38) = pd(38,38) - rrt(124) * density(20) 
  pd(02,02) = pd(02,02) - rrt(125) * density(38) 
  pd(02,38) = pd(02,38) - rrt(125) * density(02) 
  pd(23,02) = pd(23,02) + rrt(125) * density(38) 
  pd(23,38) = pd(23,38) + rrt(125) * density(02) 
  pd(37,02) = pd(37,02) + rrt(125) * density(38) 
  pd(37,38) = pd(37,38) + rrt(125) * density(02) 
  pd(38,02) = pd(38,02) - rrt(125) * density(38) 
  pd(38,38) = pd(38,38) - rrt(125) * density(02) 
  pd(02,02) = pd(02,02) - rrt(126) * density(38) 
  pd(02,38) = pd(02,38) - rrt(126) * density(02) 
  pd(18,02) = pd(18,02) + rrt(126) * density(38) 
  pd(18,38) = pd(18,38) + rrt(126) * density(02) 
  pd(38,02) = pd(38,02) - rrt(126) * density(38) 
  pd(38,38) = pd(38,38) - rrt(126) * density(02) 
  pd(41,02) = pd(41,02) + rrt(126) * density(38) 
  pd(41,38) = pd(41,38) + rrt(126) * density(02) 
  pd(02,23) = pd(02,23) + rrt(127) * density(37) 
  pd(02,37) = pd(02,37) + rrt(127) * density(23) 
  pd(23,23) = pd(23,23) - rrt(127) * density(37) 
  pd(23,37) = pd(23,37) - rrt(127) * density(23) 
  pd(37,23) = pd(37,23) - rrt(127) * density(37) 
  pd(37,37) = pd(37,37) - rrt(127) * density(23) 
  pd(38,23) = pd(38,23) + rrt(127) * density(37) 
  pd(38,37) = pd(38,37) + rrt(127) * density(23) 
  pd(02,28) = pd(02,28) + rrt(128) * density(37) 
  pd(02,37) = pd(02,37) + rrt(128) * density(28) 
  pd(28,28) = pd(28,28) - rrt(128) * density(37) 
  pd(28,37) = pd(28,37) - rrt(128) * density(28) 
  pd(37,28) = pd(37,28) - rrt(128) * density(37) 
  pd(37,37) = pd(37,37) - rrt(128) * density(28) 
  pd(41,28) = pd(41,28) + rrt(128) * density(37) 
  pd(41,37) = pd(41,37) + rrt(128) * density(28) 
  pd(15,21) = pd(15,21) + rrt(129) * density(41) 
  pd(15,41) = pd(15,41) + rrt(129) * density(21) 
  pd(21,21) = pd(21,21) - rrt(129) * density(41) 
  pd(21,41) = pd(21,41) - rrt(129) * density(21) 
  pd(36,21) = pd(36,21) + rrt(129) * density(41) 
  pd(36,41) = pd(36,41) + rrt(129) * density(21) 
  pd(37,21) = pd(37,21) + rrt(129) * density(41) 
  pd(37,41) = pd(37,41) + rrt(129) * density(21) 
  pd(41,21) = pd(41,21) - rrt(129) * density(41) 
  pd(41,41) = pd(41,41) - rrt(129) * density(21) 
  pd(01,01) = pd(01,01) + rrt(131) * density(37) 
  pd(01,37) = pd(01,37) + rrt(131) * density(01) 
  pd(37,01) = pd(37,01) - rrt(131) * density(37) 
  pd(37,37) = pd(37,37) - rrt(131) * density(01) 
  pd(38,01) = pd(38,01) + rrt(131) * density(37) 
  pd(38,37) = pd(38,37) + rrt(131) * density(01) 
  pd(37,01) = pd(37,01) - rrt(132) * density(37) 
  pd(37,37) = pd(37,37) - rrt(132) * density(01) 
  pd(39,01) = pd(39,01) + rrt(132) * density(37) 
  pd(39,37) = pd(39,37) + rrt(132) * density(01) 
  pd(37,01) = pd(37,01) + rrt(133) * density(39) 
  pd(37,39) = pd(37,39) + rrt(133) * density(01) 
  pd(39,01) = pd(39,01) - rrt(133) * density(39) 
  pd(39,39) = pd(39,39) - rrt(133) * density(01) 
  pd(01,01) = pd(01,01) + rrt(134) * density(39) 
  pd(01,39) = pd(01,39) + rrt(134) * density(01) 
  pd(38,01) = pd(38,01) + rrt(134) * density(39) 
  pd(38,39) = pd(38,39) + rrt(134) * density(01) 
  pd(39,01) = pd(39,01) - rrt(134) * density(39) 
  pd(39,39) = pd(39,39) - rrt(134) * density(01) 
  pd(37,37) = pd(37,37) + rrt(135) * density(40) 
  pd(37,40) = pd(37,40) + rrt(135) * density(37) 
  pd(38,37) = pd(38,37) + rrt(135) * density(40) 
  pd(38,40) = pd(38,40) + rrt(135) * density(37) 
  pd(40,37) = pd(40,37) - rrt(135) * density(40) 
  pd(40,40) = pd(40,40) - rrt(135) * density(37) 
  pd(01,39) = pd(01,39) + rrt(136) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(136) * density(39) * 4.0d0
  pd(40,39) = pd(40,39) + rrt(136) * density(39) * 2.0d0
  pd(37,37) = pd(37,37) + rrt(137) * density(37) * density(39) * 2.0d0
  pd(37,39) = pd(37,39) + rrt(137) * density(37)**2 
  pd(39,37) = pd(39,37) - rrt(137) * density(37) * density(39) * 2.0d0
  pd(39,39) = pd(39,39) - rrt(137) * density(37)**2 
  pd(38,38) = pd(38,38) - rrt(138) 
  pd(42,38) = pd(42,38) + rrt(138) 
  pd(40,40) = pd(40,40) - rrt(139) 
  pd(43,40) = pd(43,40) + rrt(139) 
  pd(23,23) = pd(23,23) - rrt(140) 
  pd(45,23) = pd(45,23) + rrt(140) 
  pd(25,25) = pd(25,25) - rrt(141) 
  pd(46,25) = pd(46,25) + rrt(141) 
  pd(22,22) = pd(22,22) - rrt(142) 
  pd(47,22) = pd(47,22) + rrt(142) 
  pd(01,01) = pd(01,01) - rrt(143) 
  pd(44,01) = pd(44,01) + rrt(143) 
  if( ldensity_constant ) then
    do i = 1, species_max
      if( density_constant(i) ) pd(i,:) = 0.0d0
    enddo
  endif
  if( lgas_heating ) then
    pd(48,1) = eV_to_K * ZDPlasKin_cfg(11)
    pd(48,:) = pd(48,:) * ZDPlasKin_cfg(13)
  endif
  return
end subroutine ZDPlasKin_jex
end module ZDPlasKin
!-----------------------------------------------------------------------------------------------------------------------------------
!
! reaction constant rates
!
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ZDPlasKin_reac_rates(Time)
  use ZDPlasKin, only : ZDPlasKin_bolsig_rates, bolsig_rates, bolsig_pointer, ZDPlasKin_cfg, ZDPlasKin_get_density_total, &
                        lreaction_block, rrt
  USE OPTIONS, ONLY : GAS_PRESSURE, RADIUS, GAP_LENGTH
  implicit none
  double precision, intent(in) :: Time
  double precision :: Tgas
  double precision :: EN
  double precision :: Te
  DOUBLE PRECISION :: DIFF_RATE
  call ZDPlasKin_bolsig_rates()
  Tgas = ZDPlasKin_cfg(1)
  EN  = ZDPlasKin_cfg(3)
  Te  = ZDPlasKin_cfg(4)
  DIFF_RATE = 1.52D0 * ( 760.D0 / GAS_PRESSURE) * ( TGAS / 273.16D0 ) * ( TE / 11600.D0 ) &
  * ( (2.405D0/RADIUS)**2 + (3.141D0/GAP_LENGTH)**2 )
  rrt(001) = bolsig_rates(bolsig_pointer(1))
  rrt(002) = bolsig_rates(bolsig_pointer(2))
  rrt(003) = bolsig_rates(bolsig_pointer(3))
  rrt(004) = bolsig_rates(bolsig_pointer(4))
  rrt(005) = bolsig_rates(bolsig_pointer(5))
  rrt(006) = bolsig_rates(bolsig_pointer(6))
  rrt(007) = bolsig_rates(bolsig_pointer(7))
  rrt(008) = bolsig_rates(bolsig_pointer(8))
  rrt(009) = bolsig_rates(bolsig_pointer(9))
  rrt(010) = bolsig_rates(bolsig_pointer(10))
  rrt(011) = bolsig_rates(bolsig_pointer(11))
  rrt(012) = bolsig_rates(bolsig_pointer(12))
  rrt(013) = bolsig_rates(bolsig_pointer(13))
  rrt(014) = bolsig_rates(bolsig_pointer(14))
  rrt(015) = bolsig_rates(bolsig_pointer(15))
  rrt(016) = bolsig_rates(bolsig_pointer(16))
  rrt(017) = bolsig_rates(bolsig_pointer(17))
  rrt(018) = bolsig_rates(bolsig_pointer(18))
  rrt(019) = bolsig_rates(bolsig_pointer(19))
  rrt(020) = bolsig_rates(bolsig_pointer(20))
  rrt(021) = bolsig_rates(bolsig_pointer(21))
  rrt(022) = bolsig_rates(bolsig_pointer(22))
  rrt(023) = bolsig_rates(bolsig_pointer(23))
  rrt(024) = 3.8D-15
  rrt(025) = 1.7D-15
  rrt(026) = 1.4D-15
  rrt(027) = 1.5D-16
  rrt(028) = 3.0D-15
  rrt(029) = 4.9D-16*EXP(0.36D0/TGAS)
  rrt(030) = 1.6D-16*EXP(0.44D0/TGAS)
  rrt(031) = 3.0D-15
  rrt(032) = 3.6D-15
  rrt(033) = 3.4D-15
  rrt(034) = 2.33D-15
  rrt(035) = 2.0D-16
  rrt(036) = 1.5D-16
  rrt(037) = 1.5D-16
  rrt(038) = 1.5D-16
  rrt(039) = 5.0D-21*EXP(2.692D0-0.07D0/TGAS)
  rrt(040) = 6.0D-16
  rrt(041) = 1.5D-16
  rrt(042) = 1.2D-15
  rrt(043) = 1.0D-15
  rrt(044) = 1.5D-15
  rrt(045) = 1.5D-15
  rrt(046) = 1.5D-15
  rrt(047) = 0.66D-13*(0.01D0/TE)
  rrt(048) = 0.3D-13*(0.01D0/TE)
  rrt(049) = 0.24D-13*(0.01D0/TE)
  rrt(050) = 5.0D-15/SQRT(TE)
  rrt(051) = 0.6D-14/SQRT(TE)
  rrt(052) = 1.2D-12*(0.026D0/TE)
  rrt(053) = 2.0D-14*SQRT(0.026D0/TE)
  rrt(054) = 2.7D-19*(TE)**(-0.7D0)
  rrt(055) = 4.0D-39*(TE)**(-4.5D0)
  rrt(056) = rrt(55)
  rrt(057) = 7.00D-15
  rrt(058) = 1.30D-13
  rrt(059) = 0.67D-13
  rrt(060) = 3.09D-13
  rrt(061) = 3.09D-13
  rrt(062) = 3.08D-13
  rrt(063) = 3.00D-13
  rrt(064) = 1.00D-13
  rrt(065) = 3.09D-13
  rrt(066) = 2.25D-13
  rrt(067) = 1.05D-13
  rrt(068) = 1.03D-13
  rrt(069) = 1.09D-13
  rrt(070) = 3.09D-13
  rrt(071) = 2.25D-13
  rrt(072) = 1.05D-13
  rrt(073) = 1.02D-13
  rrt(074) = 3.08D-13
  rrt(075) = 1.02D-13
  rrt(076) = 1.01D-13
  rrt(077) = 3.05D-13
  rrt(078) = 0.90D-13
  rrt(079) = 1.80D-18
  rrt(080) = 7.7D-18*EXP(-0.1806D0/TGAS)
  rrt(081) = 2.3D-17*EXP(-0.00946D0/TGAS)
  rrt(082) = 2.1D-16
  rrt(083) = 1.03D-22
  rrt(084) = 1.7D-18
  rrt(085) = 1.51D-17
  rrt(086) = 1.9D-18*EXP(-0.086D0/TGAS)
  rrt(087) = 8.0D-17
  rrt(088) = 2.0D-18
  rrt(089) = 4.8D-17*EXP(-0.0215D0/TGAS)
  rrt(090) = 3.0D-17*EXP(-0.0172D0/TGAS)
  rrt(091) = 3.08D-18
  rrt(092) = 3.01D-18
  rrt(093) = 1.4D-20*EXP(-0.0516D0/TGAS)
  rrt(094) = 1.8D-21
  rrt(095) = 5.09D-20
  rrt(096) = 1.22D-22
  rrt(097) = 2.0D-17*EXP(-0.196/TGAS)
  rrt(098) = 9.1D-24
  rrt(099) = 6.4D-43
  rrt(100) = 1.8D-18
  rrt(101) = 1.4D-16*EXP(-0.04128D0/TGAS)
  rrt(102) = 3.6D-17
  rrt(103) = 2.5D-16
  rrt(104) = 2.5D-16
  rrt(105) = 0.7D0
  rrt(106) = 0.322D8
  rrt(107) = 0.369D8
  rrt(108) = 1.879D8
  rrt(109) = 1.8D-4
  rrt(110) = 800D-9
  rrt(111) = 1.6D-9
  rrt(112) = 10.2D-9
  rrt(113) = 33.7D-9
  rrt(114) = 1.25D6
  rrt(115) = 1.25D6
  rrt(116) = 1.25D6
  rrt(117) = 1.25D6
  rrt(118) = 1.25D6
  rrt(119) = 1.25D6
  rrt(120) = 1.0D-10
  rrt(121) = 2.0D-7
  rrt(122) = 1.1D-9
  rrt(123) = 2.0D-7
  rrt(124) = 1.0D-7
  rrt(125) = 7.0D-10
  rrt(126) = 3.0D-10
  rrt(127) = 2.2D-10
  rrt(128) = 1.0D-11
  rrt(129) = 1.0D-7
  rrt(130) = bolsig_rates(bolsig_pointer(24))
  rrt(131) = bolsig_rates(bolsig_pointer(25))
  rrt(132) = bolsig_rates(bolsig_pointer(26))
  rrt(133) = bolsig_rates(bolsig_pointer(27))
  rrt(134) = bolsig_rates(bolsig_pointer(28))
  rrt(135) = 6.06D-6/301.0D0*EXP(-15130.0D0/301.0D0)
  rrt(136) = 6.0D-10
  rrt(137) = 1.4D-32
  rrt(138) = DIFF_RATE
  rrt(139) = DIFF_RATE
  rrt(140) = DIFF_RATE
  rrt(141) = DIFF_RATE
  rrt(142) = DIFF_RATE
  rrt(143) = DIFF_RATE
  where( lreaction_block(:) ) rrt(:) = 0.0d0
  return
end subroutine ZDPlasKin_reac_rates
!-----------------------------------------------------------------------------------------------------------------------------------
!
! END OF FILE
!
!-----------------------------------------------------------------------------------------------------------------------------------
