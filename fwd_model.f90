!..................................................................................85
!..                     fwd_model.f90 --- Forward model module                      !
!...................................................................................!
!.. A mostly empty forward model module.   To be modified by user. Currently does   !
!.. bare-bones power law regression. This problem is linear regression when the     !
!.. coefficients ("a parameters") are sampled in log-space and obs have log-errors. !
!..                                                                                 !
!.. A few somewhat important notes:                                                 !
!..   - things with "p" in front of them are for a version of the main program      !
!..     that does not do MCMC and instead samples from the a previously run MCMC    !
!..     chain. These can be safely ignored.                                         !
!..   - The module variables are useful for settings that will be read in and that  !
!..     will apply to all forward simulations                                       !
!..   - The meat of the module sits in run_fwd_model.                               !
!..                                                                                 !
!..                                                                                 !
!..                                                                                 !
!..                                                                                 !
!...................................................................................!
!...................................................................................!
MODULE fwd_model

use shr_kind_mod, only: r8 => shr_kind_r8
use io_util, only: assert, error_abort
use rand_util, only: random_real
!use simple_boss, only: nproc

IMPLICIT NONE

PRIVATE !..Default attribute of private to prevent strange behavior
PUBLIC ::  initialize_nl, initialize_values,  &
           run_fwd_model, load_sigma_obs, &
           load_obs_types, load_sim_obs, load_morr_obs

CHARACTER(LEN=70),PUBLIC :: namelist_flnm, ic_filename
INTEGER, PUBLIC   :: nobs, nsims
!..The following can be used to define native dimensions of the model, which are strung
!..into a vector of length Nobs
INTEGER           :: Npol, Nz, Nt, obs_nslice, t_nslice, nscalarobs, irain
                     !..,nsims
!..Variables below for run_fwd_model
INTEGER, PUBLIC, ALLOCATABLE, DIMENSION(:) :: obs_slices , t_slices!..(obs_nslice)
LOGICAL                 :: obs_in_db, alt_proc, normalized, obs_from_file!, ic_from_clm
CHARACTER(LEN=70)       :: ic_type
character(len=256)      :: obs_file = ""
character(len=256)      :: pobs_file = ""
!REAL                    :: min_obs
REAL, ALLOCATABLE       :: min_obs(:)

!real, allocatable :: nccn(:)
!real, allocatable :: rhinit(:)
!real, allocatable :: lhflux(:)
!real, allocatable :: ampw(:)

CONTAINS
!...................................................................................!
SUBROUTINE initialize_nl

 PRINT*,'Please input namelist filename...'
 READ*,namelist_flnm
 PRINT*,'Posterior namelist filename = ',namelist_flnm

END SUBROUTINE initialize_nl


!...................................................................................!
SUBROUTINE initialize_values(nobs_alt)
INTEGER, PARAMETER :: unit_nl=503, unit_lhs = 132
INTEGER, INTENT(OUT) :: nobs_alt
LOGICAL            :: lhs_flag
real, allocatable :: lhs_base(:,:)

NAMELIST /record0/ t_slices,   obs_slices, obs_in_db, min_obs, &
                   nterms, dt, output_nsteps
                   
NAMELIST /record4/ obs_nslice, t_nslice, nsims, nscalarobs, &
                   obs_from_file, obs_file
NAMELIST /record6/ lhs_flag, alt_proc
 
 OPEN ( unit_nl, FILE=namelist_flnm, FORM='formatted', ACTION='read',    &
        ACCESS='sequential', STATUS='old')
 READ ( unit_nl, NML = record4 )
 READ ( unit_nl, NML = record6 ) 
 CLOSE( unit_nl )

 call assert(.not. (obs_from_file .and. obs_file == ""), &
      "cannot read obs from file if no obs_file was provided")

 !..You can string out obs however you want...
 nobs     = nsims*(obs_nslice*t_nslice*Npol + nscalarobs*t_nslice)  !..vector of len obs_nslice for each Npol
 PRINT*,'NOBS = ',nobs
 
 !..alt_proc gives a vector to output other model variables of interest
 IF (alt_proc) THEN
   !..Process rate output
   !..        nprocesses * n_proc_out_levels * nsims
   nobs_alt =      4     *        3          * nsims
 ELSE
  !..Fully-resolved height level moments output
   nobs_alt = Npol*Nz*nsims
 ENDIF
 !..
 print*,'THE VALUE OF NZ IS',Nz
 print*,'THE VALUE OF NPOL IS',Npol
 PRINT*,'THE VALUE OF OBS_NSLICE IS',obs_nslice
 
 !..This is for if you're using Latin Hypercube to generate initial conditions
 IF (lhs_flag) THEN
   PRINT*,'generating LH sample'
   ALLOCATE( lhs_base(4,Nsims))    !..Here 4=lccn, rh, ff, wamp
   CALL make_lhs_sample(lhs_base)
   call load_from_lhs(lhs_base)
   deallocate(lhs_base)
 ELSE
   PRINT*,'NOT generating LH sample, reading from file instead'
   CALL dont_make_lhs_sample()
   if (ic_type == 'ic_file') then
      ALLOCATE( lhs_base(4,Nsims))    !..Here 4=lccn, rh, ff, wamp
      OPEN( unit_lhs, FILE=ic_filename, action='read')
      READ( unit_lhs,*) lhs_base
      CLOSE( unit_lhs)
      call load_from_lhs(lhs_base)
      deallocate(lhs_base)
   else if (ic_type == 'obs_file') then
      call load_offline_initialization(obs_file)
   else
      call error_abort("unrecognized ic_type: " // trim(ic_type))
   end if
 ENDIF

 !..These need to be allocated before being read from namelist
 ALLOCATE( obs_slices( obs_nslice) )
 ALLOCATE( t_slices( t_nslice) )
 ALLOCATE( min_obs( nscalarobs + Npol ) )

 OPEN ( unit_nl, FILE=namelist_flnm, FORM='formatted', ACTION='read',    &
        ACCESS='sequential', STATUS='old')
 READ ( unit_nl, NML = record0 )
 CLOSE( unit_nl )


contains

  subroutine load_from_lhs(lhs_base)
    real, intent(in) :: lhs_base(:,:)
    nccn = 10.**lhs_base(1,:)
    rhinit = lhs_base(2,:)
    lhflux = lhs_base(3,:)
    ampw = lhs_base(4,:)
  end subroutine load_from_lhs

END SUBROUTINE initialize_values




!...................................................................................!
!..Subroutine to load obs error values
SUBROUTINE load_sigma_obs(sigma_obs, err_from_file)

LOGICAL, INTENT(IN)               :: err_from_file
REAL, DIMENSION(:), INTENT(OUT)   :: sigma_obs
REAL, ALLOCATABLE, DIMENSION(:)   :: sigma_file
CHARACTER(LEN=200)                :: obs_err_file
INTEGER :: i,j,k,l,m
INTEGER, PARAMETER :: unit_nl=504
NAMELIST /record5/ sigma_file, obs_err_file
  
  IF (SIZE(sigma_obs) /= nobs) &
                        STOP 'Something is terribly wrong, size(sigma_obs) /= nobs'
  ALLOCATE( sigma_file(Npol+nscalarobs) )
  !..
  OPEN ( unit_nl, FILE=namelist_flnm, FORM='formatted', ACTION='read',    &
         ACCESS='sequential', STATUS='old')
  READ ( unit_nl, NML = record5 )
  CLOSE( unit_nl )
  !..
  IF (err_from_file) THEN
    STOP 'Not coded yet. Please complain to Marcus'
  ELSE
    k=1
    DO l=1,Nsims
      DO m=1,t_nslice
        DO i=1,Npol
          DO j=1,obs_nslice
            sigma_obs(k) = sigma_file(i)
            k=k+1
          ENDDO
        ENDDO
        DO i=1,nscalarobs
          sigma_obs(k) = sigma_file(Npol+i)
          k=k+1
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  !sigma_obs(nobs) = sigma_file(Npol+1)
  print*,'sigma_obs = ',sigma_obs
  DEALLOCATE(sigma_file)

END SUBROUTINE load_sigma_obs

!...................................................................................!
subroutine load_obs_types(n_obs_types, obs_types)
  !..I'm not sure what this does, leaving it in for now
  integer, intent(out) :: n_obs_types
  integer, intent(out) :: obs_types(:)
  integer :: i, j, k, l, m

  if (size(obs_types) /= nobs) &
       stop 'Something is terribly wrong, size(obs_types) /= nobs'

  n_obs_types = npol + nscalarobs

  k = 1
  do l = 1, Nsims
     do m = 1, t_nslice
        do i = 1, Npol
           do j = 1, obs_nslice
              obs_types(k) = i
              k = k + 1
           end do
        end do
        do i = 1, nscalarobs
           obs_types(k) = Npol + i
           k = k + 1
        end do
     end do
  end do

end subroutine load_obs_types


!...................................................................................!
SUBROUTINE load_sim_obs(obs,obs_alt,nparams,params,nobs_alt)
!..This loads simulated observations from some parameter value (ptrue in namelist)
!..Two options: real_obs=.TRUE.   :: we read real obs from some file
!..             real_obs=.FALSE.  :: we create synthetic obs from "true" params

INTEGER, INTENT(IN)                  :: nparams, nobs_alt
REAL, INTENT(IN), DIMENSION(nparams) :: params
REAL, INTENT(OUT), DIMENSION(:)      :: obs, obs_alt
!..
  !.
  PRINT*, 'Sim. Obs parameters = ',params
  !..
  IF (SIZE(obs) /= nobs) STOP 'Something is terribly wrong, size(obs) /= nobs'
  IF (SIZE(obs_alt) /= nobs_alt) STOP 'something is.... size(obs_alt) /= nobs_alt'
  !..
  CALL run_fwd_model(params, obs, obs_alt)
  
END SUBROUTINE load_sim_obs


!...................................................................................!
!..Even if you don't use LHS, things still need to be set
SUBROUTINE dont_make_lhs_sample
IMPLICIT NONE

INTEGER, PARAMETER                :: unit_nl = 506
REAL, DIMENSION(2)                :: qtop_sim, ntop_sim, mutop_sim, ssub

NAMELIST /record1/ qtop_sim,  ntop_sim,  mutop_sim, ssub, &
                   ic_type

  OPEN ( unit_nl, FILE=namelist_flnm, FORM='formatted', ACTION='read',    &
         ACCESS='sequential', STATUS='old')
  READ ( unit_nl, NML = record1 )
  CLOSE( unit_nl )

  call assert(.not. obs_from_file .or. ic_type == 'obs_file', &
       "if obs are taken from file, should use ic_type = obs_file")

  call assert(ic_type /= 'obs_file' .or. obs_file /= "", &
       "cannot use ic_type = obs_file if no obs_file is provided")

END SUBROUTINE dont_make_lhs_sample

!...................................................................................!
SUBROUTINE make_lhs_sample(lhs_base)
!..leaving this in to give an example of how to use LHS to set IC's
use rand_util, only: latin_random

REAL, DIMENSION(:,:),INTENT(OUT)  :: lhs_base

INTEGER, PARAMETER                :: unit_lhs=234, &
                                     unit_nl = 505,&
                                     unit_rand = 324, &
                                     kunit = 654
!INTEGER, ALLOCATABLE, DIMENSION(:):: iseed
REAL, DIMENSION(2)                :: lccn_sim,  rh_sim, ff_sim, wamp_sim
INTEGER                           :: i
LOGICAL                           :: is_okay

NAMELIST /record1/ lccn_sim,  rh_sim, ff_sim, wamp_sim,   &
                   ic_type
  
  !..There is no need to initialize random number generator because that is now
  !..done in the main program with `kissinit'. That should be sufficient for all
  !..other random number calls. Note that this means the LHS and MCMC use the 
  !..same random number sequence.
  
  !! Initialize random number generator
  !OPEN(unit_rand, FILE='random_seed', STATUS='old', FORM='formatted', ACTION='read', &
  !     ACCESS='sequential')
  !READ  ( unit_rand, * ) k
  !ALLOCATE ( iseed(k) )
  !READ  ( unit_rand, * ) iseed
  !CLOSE ( unit_rand )

  !seed = iseed(1) + (seed_num*100000)

  OPEN ( unit_nl, FILE=namelist_flnm, FORM='formatted', ACTION='read',    &
         ACCESS='sequential', STATUS='old')
  READ ( unit_nl, NML = record1 )
  CLOSE( unit_nl )

  IF (ic_type == 'lhs') THEN
    PRINT*, SHAPE(lhs_base)
    lhs_base=0.
    !..Run latin hypercube sampling to generate samples in CCN, RH, Flux Forcing, w_amp
    !CALL latin_random( 3, Nsims, seed, lhs_base)   !..seed is now ignored
    CALL latin_random( 4, Nsims, lhs_base)   !..seed is now ignored
    print*,'lhs base first'
    print*,lhs_base
    lhs_base(1,:) = (lccn_sim(2)-lccn_sim(1))*lhs_base(1,:) + lccn_sim(1)
    lhs_base(2,:) = (rh_sim(2)-rh_sim(1))*lhs_base(2,:) + rh_sim(1)
    lhs_base(3,:) = (ff_sim(2)-ff_sim(1))*lhs_base(3,:) + ff_sim(1)
    lhs_base(4,:) = (wamp_sim(2)-wamp_sim(1))*lhs_base(4,:) + wamp_sim(1)
    print*,'lhs base second'
    print*,lhs_base
  ELSE
    call error_abort("make_lhs_sample not set up for ic_type = " // trim(ic_type))
  ENDIF

  !..ADD IN ERROR HANDLING HERE -- MAKE SURE MORR DOESN'T BORK!!
  !..IF MORR borks, then stop right there
  IF (TRIM(ic_type)=='lhs') THEN
    DO i=1,nsims
      CALL test_run_morr(lhs_base(:,i), is_okay)
      IF (.NOT. (is_okay)) THEN
        print*,lhs_base(:,i)
        PRINT*,'BAD IC ALERT!'
        STOP 'IC resulted in crash for Morr'
      ENDIF
    ENDDO
  ENDIF

  !..Write out to some file
  OPEN( unit_lhs, file=ic_filename,  action='write')
  WRITE( unit_lhs,*) lhs_base
  CLOSE( unit_lhs)

END SUBROUTINE make_lhs_sample

!...................................................................................!
SUBROUTINE run_fwd_model(params, fwd_obs, fwd_obs_alt)


!..
INTEGER, PARAMETER                    :: unit_nl = 203!, &
                                         !unit_lhs= 132
REAL,DIMENSION(:),INTENT(IN)          :: params   !..Parameter values
REAL,DIMENSION(:),INTENT(OUT)         :: fwd_obs, fwd_obs_alt
REAL                                  :: pi
                                         
!INTEGER, ALLOCATABLE, DIMENSION(:)    :: obs_slices
INTEGER                               :: i,j,k,dm,dma,eflag(Nsims),ip,im

integer :: auto_ignore


  ip = 1


  pi = 4.*atan(1.)
  
  !..Cacalate bbreak1 from c parameter
  !bbreak1 = params(14)+(LOG(10.**params(13))-LOG(10.**params(15)))/LOG(params(16)**3*pi*1000.)

  fwd_obs = 0.
  dm  = 0
  dma = 1
 
  !...ADD IN CODE HERE TO RUN FWD MODEL FROM PARAMS
  !...ADD IN CODE HERE TO RUN FWD MODEL FROM PARAMS
  !...ADD IN CODE HERE TO RUN FWD MODEL FROM PARAMS
  !...ADD IN CODE HERE TO RUN FWD MODEL FROM PARAMS
  !...ADD IN CODE HERE TO RUN FWD MODEL FROM PARAMS


END SUBROUTINE run_fwd_model




end module fwd_model
