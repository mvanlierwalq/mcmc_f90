!----------------------------------------------------------------------------------85
!                           MCMC ~ ~ ~ MAIN.F90                                     !
!-----------------------------------------------------------------------------------!
!..Lightly edited version of MCMC code used in mcmc ice sticking work, now being used
!..for MCMC constraint of 1D rain shaft model with BOSS           microphysics      !
!..  start edits: Nov 2 2015                                                        !
!..                                                                                 !
!..UPDATE 12/01/2017 -- specific to 2-term BOSS                                     !
!..  Owing to the possible need to sample both the positive and negative ranges of  !
!..  a-parameter space while ALSO sampling in log- (rather than linear) space of    !
!..  those parameters, we have to do lots of funny stuff that doesn't exist in      !
!..  the "standard" version of this AM sampler. For example, "p_naught" serves as a !
!..  switch for this funny behavior, and is fed into perturbation and rangecheck    !
!..  routines...                                                                    !
!...................................................................................!


PROGRAM main
use rand_util, only: random_real, kissinit
USE fwd_model
USE netcdf
IMPLICIT NONE
!.. variables
INTEGER, PARAMETER                              ::  &
          unit_nl = 303,        &
          unit_rand = 314,      &
          unit_icov = 315,      &
          unit_mean = 316,      &
          unit_lhs= 152
REAL, PARAMETER                                 ::  &
          target_rate = 0.23,   &   !...Desired (target) acceptance rate
          acc_thresh  = 0.05        !...Convergence threshold (target /pm thresh)
!.......................................
CHARACTER (LEN=200)                             ::  &
          out_file,             &   !..
          sa_file,              &
          prior_icovfile,       &
          prior_meanfile,       &
          cov_file
INTEGER, ALLOCATABLE, DIMENSION(:)              ::  &
          pflag,                &   !..
          b_test,               &
          monotonic,            &
          obs_types
REAL, ALLOCATABLE, DIMENSION(:)                 ::  &
          params,               &   !..
          ptrue,                &   !..
          mptrue,               &   !..
          obs,                  &   !..
          fwd_obs,              &   !..
          fwd_obs_alt,          &   !..
          obs_alt,              &   !..
          fwd_obs_star,         &
          fwd_obs_save,         &
          fwd_obs_alt_star,     &
          fwd_obs_mle,          &
          p_fwd_obs,            &
          p_fwd_obs_alt,        &
          sigma_obs,            &
          p_params,             &   !..Proposal parameter values
          params_mle,           &   !..MLE parameter vector
          dp,                   &   !..Amount to perturb parameters by
          psigma,               &   !..Parameter proposal std. dev. (independent)
          psigma_anneal,        &   !..Param proposal std. dev. for sim. annealing
          psigma_running,       &   !..Running (adaptive) standard dev.
          psigma_old,           &   !..Saved for posterity
          pmean,                &   !..Mean parameter value
          pmean_new,            &   !..New mean parameter value
          pmeen,                &
          pmeen_new,            &
          pstar,                &   !..Intermediate parameter value (delay. rej.)
          pmin,                 &
          pmax,                 &
          p_sa_mean,            &
          p_sa_accmean,         &
          pnaught,              &
          p_linspace_ratio,     &
          mparams,              &
          mparams_mle,          &
          p_mparams,            &
          mpstar,               &
          priormean
DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)     ::  &          
          p_like_save,          &
          p_phi_save,           &
          sigma_mle,            &
          sigma_diag
REAL, ALLOCATABLE, DIMENSION(:,:)               ::  &
          first_guess,          &
          dp_save,              &
          p_store,              &
          icovmat,              &
          !....CHECK THESE v v v
          pcov_cut,             &   !...Cut down to only perturbed parameters
          R,                    &   !...Cholesky decomp of pcov_cut
          pcov,                 &   !...Proposal covariance matrix
          pcov_old,             &   !...Same as above, saved for posterity
          pcov_running,         &   !...Running (adaptive) inv
          pcond                     !...Conditioning matrix (for smooth inversion)
          !data_save,            &   !...Extra saved data array
          !p_data_save               !...proposal saved data array
          !....CHECK THESE ^ ^ ^
!LOGICAL, ALLOCATABLE, DIMENSION(:)              ::  &
!          gauss_prior
INTEGER                                         ::  &
          nmc,                  &   !..
          nguess,               &   !..
          nparams,              &   !..
          mc_ncid,              &   !..
          sa_ncid,              &   !..For simulated annealing
          sb_ncid,              &   !..To get pcov if skip_burn
          param_dimid,          &   !..
          par_dimid,            &   !..For simulated annealing 
          obs_dimid,            &   !..
          obs_alt_dimid,        &   !..
          mc_dimid,             &   !..
          sa_dimid,             &   !..For simulated annealing
          params_varid,         &   !..
          m_params_varid,         &   !..
          par_varid,            &   !..For simulated annealing
          pcov_varid,           &
          sb_pcov_varid,        &   !..For pcov if skip_burn
          obs_varid,            &   !..
          obso_varid,           &   !..
          obs_alt_varid,        &   !..
          obso_alt_varid,       &   !..
          fpsigma_varid,        &   !..
          mle_varid,            &
          mmle_varid,            &
          mptrue_varid,         &
          ll_varid,             &   !..log-likelihood netcdf varid
          lls_varid,            &   !..log-likelihood netcdf varid for sim annealing
          !pcov_varid,           &   !..
          !.........................MCMC VARS
          max_burn,             &   !...Max number of burn iterations = the lesser of
                                    !...this or 0.5*nmc
          nchain,               &   !...Chain number - used for random number stream
          nchain_lhs,           &   !...Chain number - for LHS random number stream
          delay_lev,            &
          delay_count,          &
          max_delay_lev,        &   !...Maximum delay level (delayed rejection > 1)
          cov_mem,              &
          cov_check,            &
          iburn_spin,           &
          cov_count,            &
          n_adapt,              &
          adapt_count,          &
          bcount,               &
          option,               &
          icount,               &
          ios,                  &
          naccept,              &
          nobs_alt,             &
          n_obs_types,          &
          diagnostic_mod,       &
          out_mod,              &   !..How often to write
          !..
          success,              &
          nprtb,                &   !..Number of perturbed parameters
          nprtb_tmp,            &   !..Number of 
          i,j,k,p,              &
          sim_ann_type,         &   !..Controls if/how simulated annealing pre-samplr
          n_sa                      !..number of simulated annealing samples
REAL(8)                                         ::  urand1
REAL                                            ::  &
          urand,                &   !...Uniform random number
          nrand1,               &   !...First 0-mean, unit std. dev rand number
          nrand2,               &   !...Second 0-mean, unit std dev rand number
          s_n,                  &   !...Factor to inflate gaussian std dev
          s_n_sq,               &   !...Factor to inflate gaussian covariance
          acc_rate,             &   !...Current acceptance rate
          inv_adapt_count,      &   !...1/adapt_count
          init_par_var,         &   !..Initial proposal variance scaling
          init_par_var_sa,      &   !..Initial proposal variance scaling sim anneal
          scale_s_n,            &   !..Scale factor for MCMC perturbation scaling
          obs_err_mult,         &   !..Scale factor to inflate obs error
          invT_scale,           &   !..Temp scaling factor for sim annealing
          invT,                 &
          prior_phi                 !..log likelihood for the prior distribution
DOUBLE PRECISION                                ::  &
          !...cost function vars:
          phi,                  &   !...Cost function value (log-likelihood)
          phi_save,             &
          p_phi,                &   !...Proposal cost function
          phi_mle,              &   !...Maximum log-likelihood
          phi_star,             &   !...Intermediate cost function
          like,                 &   !...Accepted (prior) likelihood
          p_like,               &   !...Proposal likelihood
          prop_fac,             &
          like_mle,             &   !...Maximum likelihood estimate
          like_star,            &   !...Intermediate likelihood (adaptive)
          alpha_star                !...Intermediate factor... see Haario et al '06
INTEGER, DIMENSION(2)                           ::  &
          count_params,         &   !..
          count_obs,            &   !..
          count_obs_alt,        &   !..
          start                     !..
LOGICAL                                         ::  &
          set_guess,            &   !..Use prescribed or random first guess
          sim_obs,              &   !..External 'obs' from truth params or real obs
          morr_obs,             &   !..Use simulated obs from MORR bulk scheme sim
          cov_prop,             &   !..Consider correlated proposal perturbations
          uniform_dist,         &   !..Sample parameter space uniformly
          skip_burn,            &   !..Skip burn-in procedure, move to main or adapt
          adaptive,             &   !..Use adaptive Metropolis algorithm
          dr_star,              &   !..Switch for choice of p_star
          found_pert,           &   !
          main_chain,           &   !..
          accepted,             &   !..
          delaying,             &   !..
          burn_in,              &   !..
          use_alt_obs,          &   !..
          !monotonic,            &   !..Whether to demand monotonicity in parameters
          err_from_file,        &   !..Err from file or from namelist?
          logdum,               &
          save_fwd_obs,         &   !..Whether to save fwd_obs from each sample
          save_obs_alt,         &   !..Whether to save obs_alt from each sample
          print_diag,           &   !..Whether to print diagnostics at this iteration
          scam_burn,            &   !..Whether to allow for covariance in burn-in
          sim_ann_run,          &   !..Whether we are within simulated annealing loop
          gauss_prior,          &   !..Whether to use a Gaussian prior for parameters.
          sigma_unknown           !..Whether to diagnose sigmas instead of reading them in.

! We could just declare this with "external :: dpotf2".
! But an explicit interface is better for compiler error checking.
interface
   subroutine dpotf2(uplo, n, a, lda, info)
     use shr_kind_mod, only: r8 => shr_kind_r8
     character(len=1) :: uplo
     integer :: n
     integer :: lda
     real(r8), dimension(lda,*) :: a
     integer :: info
   end subroutine dpotf2
end interface

!.......................... ~ ~ ~END DECLARATIONS ~ ~ ~ ............................!

!......................................................................!
!                         Initial Setup                                ! 
!......................................................................!


!..Primary namelist statement, open and read-in
NAMELIST /record2/ nchain, nchain_lhs, nmc, set_guess, nguess, nparams, out_file,   &
                   sim_obs, morr_obs,  cov_check, cov_mem, max_burn, scale_s_n,     &
                   init_par_var,  &
                   uniform_dist, cov_prop, n_adapt, adaptive, max_delay_lev,        &
                   skip_burn, dr_star, obs_err_mult, use_alt_obs,         &
                   err_from_file, save_fwd_obs, save_obs_alt,             &
                   diagnostic_mod, out_mod, ic_filename, scam_burn,       &
                   sim_ann_type, invT_scale, n_sa, sa_file, init_par_var_sa, &
                   cov_file, sigma_unknown
NAMELIST /record3/ ptrue, pmin, pmax, pflag, first_guess, monotonic, &
                   pnaught, p_linspace_ratio, gauss_prior, prior_icovfile, prior_meanfile
!..

PRINT*, 'STARTING MARCUS MAIN PROGRAM!!'

CALL initialize_nl

call initialize_wv_sat

!..
!OPEN  ( unit_nl, FILE='namelist.input', FORM='formatted', ACTION='read',             &
OPEN  ( unit_nl, FILE=namelist_flnm, FORM='formatted', ACTION='read',             &
       ACCESS='sequential', STATUS='old')
!..Read in MCMC parameters
READ  ( unit_nl, NML= record2 )
CLOSE ( unit_nl )


!IF (cov_prop) STOP 'Code does not allow this yet!!'
IF (max_delay_lev .GT. 2) STOP 'Maximum delay level is 2!!!'

!..Get value of nobs, however that is done. Note that fwd_model may require that
!..this is run before anything else
! Initialize random number generator
PRINT*,'LHS Random number seed = ',nchain_lhs,' x 100'
CALL kissinit(nchain_lhs*100)
DO i = 1,10
  urand1 = random_real()
  print*, 'random number ',i,' = ',urand1
ENDDO

!............................
CALL initialize_values(nobs_alt)
print*,'nobs     = ',nobs
print*,'nobs_alt = ',nobs_alt
!............................
! Initialize random number generator
PRINT*,'MCMC Random number seed = ',nchain,' x 100'
CALL kissinit(nchain*100)
DO i = 1,10
  urand1 = random_real()
  print*, 'random number ',i,' = ',urand1
ENDDO
!! Initialize random number generator
!OPEN(unit_rand, FILE='random_seed', STATUS='old', FORM='formatted', ACTION='read', &
!     ACCESS='sequential')
!READ  ( unit_rand, * ) k
!ALLOCATE ( iseed(k) )
!READ  ( unit_rand, * ) iseed
!CLOSE ( unit_rand )
!!..Intialize random number seed
!iseed = iseed + (nchain*100000)
!CALL RANDOM_SEED(put=iseed)
!............................

ALLOCATE( ptrue(nparams) )
ALLOCATE( mptrue(nparams) )
ALLOCATE( pmin(nparams) )
ALLOCATE( pmax(nparams) )
ALLOCATE( pflag(nparams) )
ALLOCATE( pmean(nparams) )
ALLOCATE( pmean_new(nparams) )
ALLOCATE( pmeen(nparams) )
ALLOCATE( pmeen_new(nparams) )
ALLOCATE( pstar(nparams) )
ALLOCATE( params(nparams) ) 
ALLOCATE( p_params(nparams) )
ALLOCATE( params_mle(nparams) )
ALLOCATE( dp(nparams))
ALLOCATE( psigma(nparams))
ALLOCATE( psigma_anneal(nparams) )
ALLOCATE( psigma_running(nparams) )
ALLOCATE( psigma_old(nparams) )
ALLOCATE( monotonic(nparams) )
ALLOCATE( p_sa_mean(nparams) )
ALLOCATE( p_sa_accmean(nparams) )
!..
ALLOCATE( pnaught(nparams) )
ALLOCATE( p_linspace_ratio(nparams) )
ALLOCATE( mparams(nparams) )
ALLOCATE( mparams_mle(nparams) )
ALLOCATE( p_mparams(nparams) )
ALLOCATE( mpstar(nparams) )
!..
ALLOCATE( obs(nobs) )
ALLOCATE( fwd_obs(nobs) )
ALLOCATE( fwd_obs_star(nobs) )
ALLOCATE( fwd_obs_save(nobs) ) 
ALLOCATE( obs_alt(nobs_alt) )
ALLOCATE( fwd_obs_alt(nobs_alt) )
ALLOCATE( fwd_obs_alt_star(nobs_alt) )
ALLOCATE( fwd_obs_mle(nobs) )
ALLOCATE( p_fwd_obs(nobs) )
ALLOCATE( p_fwd_obs_alt(nobs_alt) )
ALLOCATE( sigma_obs(nobs) )
ALLOCATE( obs_types(nobs) )
ALLOCATE( b_test(cov_mem) )
!..
ALLOCATE( pcov(nparams,nparams) )
ALLOCATE( pcov_old(nparams,nparams) )
ALLOCATE( pcov_running(nparams,nparams) )
ALLOCATE( pcond(nparams,nparams) )
!..
ALLOCATE( dp_save(nparams,max_delay_lev) )
ALLOCATE( p_like_save(max_delay_lev) )
ALLOCATE( p_phi_save(max_delay_lev) ) 
ALLOCATE( p_store(nparams, cov_mem) )


pmax=0.
print*,'pmax = ',pmax

!..Structure of output paramters (nparams,nmc) where nmc is an unlimited dimension 
!..in the netcdf file
!..Structure of output fwd obs (nvd,nzd,nmc) where nmc is an unlimited dimension in 
!..the netcdf file
CALL check(NF90_CREATE( path=out_file, cmode=NF90_CLOBBER, ncid=mc_ncid ))
!..
CALL check(NF90_DEF_DIM( mc_ncid, 'param_dim', nparams, param_dimid))
CALL check(NF90_DEF_DIM( mc_ncid, 'fwd_obs_dim', nobs, obs_dimid))
CALL check(NF90_DEF_DIM( mc_ncid, 'fwd_obs_alt_dim', nobs_alt, obs_alt_dimid))
CALL check(NF90_DEF_DIM( mc_ncid, 'mc_dim', NF90_UNLIMITED, mc_dimid))
!..
CALL check(NF90_DEF_VAR( mc_ncid, 'parameters', NF90_DOUBLE,(/param_dimid, mc_dimid/),  &
                    params_varid))
CALL check(NF90_DEF_VAR( mc_ncid, 'm_parameters', NF90_DOUBLE,(/param_dimid, mc_dimid/),  &
                    m_params_varid))
CALL check(NF90_DEF_VAR( mc_ncid, 'loglikelihood',NF90_DOUBLE, mc_dimid, ll_varid))
CALL check(NF90_DEF_VAR( mc_ncid, 'params_mle', NF90_DOUBLE,(/param_dimid/),mle_varid))
CALL check(NF90_DEF_VAR( mc_ncid, 'm_params_mle', NF90_DOUBLE,(/param_dimid/),mmle_varid))
CALL check(NF90_DEF_VAR( mc_ncid, 'm_params_true', NF90_DOUBLE,(/param_dimid/),mptrue_varid))
IF (save_fwd_obs) &
  CALL check(NF90_DEF_VAR( mc_ncid, 'fwd_obs_unrav',NF90_DOUBLE,(/obs_dimid, mc_dimid/),  &
                      obs_varid))
IF (save_obs_alt) &
  CALL check(NF90_DEF_VAR( mc_ncid, 'fwd_obs_alt_unrav',NF90_DOUBLE,(/obs_alt_dimid,mc_dimid/), &
                    obs_alt_varid))
IF (cov_prop) &
  CALL check(NF90_DEF_VAR( mc_ncid, 'pcov', NF90_DOUBLE, (/param_dimid,param_dimid/),   &
                      pcov_varid))
CALL check(NF90_DEF_VAR( mc_ncid,'obs_unrav',NF90_DOUBLE,(/obs_dimid/),obso_varid))
CALL check(NF90_DEF_VAR(mc_ncid,'obs_alt_unrav',NF90_DOUBLE,(/obs_alt_dimid/),obso_alt_varid))
CALL check(NF90_DEF_VAR( mc_ncid,'final_psigma',NF90_DOUBLE,(/param_dimid/),fpsigma_varid))
CALL check(NF90_ENDDEF( mc_ncid ))
CALL check(NF90_CLOSE(mc_ncid ))


!...................................................................!
!               Load Observations and Obs Error from ...            ! 
!...................................................................!
!..If first guess from file, read in, else generate randomly from prior dist'n
!..Read in namelist record with first guess
ALLOCATE (first_guess(nparams,nguess))
ALLOCATE (priormean(nparams) )
ALLOCATE (icovmat(nparams, nparams))
!ALLOCATE (gauss_prior(nparams))
!OPEN ( unit_nl, FILE='namelist.input', FORM='formatted', ACTION='read',           &
OPEN ( unit_nl, FILE=namelist_flnm, FORM='formatted', ACTION='read',           &
       ACCESS='sequential', STATUS='old')
READ ( unit_nl, NML = record3 )
CLOSE( unit_nl )

print*,'pmax = ',pmax
print*,'first guess =',first_guess

nprtb = SUM(pflag)
ALLOCATE( pcov_cut (nprtb,nprtb) )
ALLOCATE( R        (nprtb,nprtb) )

print*,'monotonic = ',monotonic
!..Load obs from file
option=0
!..Load/set obs error vector or matrix
!CALL load_obs(obs,sigma_obs,option)

!..we assume here that ptrue is actually a model-readable vector of parameter values
mptrue = ptrue
CALL inv_preprocess_parameters(mptrue, pmin, pmax, pnaught, p_linspace_ratio, ptrue)
IF (sim_obs) THEN

  CALL load_sim_obs(obs,obs_alt,nparams,mptrue,nobs_alt)
  CALL load_sigma_obs(sigma_obs, err_from_file)
ELSE
  IF (morr_obs) THEN
    CALL load_morr_obs(obs,obs_alt,nobs_alt)
    CALL load_sigma_obs(sigma_obs, err_from_file)
  ELSE
    STOP "No function to load real obs. FIX plz!"
  ENDIF
ENDIF
if (sigma_unknown) then
   call load_obs_types(n_obs_types, obs_types)
   allocate(sigma_diag(n_obs_types))
   allocate(sigma_mle(n_obs_types))
else
   call load_sigma_obs(sigma_obs, err_from_file)
end if

if (.not. sigma_unknown) then
   sigma_obs = obs_err_mult*sigma_obs
end if
!..Write obs vector to file
print*,'got here 00'
CALL check(NF90_OPEN(path=out_file, mode = NF90_WRITE, ncid=mc_ncid))
print*,'got here 01'
CALL check(NF90_PUT_VAR(mc_ncid,obso_varid,obs))
print*,'got here 02'
CALL check(NF90_PUT_VAR(mc_ncid,obso_alt_varid,obs_alt))
print*,'mptrue = ',mptrue
CALL check(NF90_PUT_VAR(mc_ncid,mptrue_varid,mptrue))
print*,'got here 03'
CALL check(NF90_CLOSE(mc_ncid))
print*,'got here 05'

PRINT*, 'done with obs loading'

!...................................................................!
!                      Run First Guess                              ! 
!...................................................................!

params = ptrue
p_params = ptrue
mparams = mptrue
p_mparams = mptrue

!IF(ANY(pmax<pmin)) STOP 'pmax and pmin set wrong!'

print*,'pflag = ',pflag
print*,'pmin = ',pmin
print*,'pmax = ',pmax
print*,'ptrue = ',ptrue
print*,'first_guess = ',first_guess
print*,'nparams = ',nparams

!..Read in the prior mean and prior inverse covariance matrix. 
!..These arrays are assumed to act in the space of sampled parameters, so
!..in log space if sampling accurs in log space. 
icovmat   = 0.
priormean = 0.
If (gauss_prior) THEN
  OPEN(unit_icov, FILE=prior_icovfile, ACTION='read')
  READ(unit_icov,*) icovmat
  CLOSE(unit_icov)
  print *, "icovmat = ", icovmat
  !..
  OPEN(unit_mean, FILE=prior_meanfile, ACTION='read')
  READ(unit_mean,*) priormean
  CLOSE(unit_mean)
  print *, "priormean = ", priormean
ENDIF

!..Set scale factors, according to Gelman (?)
s_n     = scale_s_n * 2.4 / SQRT(SUM(real(pflag)))
s_n_sq  = s_n**2

!..Run forward model using first_guess (nguess times)
count_params = (/nparams,1/)
count_obs    = (/nobs,1/)
count_obs_alt= (/nobs_alt,1/)
start  = (/1,0/)
IF (set_guess) THEN
  DO i=1,nguess
    print*,'Guess no. ',i
    start(2) = start(2) + 1
    mparams = first_guess(:,i)
    !..Run forward model
    !..No need to run preprocessor on first_guess, but need to run inverse to
    !..initialize value of params
    CALL inv_preprocess_parameters(mparams, pmin, pmax, pnaught, p_linspace_ratio,  params)
    print*,'params first guess  = ',params
    print*,'mparams first guess = ',mparams
    !..
    CALL run_fwd_model(mparams, fwd_obs, fwd_obs_alt)      !..has been pre-processed
    if (sigma_unknown) then
       call likelihood_ind_diag_sigma(obs, fwd_obs, n_obs_types, obs_types, phi, like, sigma_diag)
    else
       IF (use_alt_obs) THEN
          CALL likelihood_ind(obs, fwd_obs_alt, sigma_obs, phi, like)
       ELSE  
          CALL likelihood_ind(obs, fwd_obs, sigma_obs, phi, like)
       ENDIF
    end if
    PRINT*,'TEST. first guess fwd_obs = ',fwd_obs
    IF (save_obs_alt) PRINT*,'TEST. first guess fwd_obs_alt = ',fwd_obs_alt
    phi_save = phi
    PRINT*,'First guess phi = ',phi
    PRINT*,'First guess like = ',like
    if (sigma_unknown) then
       print *, 'First guess sigma_diag = ', sigma_diag
    end if
    !..Write output
    CALL check(NF90_OPEN(path=out_file, mode = NF90_WRITE, ncid=mc_ncid))
    CALL check(NF90_PUT_VAR(mc_ncid,params_varid, params, start=start,count=count_params))
    CALL check(NF90_PUT_VAR(mc_ncid,m_params_varid, mparams, start=start,count=count_params))
    CALL check(NF90_PUT_VAR(mc_ncid,ll_varid, [phi], start=[start(2)],count=[1]))
    IF (save_fwd_obs) CALL check(NF90_PUT_VAR(mc_ncid,obs_varid, fwd_obs, &
                                 start=start,count=count_obs))
    IF (save_obs_alt) CALL check(NF90_PUT_VAR(mc_ncid,obs_alt_varid,fwd_obs_alt, &
                                              start=start,count=count_obs_alt))
    CALL check(NF90_CLOSE(mc_ncid))
  ENDDO
ELSE
  print*,'doing random first guess'
  CALL perturb_uni (pflag,pmin,pmax,pnaught,params)
  CALL preprocess_parameters(params, pmin, pmax, pnaught, p_linspace_ratio, mparams)
  print*,'params first guess  = ',params
  print*,'mparams first guess = ',mparams
  start(2) = start(2) + 1
  CALL run_fwd_model(mparams, fwd_obs, fwd_obs_alt)
  if (sigma_unknown) then
     call likelihood_ind_diag_sigma(obs, fwd_obs, n_obs_types, obs_types, phi, like, sigma_diag)
  else
     IF (use_alt_obs) THEN
        CALL likelihood_ind(obs, fwd_obs_alt, sigma_obs, phi, like)
     ELSE
        CALL likelihood_ind(obs, fwd_obs, sigma_obs, phi, like)
     ENDIF
  end if
  PRINT*,'First guess phi = ',phi
  PRINT*,'First guess like = ',like
  if (sigma_unknown) then
     print *, 'First guess sigma_diag = ', sigma_diag
  end if
  !..Write output
  CALL check(NF90_OPEN(path=out_file, mode = NF90_WRITE, ncid=mc_ncid))
  CALL check(NF90_PUT_VAR(mc_ncid,params_varid, params, start=start,count=count_params))
  CALL check(NF90_PUT_VAR(mc_ncid,m_params_varid, mparams, start=start,count=count_params))
  CALL check(NF90_PUT_VAR(mc_ncid,ll_varid, [phi], start=[start(2)],count=[1]))
  IF (save_fwd_obs) CALL check(NF90_PUT_VAR(mc_ncid,obs_varid, fwd_obs, start=start,count=count_obs))
  IF (save_obs_alt) CALL check(NF90_PUT_VAR(mc_ncid,obs_alt_varid,fwd_obs_alt,start=start,count=count_obs_alt))
  CALL check(NF90_CLOSE(mc_ncid))
ENDIF

prior_phi = -50.
!..Figure out the prior (relative) probability of the first guess
if (gauss_prior) prior_phi = return_prior_phi(params, priormean, icovmat)
!..Note that while this will apply the gaussian prior, it should have a 100% accept probability
IF ( .NOT. (rangecheck(pmin,pmax,pnaught,params,pflag,monotonic, &
                       gauss_prior,priormean,icovmat,prior_phi))) THEN
  PRINT*,'First guess is not in range you dummy!'
  PRINT*,'PMIN                PMAX               PARAMS           MPARAMS'
  DO i=1,nparams
    IF (pnaught(i) < 900.) THEN
      IF(params(i) .GT. 1.0) THEN
        PRINT*,'w',i,0.,1.,params(i),mparams(i)
      ELSEIF(params(i) .LT. 0.0) THEN
        PRINT*,'w',i,0.,1.,params(i),mparams(i)
      ENDIF
    ELSEIF (params(i) .GT. pmax(i)) THEN
      PRINT*,i,pmin(i),pmax(i),params(i),mparams(i)
    ELSEIF (params(i) .LT. pmin(i)) THEN
      PRINT*,i,pmin(i),pmax(i),params(i),mparams(i)
    ENDIF
  ENDDO
  STOP 'you done messed up'
ENDIF
!..Initialize maximum likelihood values
phi_mle     = phi
like_mle    = like
params_mle  = params
mparams_mle = mparams
fwd_obs_mle = fwd_obs
if (sigma_unknown) sigma_mle = sigma_diag

!..Set spinup in burn-in to 10x var checks (make sure burn-in continues for a while)
iburn_spin = cov_check * 10

! Set scale parameter (depends on number of perturbed parameters)
! s_n = scale_s_n * 2.4 / sqrt(float(nparams))
! Bugfix here--make sure we are scaling by actual number of parameters perturbed,
! not the total number of parameters...
s_n     = scale_s_n * 2.4 / SQRT(SUM(real(pflag)))
s_n_sq  = s_n**2
PRINT '(a,f12.3)','Scale factor on Gaussian perturbation: ',s_n
pcond = 0.
pcov_cut = 0.
pcov = 0.
!..Prepare initial proposal variance
IF (cov_prop) THEN
  nprtb = 0
  DO i = 1,nparams
    !...If skip_burn is true, load covariance from file. Otherwise, set.
    IF ( skip_burn ) THEN
      CALL check(NF90_OPEN(path=cov_file, mode=NF90_NOWRITE, ncid = sb_ncid) )
      CALL check(NF90_INQ_VARID(sb_ncid, "pcov", sb_pcov_varid) )
      CALL check(NF90_GET_VAR(sb_ncid, sb_pcov_varid, pcov) )
      CALL check(NF90_CLOSE(sb_ncid) ) 
      !OPEN ( unit_pcov, FILE = TRIM('parameter_cov_04.dat'), STATUS = 'old', &
      !      FORM = 'unformatted', ACTION = 'read', ACCESS = 'stream')
      !READ ( unit_pcov ) pcov
      !CLOSE( unit_pcov )
      !pcov = s_n_sq*pcov
    ELSE
      IF (pnaught(i) > 900.) THEN
        pcov(i,i)   = (init_par_var*( pmax(i) - pmin(i)))**2
        !pcov(i,i)   = init_par_var*( pmax(i) - pmin(i)) 
      ELSE
        pcov(i,i)   = init_par_var**2  !..bc the range of these is 0-1
      ENDIF
    ENDIF

    pcond(i,i)  = 0.01*pcov(i,i)
    !pcond(i,i)  = 1e-10
    !pcond(i,i)  = 1.d-15
    IF ( pflag(i) .EQ. 1) THEN
      nprtb = nprtb + 1
      nprtb_tmp = 0
      DO j = 1,nparams
        IF ( pflag(j) .EQ. 1) THEN
          nprtb_tmp = nprtb_tmp + 1
          pcov_cut(nprtb,nprtb_tmp) = pcov(i,j)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !...Perform Cholesky decomposition, test for success/failure
  R = pcov_cut
  !CALL spotf2('L', nprtb, R, nprtb, success)
  CALL dpotf2('L', nprtb, R, nprtb, success)
  !CALL dpotrf_f95(R, 'U', info=success)
  IF ( success .NE. 0 ) THEN
    STOP 'Cholesky Decomposition Failed'
  ENDIF
  !..Zero non-diagonal elements
  DO i = 1,nprtb
    DO j = 1, i-1
      R(j,i) = 0.d0  !..If dpotf2('L',...)
    ENDDO
  ENDDO
ELSE
  DO i=1,nparams
    IF (pnaught(i) > 900.) THEN
      psigma(i) = init_par_var*( pmax(i) - pmin(i) )
    ELSE
      psigma(i) = init_par_var
    ENDIF
  ENDDO
  psigma = s_n*psigma
  PRINT*,'psigma = ', psigma
ENDIF


!.................................................................................!
!                SIMULATED ANNEALING PRE-SAMPLER (IF REQUESTED)                   !
!.................................................................................!
IF (sim_ann_type /= 0) THEN

  !..Set up output just for Simulated Annealing pre-sampler
  CALL check(NF90_CREATE( path = sa_file, cmode=NF90_CLOBBER, ncid = sa_ncid ) )
  CALL check(NF90_DEF_DIM( sa_ncid, 'param_dim', nparams, par_dimid) )
  CALL check(NF90_DEF_DIM( sa_ncid, 'mc_dim', NF90_UNLIMITED, sa_dimid) )
  !..
  CALL check(NF90_DEF_VAR( sa_ncid, 'parameters', NF90_DOUBLE, (/par_dimid, sa_dimid/), &
                           par_varid))
  CALL check(NF90_DEF_VAR( sa_ncid, 'loglikelihood', NF90_DOUBLE, sa_dimid, lls_varid))
  !..
  CALL check(NF90_ENDDEF( sa_ncid ) )
  CALL check(NF90_CLOSE( sa_ncid) )


  !..NB: Simulated Annealing here will ONLY use uncorrelated proposal
  !..Take care here with weird parameters!
  DO i=1,nparams
    IF (pnaught(i) > 900.) THEN
      psigma_anneal(i) = init_par_var_sa*( pmax(i) - pmin(i) )
    ELSEIF (pnaught(i) > -900.) THEN
      PRINT*,'THIS IS A FUNNY VAR, var no. = ',i
      psigma_anneal(i) = init_par_var_sa
    ENDIF
  ENDDO
  !..Initialize sigmulated annealing statistics
  p_sa_mean    = params  !..The full mean
  p_sa_accmean = params  !..The mean of only accepted points (no duplicates)
  
  !..Se anealing parameters and initialize variables
  k       = 1       !
  naccept = 0       !
  
  count_params = (/nparams,1/)
  start = (/1,0/)
  CALL check(NF90_OPEN(path=sa_file, mode=NF90_WRITE, ncid = sa_ncid) )
  !DO WHILE (k < n_sa )
  sim_ann_run = .TRUE.
  DO WHILE (sim_ann_run)
  
    start(2) = start(2) + 1
    !
    IF (sim_ann_type == 1) THEN 
      invT     = 1.
    ELSE IF (sim_ann_type == 2) THEN
      invT    = REAL(k+1)*invT_scale
    ELSE IF (sim_ann_type == 3) THEN
      invT    = LOG(REAL(k+1))*invT_scale
    ELSE IF (sim_ann_type == 4) THEN
      invT    = SQRT(REAL(k+1))*invT_scale
    ELSE
      STOP 'THIS mcmc_type OPTION NOT YET CODED'
    ENDIF
    
    IF (k>1) THEN
      !..Generate perturbation parameter
      DO i = 1,nparams
        IF (pflag(i) == 1) THEN
          CALL box_muller_polar(nrand1, nrand2)
          dp(i) = nrand1*psigma_anneal(i)
        ELSE
          dp(i) = 0.
        ENDIF
      ENDDO

      p_params = params + dp
    ENDIF
    !..Note that this rangecheck will enforce the prior, but that prior will NOT be 
    !..scaled by the simulated annealing temperature
    IF ( .NOT. (rangecheck(pmin,pmax,pnaught,p_params,pflag,monotonic, &
                           gauss_prior,priormean,icovmat,prior_phi))) THEN
      PRINT*,'OUT OF RANGE, SIMULATED ANNEALING, OUT OF RANGE ! ! ! ! ! ! '
      !print*,'p_params = ',p_params
      !p_like = 0.
      p_like = 1.e-40
      p_phi  = -9.e10
    ELSE
      CALL preprocess_parameters(p_params, pmin, pmax, pnaught, p_linspace_ratio, p_mparams)
      !print*,'params SA  = ',p_params
      !print*,'mparams SA = ',p_mparams
      CALL run_fwd_model(p_mparams, p_fwd_obs, p_fwd_obs_alt)
      if (sigma_unknown) then
         call likelihood_ind_diag_sigma(obs, p_fwd_obs, n_obs_types, obs_types, p_phi, p_like, sigma_diag)
      else
         IF (use_alt_obs) THEN
            CALL likelihood_ind(obs, p_fwd_obs_alt, sigma_obs, p_phi, p_like)
         ELSE  
            CALL likelihood_ind(obs, p_fwd_obs, sigma_obs, p_phi, p_like)
         ENDIF
      end if
    ENDIF
    !..Check to see if accepted proposal is a new MLE
    IF (p_phi > phi_mle) THEN
      like_mle    = p_like
      params_mle  = p_params
      mparams_mle = p_mparams
      phi_mle     = p_phi
      fwd_obs_mle = p_fwd_obs
      if (sigma_unknown) sigma_mle = sigma_diag
    ENDIF
    
    !..Apply accept/reject criterion
    urand1 = random_real()
    urand = REAL(urand1)
    IF ( urand <= EXP((p_phi - phi)*invT) ) THEN
      IF (MOD(k,100) == 0) THEN
        PRINT*,'DIAGNOSTIC: Sim. Annealing Proposal Accepted!!'
        PRINT*,'Iter no. ',k,' of ',n_sa
        PRINT*,'Accpted total = ',naccept
        PRINT*,'urand = ',urand,'  ratio = ',EXP(p_phi-phi),'  sa_ratio = ',EXP((p_phi - phi)*invT) 
        PRINT*,'phi = ',phi,' p_phi = ',p_phi
        if (sigma_unknown) then
           print *, 'sigma_diag = ', sigma_diag
        end if
      ENDIF
      naccept     = naccept + 1
      params      = p_params  
      mparams     = p_mparams  
      like        = p_like
      phi         = p_phi
      fwd_obs     = p_fwd_obs
      fwd_obs_alt = p_fwd_obs_alt
      fwd_obs_save= fwd_obs
      !..ADD IO step here?
      IF (gauss_prior) prior_phi = return_prior_phi(params, priormean, icovmat)

    ELSE
      IF (MOD(k,100) == 0) THEN
        PRINT*,'DIAGNOSTIC: Sim Annealing proposal Rejected'
        PRINT*,'Iter no. ',k,' of ',n_sa
        PRINT*,'Accpted total = ',naccept
        PRINT*,'urand = ',urand,'  ratio = ',EXP(p_phi-phi),'  sa_ratio = ',EXP((p_phi - phi)*invT) 
        PRINT*,'phi = ',phi,' p_phi = ',p_phi
        if (sigma_unknown) then
           print *, 'sigma_diag = ', sigma_diag
        end if
      ENDIF
    ENDIF

    p_sa_mean = p_sa_mean*(k-1)/k +params/k

    !..IO step
    CALL check(NF90_PUT_VAR(sa_ncid, par_varid, mparams, start = start, count= count_params))
    CALL check(NF90_PUT_VAR(sa_ncid, lls_varid, [phi], start = [start(2)], count = [1]))

    k = k + 1
    
    IF (k > n_sa+1) THEN
      sim_ann_run = .FALSE.
      !..If the log-likelihood is too low, then stop because sampling will likely fail
      !..Minimum phi = 2xobs standard deviation = -2*nobs
      IF ( phi_mle < MIN(phi_save,-2.*(REAL(nobs)) )) THEN
        PRINT*, '---------------------------------------------------------------'
        PRINT*, '---------------------------------------------------------------'
        PRINT*, ' SIMULATED ANNEALING SAMPLER PHI_MLE TOO LOW! HAS NOT FOUND'
        PRINT*, ' A HIGH-PROBABILITY REGION. TRY A DIFFERENT RANDOM SEED!'
        PRINT*, ' phi_mle  = ',phi_mle,',   phi_min = ',-2.*(REAL(nobs))
        PRINT*, ' p_params = ',p_params
        PRINT*, ' p_mparams= ',p_mparams
        PRINT*, 'p_fwd_obs = ',p_fwd_obs
        if (sigma_unknown) print *, 'sigma_mle = ',sigma_mle
        PRINT*, '---------------------------------------------------------------'
        PRINT*, '---------------------------------------------------------------'
        STOP 'MCMC chain stopped owing to poor simulated annealing performance'
      ENDIF
    ENDIF

  ENDDO !..while k < n_sa
  CALL check(NF90_CLOSE(sa_ncid) )
  !..Start annealing loop, to conclude when k == n_sa
  
  !..Print diagnostic outputs
  PRINT*,'--------------------------------------------'
  PRINT*,'SIMULATED ANNEALING PRE-SAMPLER COMPLETE!!!'
  PRINT*,'Number of accepted samples = ', naccept
  PRINT*,'Accept percentage          = ', 100.*naccept/n_sa
  PRINT*,'MLE log-likelihood         = ',phi_mle
  PRINT*,'MLE likelihood             = ',like_mle
  PRINT*,'MLE params                 = ',params_mle
  PRINT*,'mean params                = ',p_sa_mean
  if (sigma_unknown) print *, 'sigma_mle = ',sigma_mle
  PRINT*,'--------------------------------------------'
  PRINT*,'Last proposed parameter values, etc:'
  PRINT*, ' p_params = ',p_params
  PRINT*, ' p_mparams= ',p_mparams
  PRINT*, 'p_fwd_obs = ',p_fwd_obs
  PRINT*,'--------------------------------------------'
  !..Set params to params_mle
  params  = params_mle
  mparams = mparams_mle
  phi     = phi_mle
  like    = like_mle
  fwd_obs = fwd_obs_mle
  if (gauss_prior) prior_phi = return_prior_phi(params, priormean, icovmat)

ENDIF !..sim_ann_type /= 0

!.................................................................................!
!                       ~ ~ ~ MAIN MCMC LOOP ~ ~ ~                                !
!.................................................................................!
!..Initialize counters
icount      = 1
bcount      = 1
cov_count   = 1
naccept     = 0
main_chain  = .TRUE.
delay_lev   = 1
delay_count = 0
b_test      = 0
!..Set burn-in flag to true -- NB: to skip burn-in, set to false
burn_in     = .TRUE.
adapt_count = cov_mem
!n_adapt     = 100       !..Proposal variance adapts every n_adapt samples
acc_rate    = 0.

!..Initialize random number generator (if this hasn't been done in set_guess)
!IF (set_guess) CALL RANDOM_NUMBER(urand)
IF (set_guess) urand1 = random_real()
urand = REAL(urand1)

print_diag = .TRUE.
PRINT*, 'STARTING MAIN MCMC LOOP!'

CALL check(NF90_OPEN(path=out_file, mode = NF90_WRITE, ncid=mc_ncid) )
    
!...Allow for namelist specified skipping of burn-in
IF ( skip_burn ) THEN
  burn_in = .false.
  naccept = 0
  IF ( adaptive ) THEN
    psigma_running = psigma
    pcov_running   = pcov
    adapt_count    = cov_mem
  ENDIF
ENDIF

DO WHILE ( main_chain )
 
  IF (MOD(icount,diagnostic_mod)==0) print_diag=.TRUE.

  !..Draw proposal perturbation
  found_pert  = .FALSE.
  accepted    = .FALSE.
  
  !..Generate parameter perturbation vector & parameter proposal vector
  IF (uniform_dist) THEN
    CALL perturb_uni(pflag,pmin,pmax,pnaught,p_params)
  ELSE
    IF ( cov_prop ) THEN
      CALL perturb(nparams, nprtb, 1, pflag, R, dp)
    ELSE
      dp = 0.
      DO i = 1,nparams
        IF (pflag(i) .EQ. 1) THEN
          CALL Box_Muller_polar(nrand1,nrand2)
          dp(i)  = nrand1 *  psigma(i)
        ENDIF
      ENDDO
    ENDIF 
    !..Scale parameter perturbation by delayed rejection level
    p_params    = params + dp/SQRT(REAL(delay_lev))
  ENDIF
  
  logdum = rangecheck(pmin,pmax,pnaught,p_params,pflag,monotonic, &
                      gauss_prior,priormean,icovmat,prior_phi)
  !..Determine how to handle this proposal..........................!
  !..Out of range, reject
  IF ( .NOT. (logdum) .AND. &
      (delay_lev .GE. max_delay_lev)) THEN
    IF (print_diag) THEN
      PRINT*,'DIAGNOSTIC: option 1 - out of range, reject'
      PRINT*,'Params = ',params
      PRINT*,'P_params = ',p_params
    ENDIF
    accepted  = .false.
    delaying  = .false.
  !..Out of range, delay rejection
  ELSEIF ( .NOT. (logdum) .AND. &
      (delay_lev .LT. max_delay_lev)) THEN
    IF (print_diag) PRINT*,'DIAGNOSTIC: option 2 - out of range, delay'
    accepted  = .false.
    delaying  = .true.
    p_like_save(delay_lev)  = 0.
    p_phi_save(delay_lev)   = -9999.
    dp_save(:,delay_lev)    = dp
  !..In range, run fwd model, calculate likelihood, accept probabilistically
  ELSE
    IF (print_diag) THEN 
      PRINT*,'DIAGNOSTIC: option 3 - in range'
      PRINT*,'Params = ',params
      PRINT*,'P_params = ',p_params
    ENDIF
    CALL preprocess_parameters(p_params, pmin, pmax, pnaught, p_linspace_ratio, p_mparams)
    CALL run_fwd_model(p_mparams, p_fwd_obs, p_fwd_obs_alt)
    !CALL likelihood_cov(nobs, obs, p_fwd_obs, cov_obs, p_phi, p_like)
    if (sigma_unknown) then
       call likelihood_ind_diag_sigma(obs, p_fwd_obs, n_obs_types, obs_types, p_phi, p_like, sigma_diag)
    else
       IF (use_alt_obs) THEN
          CALL likelihood_ind(obs, p_fwd_obs_alt, sigma_obs, p_phi, p_like)
       ELSE  
          CALL likelihood_ind(obs, p_fwd_obs, sigma_obs, p_phi, p_like)
       ENDIF
    end if
    IF (print_diag) PRINT*,'DIAGNOSTIC: like, p_like ',like, p_like
    IF (print_diag) PRINT*,'------------phi, p_phi ', phi, p_phi
    !...Test for accept/reject, make special consideration for delayed rejection
    IF ( delay_lev .EQ. 1 ) THEN
      IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.1 - delay lev=1'
      !CALL RANDOM_NUMBER(urand)
      urand1 = random_real()
      urand = REAL(urand1)
      !...Accept
      !IF ((urand .LE. (p_like/like)) .OR. uniform_dist) THEN
      IF ((urand .LE. REAL(EXP(p_phi-phi))) .OR. uniform_dist) THEN
        IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.1.1 - accept'
        IF (print_diag) PRINT*,'----------- urand = ',urand
        !IF (print_diag) PRINT*,'----------- p_like/like = ', (p_like/like)
        IF (print_diag) PRINT*,'----------- EXP(p_phi-phi) = ',(EXP(p_phi-phi))
        if (print_diag .and. sigma_unknown) then
           print *, 'sigma_diag = ', sigma_diag
        end if
        accepted  = .TRUE.
        delaying  = .FALSE.
      !...Reject probabilistically, reject
      ELSEIF (delay_lev .GE. max_delay_lev) THEN
        IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.1.2 - reject'
        IF (print_diag) PRINT*,'----------- urand = ',urand
        !IF (print_diag) PRINT*,'----------- p_like/like = ', (p_like/like)
        IF (print_diag) PRINT*,'----------- EXP(p_phi-phi) = ',(EXP(p_phi-phi))
        if (print_diag .and. sigma_unknown) then
           print *, 'sigma_diag = ', sigma_diag
        end if
        accepted  = .false.
        delaying  = .false.
      !...Reject probabilistically, delay rejection
      ELSE
        IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.1.3 - reject, delay'
        IF (print_diag) PRINT*,'----------- urand = ',urand
        !IF (print_diag) PRINT*,'----------- p_like/like = ', (p_like/like)
        IF (print_diag) PRINT*,'----------- EXP(p_phi-phi) = ',(EXP(p_phi-phi))
        if (print_diag .and. sigma_unknown) then
           print *, 'sigma_diag = ', sigma_diag
        end if
        accepted  = .false.
        delaying  = .true.
        p_like_save(delay_lev)  = p_like
        p_phi_save(delay_lev)   = p_phi
        dp_save(:,delay_lev)    = dp
      ENDIF
    !...2nd proposal, consider 'star' point in accept/reject criterion
    ELSEIF ( delay_lev .EQ. 2 ) THEN
      IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.2 - delay lev=2'
      delaying  = .false.
      !...Choose between y* = y  or y* = z-(y-x); (dr_star=false, true, resp.)
      IF (.NOT. (dr_star)) THEN
        !pstar  =
        phi_star = p_phi_save(delay_lev-1)
        like_star= p_like_save(delay_lev-1)
        !CALL RANDOM_NUMBER(urand)
        urand1 = random_real()
        urand = REAL(urand1)
        !IF ( urand <= MIN(1.,EXP(p_phi - phi)*(1.-MIN(1.,EXP(phi_star-p_phi))) /  &
        !                                    (1.-MIN(1.,EXP(p_phi_save(delay_lev-1)&
        !                                    -phi))))) THEN
        prop_fac = 4.
        !prop_fac = GAUSS_PROB(sig,dx)/GAUSS_PROB(sig,dx)
        !prop_fac = MULTIGAUSS_PROB(cov,dx)/MULTIGAUSS_PROB(cov,dx)
        IF ( urand .LE. REAL(MIN(1., prop_fac * EXP(p_phi-phi) *         &
                            (1.-MIN(1.,EXP(phi_star-p_phi))) /       &
                            (1.-EXP(phi_star-phi))))) THEN
          IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.2.1 - delayed, accept'
          accepted = .true.
        ELSE
          IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.2.2 - delayed, reject'
          accepted = .false.
        ENDIF
      ELSE
        pstar = p_params - dp_save(:,delay_lev-1)
        logdum = rangecheck(pmin,pmax,pnaught,pstar,pflag, monotonic, &
                            gauss_prior,priormean,icovmat,prior_phi)
        IF (logdum) THEN
          !..If pstar is in range, cacalate like_star
          CALL preprocess_parameters(pstar, pmin, pmax, pnaught, p_linspace_ratio, mpstar)
          CALL run_fwd_model(mpstar, fwd_obs_star, fwd_obs_alt_star)
          if (sigma_unknown) then
             call likelihood_ind_diag_sigma(obs, fwd_obs_star, n_obs_types, obs_types, phi_star, like_star, sigma_diag)
          else
             IF (use_alt_obs) THEN
                CALL likelihood_ind(obs, fwd_obs_alt_star, sigma_obs, phi_star, like_star)
             ELSE
                CALL likelihood_ind(obs, fwd_obs_star, sigma_obs, phi_star, like_star)
             ENDIF
          end if
          alpha_star = MIN(1., EXP(phi_star-p_phi))
          !CALL RANDOM_NUMBER(urand)
          urand1 = random_real()
          urand = REAL(urand1)
          IF ( urand .LE. REAL((EXP(p_phi-phi))*(1.-alpha_star)/                     &
                                (1.-MIN(1.,EXP(p_phi_save(delay_lev-1)-phi))))) THEN
            PRINT*,'DIAGNOSTIC: option 3.2.1 - delayed, accept'
            accepted  = .true.
          ELSE
            PRINT*,'DIAGNOSTIC: option 3.2.2 - delayed, reject'
            accepted  = .false.
          ENDIF
        !..If pstar is out of range, assume like_star = 0.
        ELSE
          !CALL RANDOM_NUMBER(urand)
          urand1 = random_real()
          urand = REAL(urand1)
          !IF (urand .LE. (p_like/like)/(1.-MIN(1.,p_like_save(delay_lev-1)/like))) &
          IF (urand.LE.REAL((EXP(p_phi-phi))/(1.-MIN(1.,p_like_save(delay_lev-1)/like)))) &
          THEN
            IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.2.3 - delayed, accept'
            accepted  = .true.
          ELSE
            IF (print_diag) PRINT*,'DIAGNOSTIC: option 3.2.4 - delayed, reject'
            accepted  = .false.
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  !...If delaying, skip main increments, skip writing to file, burn-in and adapt
  IF ( delaying ) THEN
    delay_lev   = delay_lev + 1
    delay_count = delay_count + 1
  ELSE
    IF (print_diag) PRINT*,'DIAGNOSTIC: completed iteration number..................... ', icount
    icount  = icount + 1
    !...If accepted, then interate counters, set prior = proposal, check MLE
    IF ( accepted ) THEN
      !...Reset counters, etc
      naccept       = naccept + 1
      mparams       = p_mparams
      params        = p_params
      like          = p_like
      phi           = p_phi
      fwd_obs       = p_fwd_obs
      fwd_obs_alt   = p_fwd_obs_alt
      fwd_obs_save  = fwd_obs
      !data_save     = p_data_save
      !ref_full      = p_ref_full
      IF (gauss_prior) prior_phi = return_prior_phi(params, priormean, icovmat)
      !...Check to see if accepted proposal is a new MLE
      IF (p_phi .GT. phi_mle ) THEN
        params_mle  = p_params
        mparams_mle = p_mparams
        like_mle    = p_like
        phi_mle     = p_phi
        fwd_obs_mle = p_fwd_obs
        if (sigma_unknown) sigma_mle = sigma_diag
        CALL check(NF90_PUT_VAR(mc_ncid,mle_varid,params_mle))
        CALL check(NF90_PUT_VAR(mc_ncid,mmle_varid,mparams_mle))
        CALL check(NF90_SYNC(mc_ncid))
      ENDIF
    ENDIF
    
    !..Iterate position of counter in file
    IF(MOD(icount,out_mod)==0) THEN
      start(2) = start(2) + 1
      !..Write data to file
      ios = NF90_PUT_VAR(mc_ncid,params_varid, params, start=start,count=count_params)
      ios = NF90_PUT_VAR(mc_ncid,m_params_varid, mparams, start=start,count=count_params)
      ios = NF90_PUT_VAR(mc_ncid,ll_varid, [phi], start=[start(2)],count=[1])
      IF (save_fwd_obs) ios = NF90_PUT_VAR(mc_ncid,obs_varid, fwd_obs, start=start,count=count_obs)
      IF (save_obs_alt) ios = NF90_PUT_VAR(mc_ncid,obs_alt_varid,fwd_obs_alt,start=start,count=count_obs_alt)
    ENDIF
    delay_lev = 1

    !.................................................................!
    !                      Burn - In                                  ! 
    !.................................................................!
    b_test(1:cov_check-1) = b_test(2:cov_check)
    IF ( accepted ) THEN
      b_test(cov_check) = 1
    ELSE
      b_test(cov_check) = 0
    ENDIF

    !...Compute new acceptance rate. Note: acc_rate for the first 99 model runs
    !...will be inaccurate (before the array is filled).
    acc_rate = REAL(SUM(b_test)) / REAL(cov_check)
    IF (print_diag) PRINT*,'DIAGNOSTIC: Acceptance rate =  ',acc_rate
    IF (print_diag) PRINT*,'DIAGNOSTIC: Acc rate (full) =  ', REAL(naccept)/REAL(icount)
    
    !!...Allow for namelist specified skipping of burn-in
    !IF ( skip_burn ) THEN
    !  burn_in = .false.
    !  naccept = 0.
    !  IF ( adaptive ) THEN
    !    psigma_running = psigma
    !    pcov_running   = pcov
    !    adapt_count    = cov_mem
    !  ENDIF
    !ENDIF
    
    !...START BURN-IN
    IF ( burn_in ) THEN
      
      !...Fill burn-in storage
      DO i = 1,nparams
        DO j = 1,cov_mem-1
          p_store(i,j) = p_store(i,j+1)
        ENDDO
        p_store(i,cov_mem) = params(i)
      ENDDO

      !...Check if we are at an update iteration
      IF ((bcount .GE. cov_mem) .AND. (cov_count .GE. cov_check)) THEN
        
        !...Check if burn-in is complete - if so, set burn_in to false
        IF (( icount .GT. iburn_spin ) .AND. &
            ( ABS(acc_rate - target_rate) .LE. acc_thresh)) THEN
          burn_in   = .false.
          naccept  = 0
          !......EXPERIMENTAL SECTION -----------
          !..if requested (scam_burn == .TRUE.) then set off-diagonal
          !..elements to zero. 
          IF (scam_burn) THEN
            DO i = 1,nparams
              DO j =1,nparams
                IF (i /= j) pcov(j,i) = 0.
              ENDDO
            ENDDO
            DO i = 1,nprtb
              DO j = 1,nprtb
                IF (i /= j) pcov_cut(j,i) = 0.
              ENDDO
            ENDDO
          ENDIF
          IF ( adaptive ) THEN
            psigma_running  = psigma
            pcov_running    = pcov
            adapt_count     = cov_mem
          ENDIF
        !...Make sure burn-in isn't taking too long
        ELSEIF ( (bcount*2 .GE. nmc) .OR. (bcount .GE. max_burn) ) THEN
          burn_in   = .false.
          naccept  = 0
          !......EXPERIMENTAL SECTION -----------
          !..if requested (scam_burn == .TRUE.) then set off-diagonal
          !..elements to zero. 
          IF (scam_burn) THEN
            DO i = 1,nparams
              DO j =1,nparams
                IF (i /= j) pcov(j,i) = 0.
              ENDDO
            ENDDO
            DO i = 1,nprtb
              DO j = 1,nprtb
                IF (i /= j) pcov_cut(j,i) = 0.
              ENDDO
            ENDDO
          ENDIF
          IF ( adaptive ) THEN
            psigma_running  = psigma
            pcov_running    = pcov
            adapt_count     = cov_mem
          ENDIF
        !...Otherwise, update proposal variance
        ELSE  
          cov_count = 1
          DO i = 1,nparams
            pmean(i) = SUM(p_store(i,:))/REAL(cov_mem)
          ENDDO

          !...Store old psigma values
          psigma_old  = psigma
          psigma      = 0.
          pcov_old    = pcov
          pcov        = 0.
          IF ( cov_prop )  THEN
            !...Update proposal covariance and R-matrix for cov_prop
            DO j = 1,nparams
              DO k = 1, nparams
                DO i = 1,cov_mem
                  pcov(j,k) = pcov(j,k) + (p_store(j,i)-pmean(j))*  &
                                                (p_store(k,i) - pmean(k))
                ENDDO
              ENDDO
            ENDDO
            pcov  = pcov/REAL(cov_mem - 1)
            nprtb   = 0
            DO i = 1,nparams
              IF ( pflag(i) .EQ. 1 ) THEN
                nprtb = nprtb + 1
                nprtb_tmp = 0
                DO j = 1,nparams
                  IF ( pflag(j) .EQ. 1 ) THEN
                    nprtb_tmp = nprtb_tmp + 1
                    !pcov_cut(nprtb,nprtb_tmp) = s_n_sq*pcov(i,j) !+ s_n_sq*pcond(i,j)
                    pcov_cut(nprtb,nprtb_tmp) = s_n_sq*pcov(i,j) + s_n_sq*pcond(i,j)
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
            R = pcov_cut
            !CALL spotf2('L', nprtb, R, nprtb, success)
            CALL dpotf2('L', nprtb, R, nprtb, success)
            !CALL dpotrf_f95('U', nprtb, R, nprtb, success)
            !CALL dpotrf_f95(R, 'U', info=success)
            IF ( success .NE. 0 ) STOP 'Cholesky Decomposition Failed'
            !..Zero non-diagonal elements
            DO i = 1,nprtb
              DO j = 1, i-1
                R(j,i) = 0.d0  !..If dpotf2('L',...)
              ENDDO
            ENDDO
          ELSE
            !..If not cov_prop, calculate parameter standard deviation
            DO p = 1,nparams
              DO i = 1,cov_mem
                psigma(p) = psigma(p) + (p_store(p,i) - pmean(p))**2
              ENDDO
            ENDDO
            psigma = s_n*(SQRT(psigma/REAL(cov_mem -1)))
            PRINT*,'DIAGNOSTIC, burn in! old psigma = ',psigma_old
            PRINT*,'DIAGNOSTIC, burn in! new psigma = ',psigma
          ENDIF
          !...ADD SOME DIAGNOSTIC OUTPUT!!
        ENDIF
      ELSE
        cov_count = cov_count + 1
      ENDIF

      bcount  = bcount + 1
    
    !.................................................................!
    !                   Adapt (Adaptive Metropolis)                   ! 
    !.................................................................!
    ELSEIF ( adaptive ) THEN
      !PRINT*,'DIAGNOSTIC: Adapting NAO!!!'
      !...Iterate adaptive counter
      adapt_count     = adapt_count + 1
      inv_adapt_count = 1./REAL(adapt_count)
      !...Update running mean
      pmean_new       = REAL(adapt_count-1)*inv_adapt_count*pmean + &
                        inv_adapt_count*params
      !...Update the not-quite-mean meen
      IF (adapt_count == (cov_count +1)) pmeen = pmean
      pmeen_new       = pmeen + (params -pmeen)/REAL(adapt_count)
      !..
      IF ( cov_prop ) THEN
        !PRINT*,'DIAGNOSTIC, BEFORE pcov = ',  pcov_running
        !PRINT*,'DIAGNOSTIC, old mean = ',     pmean
        !PRINT*,'DIAGNOSTIC, param value = ',  params
        !PRINT*,'DIAGNOSTIC, new mean = ',     pmean_new
        !WRITE( unit_prop ) pcov_running
        !...Update running covariance
        DO i = 1,nparams
          DO j = 1,nparams
            !...From wikipedia (yes, really)
            pcov_running(j,i) = inv_adapt_count*( (adapt_count-1.)*pcov_running(j,i) + &
                                s_n_sq*(params(j)-pmean_new(j))*(params(i)-pmean(i)) )
                                !(params(j)-pmean_new(j))*(params(i)-pmean(i)) )
            !pcov_running(j,i) = pcov_running(j,i) + (s_n_sq*                    &
            !                    (params(j)-pmean_new(j))*(params(i)-pmean(i)))
            !...not sure whats different here
            !pcov_running(j,i) = pcov_running(j,i) + 1./REAL(adapt_count-1)*( &
            !   REAL(adapt_count-1)*inv_adapt_count*(params(j)-pmean(i))*         &
            !   (params(i)-pmean(j))-pcov_running(j,i))
            !...Original (from AM Haario paper)
            !pcov_running(j,i)=REAL(adapt_count-1)*inv_adapt_count*            &
            !   pcov_running(j,i) + s_n_sq*inv_adapt_count*(REAL(adapt_count)* &
            !   (pmean(j)*pmean(i)) - REAL(adapt_count+1)*pmean_new(j)*        &
            !   pmean_new(i) + params(j)*params(i))
          ENDDO
        ENDDO
        !PRINT*,'DIAGNOSTIC, AFTER pcov = ',pcov_running
        pmean = pmean_new
        pmeen = pmeen_new
        !...If at an adapt cycle, set psigma=psigma_running
        IF ( MOD(adapt_count,n_adapt) .EQ. 0 ) THEN
          CALL check(NF90_SYNC(mc_ncid))
          PRINT*,'DIAGNOSTIC, BEFORE pcov = ',  (pcov(i,i)**2, i=1,nparams)
          !CALL disp('DIAGNOSTIC, BEFORE pcov = ',pcov)
          PRINT*,'DIAGNOSTIC, old mean = ',     pmean
          PRINT*,'DIAGNOSTIC, param value = ',  params
          PRINT*,'DIAGNOSTIC, new mean = ',     pmean_new
          PRINT*,'DIAGNOSTIC, AFTER pcov = ',   (pcov_running(i,i)**2, i=1,nparams)
          !CALL disp('DIAGNOSTIC, AFTER pcov = ',pcov_running)
          !..update pcond
          DO i = 1,nparams
            pcond(i,i) = 0.010*pcov_running(i,i)
          ENDDO
          pcov = pcov_running + pcond
          nprtb = 0
          DO i = 1,nparams
            IF ( pflag(i) .EQ. 1 ) THEN
              nprtb     = nprtb + 1
              nprtb_tmp = 0
              DO j = 1,nparams
                IF ( pflag(j) .EQ. 1 ) THEN
                  nprtb_tmp = nprtb_tmp + 1
                  pcov_cut(nprtb,nprtb_tmp) = pcov(i,j)
                ENDIF
              ENDDO
            ENDIF
          ENDDO
          R = pcov_cut
          !PRINT*,'DIAGNOSTIC, pcov_cut = ',pcov_cut
          !CALL spotf2('L', nprtb, R, nprtb, success)
          CALL dpotf2('L', nprtb, R, nprtb, success)
          !CALL dpotrf_f95('U', nprtb, R, nprtb, success)
          !CALL dpotrf_f95(R, 'U', info=success)
          IF ( success .NE. 0 ) THEN
            PRINT*,'Error in Cholesky Decomposition, Code = ',success
            STOP 'Cholesky Decomposition Failed'
          ENDIF
          !..Zero non-diagonal elements
          DO i = 1,nprtb
            DO j = 1, i-1
              R(j,i) = 0.d0  !..If dpotf2('L',...)
            ENDDO
          ENDDO
        ENDIF
      !...If we are using uncorrelated proposal... this may not work
      ELSE
        !...Update running standard deviation - NOT SURE THIS IS RIGHT!
        !...Update running standard deviation - NOT SURE THIS IS RIGHT!
        IF ( MOD(adapt_count,n_adapt) .EQ. 0 ) THEN
          PRINT*, 'DIAGNOSTIC, BEFORE psigma = ', psigma_running
          PRINT*, 'DIAGNOSTIC, old mean = ', pmean
          PRINT*, 'DIAGNOSTIC, param value ',params
          PRINT*, 'DIAGNOSTIC, new mean = ', pmean_new
        ENDIF
        !WRITE( unit_prop ) psigma_running
        !...New method for calculating std dev (Welford method)
        psigma_running = SQRT( ( (psigma_running**2)*REAL(adapt_count-2) +       &
                   s_n_sq*(params-pmeen)*(params-pmeen_new) )/REAL(adapt_count -1))
        pmean = pmean_new
        pmeen = pmeen_new
        IF ( MOD(adapt_count,n_adapt) .EQ. 0 ) THEN
          psigma = psigma_running
          PRINT*, 'DIAGNOSTIC, AFTER psigma = ', psigma_running
          !..Write this psigma to file
          !ios = NF90_OPEN(path=out_file, mode = NF90_WRITE, ncid=mc_ncid)
          ios = NF90_PUT_VAR(mc_ncid,fpsigma_varid,psigma)
          !ios = NF90_CLOSE(mc_ncid)
        ENDIF
      ENDIF
    ENDIF
    !...TEST FOR FINISHED SAMPLE.................................................!
    IF ( (icount - bcount + 1) .GE. nmc ) THEN
      main_chain = .false.
      PRINT*, 'Main Chain complete! Sampling done!'
    ENDIF
  ENDIF
  print_diag = .FALSE.

ENDDO

!ios = NF90_OPEN(path=out_file, mode = NF90_WRITE, ncid=mc_ncid)
CALL check(NF90_PUT_VAR(mc_ncid,mle_varid,params_mle))
CALL check(NF90_PUT_VAR(mc_ncid,mmle_varid,mparams_mle))
CALL check(NF90_PUT_VAR(mc_ncid,pcov_varid,pcov))
CALL check(NF90_CLOSE(mc_ncid))

print*,'done with main program! aha!'
PRINT*,'---------------------------------'
PRINT*,'FINAL DIAGNOSTICS:'
PRINT*,'---------------------------------'
PRINT*,'Acc rate (full) =  ', REAL(naccept)/REAL(icount)
IF (cov_prop) THEN
  DO i=1,nparams
    psigma(i) = SQRT(pcov(i,i))
  ENDDO
ENDIF
PRINT*,'psigma = ', psigma
PRINT*,'params_mle = ',params_mle
PRINT*,'phi_mle = ',phi_mle
if (sigma_unknown) print *, 'sigma_mle = ',sigma_mle
PRINT*,'---------------------------------'


CONTAINS

subroutine initialize_wv_sat()
  use shr_kind_mod, only: shr_kind_cm
  use shr_const_mod, only: &
       tmelt => shr_const_tkfrz, &
       h2otrip => shr_const_tktrip, &
       shr_const_mwwv, &
       shr_const_mwdair
  use shr_wv_sat_mod

  character(len=*), parameter :: wv_sat_scheme = "Flatau"

  ! Width of temperature range for "mixed phase" saturation calculations.
  ! We don't actually use these mixed phase calculations, so just use 0.
  real(8), parameter :: ttrice = 0._8

  ! Ratio of H2O to dry air mean molecular weight.
  ! Labeled "epsilo" to avoid a naming conflict with the Fortran intrinsic
  ! "epsilon", which refers to machine epsilon for a given real type.
  real(8), parameter :: epsilo = shr_const_mwwv / shr_const_mwdair

  ! Error message from initialization.
  character(len=shr_kind_cm) :: errstring

  call shr_wv_sat_init(tmelt, h2otrip, ttrice, epsilo, errstring)

  if (errstring /= "") then
     print *, "Error: ", errstring
     stop "error initializing shr_wv_sat_mod"
  end if

  if (.not. shr_wv_sat_set_default(wv_sat_scheme)) then
     stop "water vapor scheme name not recognized"
  end if

end subroutine initialize_wv_sat

!...................................................................................!
subroutine Box_Muller_polar(y1,y2)
  
  real :: y1,y2
  real :: x1,x2,ww
  real,dimension(2) :: u
  real(8) :: urand1, urand2
  
  ww=2.0
  
  do while( ww >= 1.0 ) ;
    !call random_number(u(1:2))
    urand1 = random_real()
    urand2 = random_real()
    u(1) = REAL(urand1)
    u(2) = REAL(urand2)
    x1 = 2.0 * u(1) - 1.0;
    x2 = 2.0 * u(2) - 1.0;
    ww = x1 * x1 + x2 * x2;
  end do
  
  ww = sqrt( (-2.0 * log( ww ) ) / ww );
  y1 = x1 * ww;
  y2 = x2 * ww;
  
end subroutine Box_Muller_polar
!...................................................................................!
SUBROUTINE likelihood_ind(obs, for_obs, sigma_obs, phi, like)

INTEGER                                 ::  nob           !..
REAL, INTENT(IN), DIMENSION(:)          ::  obs, for_obs  !..Obs and forward obs
REAL, INTENT(IN), DIMENSION(:)          ::  sigma_obs     !..Obs std. dev (indep.)

DOUBLE PRECISION, INTENT(OUT)                       ::  like, phi     !..Likelihood (total)
REAL                                    ::  pi
INTEGER                                 ::  j

pi = 4.*ATAN(1.)
nob = SIZE(obs)

phi = 0.
like = 0.
DO j = 1,nob
  phi  = phi  - 0.5 * (DBLE(ABS( obs(j) - for_obs(j)) )**2 / DBLE(sigma_obs(j)*sigma_obs(j)))
ENDDO

like = 0.!EXP(phi)*(1./DBLE(PRODUCT(sigma_obs)))*DBLE((1./SQRT(2.*pi))**nob)

END SUBROUTINE
!...................................................................................!
SUBROUTINE likelihood_ind_diag_sigma(obs, for_obs, n_obs_types, obs_types, phi, like, &
                                     sigma_diag)

INTEGER                                 ::  nobs          !..Total number of observations
REAL, INTENT(IN), DIMENSION(:)          ::  obs, for_obs  !..Obs and forward obs
INTEGER, INTENT(IN)                     ::  n_obs_types   !..Number of different obs types
INTEGER, INTENT(IN), DIMENSION(:)       ::  obs_types     !..Obs types

DOUBLE PRECISION, INTENT(OUT)           ::  like, phi     !..Likelihood (total)
DOUBLE PRECISION, INTENT(OUT)           ::  sigma_diag(:)  !..MLE of sigma, i.e. RMS error
INTEGER                                 ::  i

double precision :: sum_square_diffs(n_obs_types)
integer :: nobs_of_type(n_obs_types)
integer :: it

nobs = SIZE(obs)

sum_square_diffs = 0.d0
nobs_of_type = 0

do i = 1, nobs
   it = obs_types(i)
   nobs_of_type(it) = nobs_of_type(it) + 1
   sum_square_diffs(it) = sum_square_diffs(it) + (DBLE(obs(i)) - DBLE(for_obs(i)))**2
end do

phi = 0.
do i = 1, n_obs_types
   phi = phi - 0.5 * (nobs_of_type(i)-1) * log(sum_square_diffs(i))
end do

like = 0.!EXP(phi)*(1./DBLE(PRODUCT(sigma_obs)))*DBLE((1./SQRT(2.*pi))**nob)

sigma_diag = sqrt(sum_square_diffs / nobs_of_type)

END SUBROUTINE

!...................................................................................!
SUBROUTINE perturb (nx,nprt,fac,xflag,rx,dx)

INTEGER, INTENT(IN)                       :: fac, nx, nprt
INTEGER, INTENT(IN), DIMENSION(nx)        :: xflag
REAL,    INTENT(IN), DIMENSION(nprt,nprt) :: rx
REAL,    INTENT(OUT), DIMENSION(nx)       :: dx

REAL,    DIMENSION(nprt)                  :: dx_cut
REAL                                      :: nrand1,nrand2
INTEGER                                   :: counter, i

dx_cut  = 0.
dx      = 0.
counter = 1
DO i=1,nprt
  CALL Box_Muller_polar(nrand1,nrand2)
  dx_cut(i) = nrand1
ENDDO
DO i = 1, nx
  IF (xflag(i) .EQ. 1 ) THEN
    dx(i)=SUM(rx(counter,:)*dx_cut(:))/SQRT(REAL(fac))
    counter = counter + 1
  ELSE
    dx(i) = 0.
  ENDIF
ENDDO
    
END SUBROUTINE

!...................................................................................!
!SUBROUTINE perturb (nx,nprt,fac,xflag,rx,dx)
!
!INTEGER, INTENT(IN)                       :: fac, nx, nprt
!INTEGER, INTENT(IN), DIMENSION(nx)        :: xflag
!REAL,    INTENT(IN), DIMENSION(nprt,nprt) :: rx
!REAL,    INTENT(OUT), DIMENSION(nx)       :: dx
!REAL,    INTENT(OUT)                      :: rho    !..proposal probability
!
!REAL,    DIMENSION(nprt)                  :: dx_cut
!REAL                                      :: nrand1,nrand2
!INTEGER                                   :: counter, i
!
!dx_cut  = 0.
!dx      = 0.
!counter = 1
!DO i=1,nprt
!  CALL Box_Muller_polar(nrand1,nrand2)
!  dx_cut(i) = nrand1
!ENDDO
!
!DO i = 1, nx
!  IF (xflag(i) .EQ. 1 ) THEN
!    dx(i)=SUM(rx(counter,:)*dx_cut(:))/SQRT(REAL(fac))
!    counter = counter + 1
!  ELSE
!    dx(i) = 0.
!  ENDIF
!ENDDO
!
!!..INSERT HERE CODE TO CACALATE PROPOSAL PROBABILITY
!rho = 0.
!
!END SUBROUTINE
!

!...................................................................................!
LOGICAL FUNCTION rangecheck(xminin,xmaxin,xnaught,xarr,xflag,monotonic, &
                            gauss_prior,xmean,icovmat,phi)

REAL, ALLOCATABLE,  DIMENSION(:)    ::  xmin, xmax
REAL, INTENT(IN), DIMENSION(:)    ::  xarr, xminin, xmaxin, xnaught,xmean
REAL, INTENT(IN), DIMENSION(:,:)  ::  icovmat
INTEGER, INTENT(IN), DIMENSION(:) ::  xflag, monotonic
INTEGER                           ::  i, nx
!LOGICAL, INTENT(IN)       ::  monotonic
LOGICAL                   ::  log_dum
LOGICAL, INTENT(IN)       ::  gauss_prior
REAL, INTENT(IN)          ::  phi
REAL                      ::  urand

!REAL, ALLOCATABLE, DIMENSION(:)         :: diff
!!..This is the negative difference array. Negative because high temperatures are
!!..at low index values, and we are going to look for monotonic increase with
!!..temperature
!ALLOCATE( diff (nx-1))
!diff = 0.
!IF (monotonic) THEN
!  diff = xarr(1:nx-1)-xarr(2:nx)
!ENDIF
nx = SIZE(xminin)
ALLOCATE(xmin(nx))
ALLOCATE(xmax(nx))


!..For "weird variables" flagged by xnaught, range is 0-1
xmin = xminin
xmax = xmaxin
DO i=1,nx
  IF (xnaught(i) < 900.) THEN
    xmin(i) = 0.
    xmax(i) = 1.
  ENDIF
ENDDO
    

rangecheck = .TRUE.
log_dum = .TRUE.
DO  i = 1,nx
  IF ( xflag(i) .EQ. 1 ) THEN
    IF ((xarr(i) .GT. xmax(i)) .OR. (xarr(i) .LT. xmin(i))) THEN
      rangecheck = .FALSE.
      RETURN
    ENDIF
  ENDIF
ENDDO

!..We are improving the monotonic check to only look at adjacent parameters 
!..whose relative value  (smaller/larger) matches that which is specified in
!..a namelist variable.
IF (ANY(monotonic .NE. 0)) THEN
  DO i=2,nx
    IF ((monotonic(i) .NE. 0) .AND. (monotonic(i-1) .NE. 0)) THEN
      IF (monotonic(i) .GT. monotonic(i-1)) THEN
        IF (xarr(i) .LT. xarr(i-1)) THEN
          rangecheck = .FALSE.
          RETURN
        ENDIF
      ELSEIF (monotonic(i) .LT. monotonic(i-1)) THEN
        IF (xarr(i) .GT. xarr(i-1)) THEN
          rangecheck = .FALSE.
          RETURN
        ENDIF
      ENDIF
    ENDIF
  ENDDO
ENDIF

!..Here the gaussian prior is enforced. NOTE! We probably do not need to worry about 
!..xflag here because taking the ratio of probabilities should cancel out unperturbed
!..parameters
IF (gauss_prior .AND. rangecheck) THEN
  p_phi = 0.
  !..Do matrix multiplication n stuff
  DO i = 1,nx
    DO j = 1,nx
      p_phi = p_phi - 0.5*ABS((xmean(j) - xarr(j))*icovmat(j,i)*(xmean(i) - xarr(i)))
    ENDDO
  ENDDO

  urand = REAL(random_real())
  rangecheck = .FALSE.
  IF (urand .LE. REAL(EXP(p_phi - phi))) rangecheck = .TRUE.
ENDIF

!!..The following checks if parameters show a monotonic increase with temperature
!!..NOTE! It ignores the last parameter. This is because we don't care about 
!!..monotonic increase there. 
!IF (monotonic) THEN
!  IF (log_dum .and. ANY(diff(1:nx-2)<0.)) log_dum= .false.
!ENDIF
!rangecheck = log_dum
DEALLOCATE(xmin)
DEALLOCATE(xmax)

RETURN
END FUNCTION rangecheck

!...................................................................................!
REAL FUNCTION return_prior_phi(xarr,xmean,icovmat)

REAL, INTENT(IN), DIMENSION(:)    :: xarr, xmean
REAL, INTENT(IN), DIMENSION(:,:)  :: icovmat
!LOGICAL, INTENT(IN), DIMENSION(:) :: xflag
!..The rest are "private to this function
!REAL    :: phi_i
INTEGER :: i,j,np

np = SIZE(xarr)

return_prior_phi = 0.
!phi_i            = 0.

!p_phi = 0.
!..Do matrix multiplication n stuff
DO i = 1,np
  DO j = 1,np
    return_prior_phi = return_prior_phi - &
                       0.5*ABS((xmean(j) - xarr(j))*icovmat(j,i)*(xmean(i) - xarr(i)))
  ENDDO
ENDDO

!DO i = 1,np
!  IF (xflag(i)) THEN
!    DO j = 1,np
!      IF (xflag(j)) THEN
!        phi_i = -0.5 * ABS((xmean(j) - xarr(j)) * icovmat(j,i) * &
!                                  (xmean(i) - xarr(i)))
!        return_prior_phi = return_prior_phi + phi_i
!      ENDIF
!    ENDDO
!  ENDIF
!ENDDO

END


!...................................................................................!
SUBROUTINE perturb_uni (xflag,xminin,xmaxin,xnaught,x)

INTEGER, INTENT(IN), DIMENSION(:)       :: xflag
REAL,    INTENT(IN), DIMENSION(:)       :: xminin, xmaxin, xnaught
REAL,    INTENT(INOUT), DIMENSION(:)    :: x
REAL, ALLOCATABLE, DIMENSION(:)         :: xmin, xmax

REAL                                    :: urand
REAL(8)                                 :: urand1
INTEGER                                 :: i, nx

nx = SIZE(xflag)
ALLOCATE(xmin(nx))
ALLOCATE(xmax(nx))

!..For "weird variables" flagged by xnaught, range is 0-1
xmin = xminin
xmax = xmaxin
DO i=1,nx
  IF (xnaught(i) < 900.) THEN
    xmin(i) = 0.
    xmax(i) = 1.
  ENDIF
ENDDO

DO i = 1,nx
  IF (xflag(i) == 1) THEN
    urand1 = random_real()
    urand = REAL(urand1)
    !CALL RANDOM_NUMBER(urand)
    x(i) = xmin(i) + (xmax(i)-xmin(i))*urand
  ELSE
  ENDIF
ENDDO

DEALLOCATE(xmin)
DEALLOCATE(xmax)

END SUBROUTINE


!...................................................................................!
SUBROUTINE preprocess_parameters(x, xmin, xmax, xnaught, x_linspace_ratio, mx)
!..If bidded to by a non 999. value of xnaught, the assumed parameter range of 0-1 is
!..mapped to a positive and negative log-range, with a funny linear range spanning 0.
!..This range is divided into three sections, and defined by four values of x
!..     x1---------------x2--------------x3-----------x4
!..     x1 = 0, mx1 = - 10**(xmin) 
!..     x2 = ?, mx2 = - 10**(xnaught)
!..     x3 = ?, mx3 = + 10**(xnaught)
!..     x4 = 1, mx4 = + 10**(xmax)
!..
!.. OK????????????
REAL, INTENT(IN), DIMENSION(:)  :: x, xmin, xmax, xnaught, x_linspace_ratio
REAL, INTENT(OUT), DIMENSION(:) :: mx
REAL                            :: xtemp, x2, x3
INTEGER :: i,nx

nx = SIZE(x)

DO i=1,nx
  IF (xnaught(i) > 9000.) THEN
    !..This is a "regular" (not weird) log-a parameter. Should be converted to linear
    !..units. Das it!
    ! SPS: if xnaught(i) > 11111., it's negative, because I'm being weird.
    if (xnaught(i) > 11111.) then
       mx(i) = -10.**x(i)
    else
       mx(i) = 10.**x(i)
    end if
  ELSEIF (xnaught(i) < -900.) THEN
    !..we set model parameters to zero (or very very small)
    !mx(i) = -100.
    mx(i) = 0.0
  ELSEIF (xnaught(i) < 900.) THEN
    IF (xmin(i) < xnaught(i)) STOP 'xnaught cannot be > xmin!!'
    IF (xmax(i) < xnaught(i)) STOP 'xnaught cannot be > xmax!!'
    IF (x_linspace_ratio(i) > 0.5) STOP 'x_linspace_ratio should not exceed 0.5, really'
    
    x2 = (1. - x_linspace_ratio(i))*(xmin(i) - xnaught(i))/(xmin(i) + xmax(i) - 2.*xnaught(i))
    x3 = x2 + x_linspace_ratio(i)

    IF (x(i) < x2) THEN
      !..In negative logzone, scale variable to backwards logrithmic scale from pnaught->pmin
      xtemp = xmin(i) - (x(i)/x2)*(xmin(i) - xnaught(i))
      !..Convert to linear value
      xtemp = -1.*10.**xtemp
    ELSEIF (x(i) < x3) THEN
      !..In linear zone, scale between -10**pnaught and 10**pnaught
      xtemp = 10.**pnaught(i)*(2.*(x(i) - x2)/x_linspace_ratio(i) - 1.)
    ELSEIF (x(i) > x3) THEN
      !..In positive logzone, scale variable to logaritihmic scale from pnaught->pmax
      xtemp = pnaught(i) + (pmax(i) - pnaught(i))*(x(i) - x3)/(1. - x3)
      !..Convert to a linear value
      xtemp = 10.**xtemp
    ELSE
      STOP 'something is wrong in preprocess parameters!'
    ENDIF
    mx(i) = xtemp
  ELSE
    !..This is probably a b-parameter, so leave it alone
    mx(i) = x(i)
  ENDIF
ENDDO

!..Convert other a-parameters to linear space??

END SUBROUTINE

!...................................................................................!
SUBROUTINE inv_preprocess_parameters(mx, xmin, xmax, xnaught, x_linspace_ratio, x)
!..This does the inverse mapping, compared with preprocess_parameters, that is, it 
!..takes parameter values on the funky scale and converts them back to the scale 
!..upon which sampling occurs
!..     x1---------------x2--------------x3-----------x4
!..     x1 = 0, mx1 = - 10**(xmin) 
!..     x2 = ?, mx2 = - 10**(xnaught)
!..     x3 = ?, mx3 = + 10**(xnaught)
!..     x4 = 1, mx4 = + 10**(xmax)
REAL, INTENT(IN), DIMENSION(:)  :: mx, xmin, xmax, xnaught, x_linspace_ratio
REAL, INTENT(OUT), DIMENSION(:) :: x
REAL                            :: xtemp, x2, x3
INTEGER :: i,nx

nx = SIZE(x)

DO i=1,nx
  IF (xnaught(i) > 9000.) THEN 
    !..This is a "regular" log-a parameter. Convert from linear to log_10
    if (xnaught(i) > 11111.) then
       x(i) = log10(-mx(i))
    else
       x(i) = log10(mx(i))
    end if
  ELSEIF (xnaught(i) < -900.) THEN
    x(i) = 0.01
  ELSEIF (xnaught(i) < 900.) THEN
    x2 = (1. - x_linspace_ratio(i))*(xmin(i) - xnaught(i))/(xmin(i) + xmax(i) - 2.*xnaught(i))
    x3 = x2 + x_linspace_ratio(i)

    IF (mx(i) < -1.*10.**(xmin(i))) THEN
      !..Something is wrong
      print*,'i = ',i
      print*,'mx(i) = ',mx(i)
      print*,'-1*10.**(xmin(i)) = ',-1.*10.**(xmin(i))
      STOP 'value is lower than the minimum'
    ELSEIF (mx(i) < -1.*10.**(xnaught(i))) THEN
      !..We're in the negative logzone -- convert from linear to log
      xtemp = LOG10(-1.*mx(i))
      !..express as ratio of negative logzone
      xtemp = x2*(xmin(i)-xtemp)/(xmin(i)-xnaught(i))
    ELSEIF (mx(i) < 10.**(xnaught(i))) THEN
      !..We're in the linear zone
      xtemp = x2 + x_linspace_ratio(i)*(mx(i)+10.**(xnaught(i)))/(2.*10.**xnaught(i))
    ELSEIF (mx(i) < 10.**(xmax(i))) THEN
      !..Convert from linear to log
      xtemp = LOG10(mx(i))
      !..We're in the positive linear zone
      xtemp = x3 + (1. - x3)*(xtemp - xnaught(i))/(xmax(i) - xnaught(i))
    ELSE
      STOP 'value is too high!'
    ENDIF
    x(i) = xtemp
  ELSE
    x(i) = mx(i)
  ENDIF
ENDDO

END SUBROUTINE

!.......................................................................!
subroutine check(status)
integer, intent ( in) :: status

if(status /= nf90_noerr) then
  print *, trim(nf90_strerror(status))
  stop "Stopped"
end if
end subroutine check

!-----------------------------------------------------------------------------------!
end program main
