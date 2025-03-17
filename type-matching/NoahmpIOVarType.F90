module NoahmpIOVarType

  use, intrinsic :: iso_c_binding, only: C_INT, C_DOUBLE, C_PTR, C_CHAR

  implicit none
  save
  private

  type, public :: NoahmpIO_type

    integer(C_INT), pointer                                ::  ids,ide, &          ! d -> domain 
                                                               jds,jde, &          ! d -> domain
                                                               kds,kde, &          ! d -> domain
                                                               ims,ime, &          ! m -> memory
                                                               jms,jme, &          ! m -> memory
                                                               kms,kme, &          ! m -> memory
                                                               its,ite, &          ! t -> tile
                                                               jts,jte, &          ! t -> tile
                                                               kts,kte             ! t -> tile
    integer(C_INT), pointer                                ::  ITIMESTEP           ! timestep number
    integer(C_INT), pointer                                ::  NSOIL               ! number of soil layers
    real(kind=C_DOUBLE),    allocatable, dimension(:,:)    ::  XLAT                ! latitude [rad]
    real(kind=C_DOUBLE),    allocatable, dimension(:,:,:)  ::  T_PHY               ! 3D atmospheric temperature valid at mid-levels [K]
    real(kind=C_DOUBLE),    allocatable, dimension(:,:,:)  ::  QV_CURR             ! 3D water vapor mixing ratio [kg/kg_dry]
    real(kind=C_DOUBLE),    allocatable, dimension(:,:,:)  ::  U_PHY               ! 3D U wind component [m/s]
    real(kind=C_DOUBLE),    allocatable, dimension(:,:,:)  ::  V_PHY               ! 3D V wind component [m/s]
    real(kind=C_DOUBLE),    allocatable, dimension(:,:)    ::  WSLAKEXY            ! lake water storage [mm]
    real(kind=C_DOUBLE),    allocatable, dimension(:,:)    ::  SHBXY               ! bare sensible heat [w/m2]  [+ to atm]
    real(kind=C_DOUBLE),    allocatable, dimension(:,:)    ::  EVBXY               ! bare soil evaporation [w/m2]  [+ to atm]
    CHARACTER(kind=C_CHAR, LEN=:), pointer                 ::  LLANDUSE            ! (=USGS, using USGS landuse classification)
    integer(C_INT), pointer                                ::  NTIME               ! timesteps
    integer(C_INT), pointer                                ::  NSNOW               ! number of snow layers fixed to 3
  end type NoahmpIO_type

end module NoahmpIOVarType
