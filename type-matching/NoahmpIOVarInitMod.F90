module NoahmpIOVarInitMod

  use NoahmpIOVarType
  implicit none

contains

  subroutine NoahmpIOVarInitDefault(NoahmpIO)

    implicit none
    type(NoahmpIO_type), intent(inout) :: NoahmpIO

    associate(                                &
              XSTART  =>  NoahmpIO%XSTART    ,&
              XEND    =>  NoahmpIO%XEND      ,&
              YSTART  =>  NoahmpIO%YSTART    ,&
              YEND    =>  NoahmpIO%YEND      ,&
              KDS     =>  NoahmpIO%KDS       ,&
              KDE     =>  NoahmpIO%KDE       ,&
              NSOIL   =>  NoahmpIO%NSOIL     ,&
              NSNOW   =>  NoahmpIO%NSNOW      &
             )
    if ( .not. allocated (NoahmpIO%XLAT)      ) allocate ( NoahmpIO%XLAT      (XSTART:XEND, YSTART:YEND) ) ! latitude [radians] 
    if ( .not. allocated (NoahmpIO%WSLAKEXY)  ) allocate ( NoahmpIO%WSLAKEXY  (XSTART:XEND, YSTART:YEND) ) ! lake water storage [mm]
    if ( .not. allocated (NoahmpIO%QV_CURR)   ) allocate ( NoahmpIO%QV_CURR   (XSTART:XEND,KDS:KDE,YSTART:YEND) ) ! 3D water vapor mixing ratio [kg/kg_dry]
    if ( .not. allocated (NoahmpIO%U_PHY)     ) allocate ( NoahmpIO%U_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) ) ! 3D U wind component [m/s]
    if ( .not. allocated (NoahmpIO%V_PHY)     ) allocate ( NoahmpIO%V_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) ) ! 3D V wind component [m/s]
    if ( .not. allocated (NoahmpIO%T_PHY)     ) allocate ( NoahmpIO%T_PHY     (XSTART:XEND,KDS:KDE,YSTART:YEND) ) ! 3D atmospheric temperature valid at mid-levels [K] 
    if ( .not. allocated (NoahmpIO%SHBXY)     ) allocate ( NoahmpIO%SHBXY     (XSTART:XEND,YSTART:YEND)  ) ! bare sensible heat [w/m2]  [+ to atm]
    if ( .not. allocated (NoahmpIO%EVBXY)     ) allocate ( NoahmpIO%EVBXY      (XSTART:XEND,YSTART:YEND) ) ! bare soil evaporation [w/m2]  [+ to atm]
 
    NoahmpIO%XLAT            = undefined_real
    NoahmpIO%WSLAKEY         = undefined_real
    NoahmpIO%QV_CURR         = undefined_real
    NoahmpIO%U_PHY           = undefined_real
    NoahmpIO%V_PHY           = undefined_real
    NoahmpIO%T_PHY           = undefined_real
    NoahmpIO%SHBXY           = undefined_real
    NoahmpIO%EVBXY           = undefined_real

    end associate

  end subroutine NoahmpIOVarInitDefault

end module NoahmpIOVarInitMod
