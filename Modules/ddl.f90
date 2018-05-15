! AUTOMATICALLY GENERATED FILE. DO NOT EDIT!

subroutine fts_fill_cache_from_ddl(handle)
  type(fts_handle), intent(inout) :: handle

  if (handle%type=='test.LS_testtype') then
    allocate(handle%columns(0:7))
    handle%columns(0)=fts_column('booldata','[arbitrary]',PLANCK_BOOL,2,1,1,0)
    handle%columns(1)=fts_column('bytedata','[arbitrary]',PLANCK_INT8,2,2,1,0)
    handle%columns(2)=fts_column('int16data','[arbitrary]',PLANCK_INT16,2,3,1,0)
    handle%columns(3)=fts_column('int32data','[arbitrary]',PLANCK_INT32,2,4,1,0)
    handle%columns(4)=fts_column('int64data','[arbitrary]',PLANCK_INT64,2,5,1,0)
    handle%columns(5)=fts_column('float32data','[arbitrary]',PLANCK_FLOAT32,2,6,1,0)
    handle%columns(6)=fts_column('float64data','[arbitrary]',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('stringdata','[arbitrary]',PLANCK_STRING,2,8,17,0)
    return
  endif
  if (handle%type=='blob.LS_blob') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('blobdata','[arbitrary]',PLANCK_INT8,2,1,1,0)
    return
  endif
  if (handle%type=='beam.LS_beammap') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('Beamdata','[arbitrary]',PLANCK_FLOAT32,2,1,1,0)
    return
  endif
  if (handle%type=='beam.LS_beammap_pol') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('Beamdata','[arbitrary]',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('BeamdataQ','[arbitrary]',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('BeamdataU','[arbitrary]',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('BeamdataV','[arbitrary]',PLANCK_FLOAT32,2,4,1,0)
    return
  endif
  if (handle%type=='beam.LS_cart_beammap') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('Beamdata','[arbitrary]',PLANCK_FLOAT32,2,1,1,0)
    return
  endif
  if (handle%type=='beam.LS_cart_beammap_pol') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('Beamdata','[arbitrary]',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('BeamdataQ','[arbitrary]',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('BeamdataU','[arbitrary]',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('BeamdataV','[arbitrary]',PLANCK_FLOAT32,2,4,1,0)
    return
  endif
  if (handle%type=='flag.LS_flag') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('flag','[arbitrary]',PLANCK_BOOL,2,1,1,0)
    return
  endif
  if (handle%type=='constant.LS_pixwin') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('Temperature','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('Polarisation','none',PLANCK_FLOAT64,2,2,1,0)
    return
  endif
  if (handle%type=='constant.LS_ringweight') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('Weight','none',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='map.LS_map') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('I_Stokes','[arbitrary]',PLANCK_FLOAT32,2,1,1,0)
    return
  endif
  if (handle%type=='map.LS_map_dp') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('I_Stokes','[arbitrary]',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='map.LS_hitmap') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('Hits','number',PLANCK_INT32,2,1,1,0)
    return
  endif
  if (handle%type=='map.LS_map_pol') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('I_Stokes','[arbitrary]',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('Q_Stokes','[arbitrary]',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('U_Stokes','[arbitrary]',PLANCK_FLOAT32,2,3,1,0)
    return
  endif
  if (handle%type=='map.LS_map_pol_dp') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('I_Stokes','[arbitrary]',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('Q_Stokes','[arbitrary]',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('U_Stokes','[arbitrary]',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='powerspectrum.LS_powspec') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('TT','[arbitrary]',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='powerspectrum.LS_powspec_pol') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('TT','[arbitrary]',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('EE','[arbitrary]',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('BB','[arbitrary]',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('TE','[arbitrary]',PLANCK_FLOAT64,2,4,1,0)
    return
  endif
  if (handle%type=='powerspectrum.LS_powspec_pol_full') then
    allocate(handle%columns(0:5))
    handle%columns(0)=fts_column('TT','[arbitrary]',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('EE','[arbitrary]',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('BB','[arbitrary]',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('TE','[arbitrary]',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('TB','[arbitrary]',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('EB','[arbitrary]',PLANCK_FLOAT64,2,6,1,0)
    return
  endif
  if (handle%type=='alm.LS_alm') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('IndexT','number',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('RealT','[arbitrary]',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('ImagT','[arbitrary]',PLANCK_FLOAT32,2,3,1,0)
    return
  endif
  if (handle%type=='alm.LS_alm_dp') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('IndexT','number',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('RealT','[arbitrary]',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('ImagT','[arbitrary]',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='alm.LS_alm_pol') then
    allocate(handle%columns(0:8))
    handle%columns(0)=fts_column('IndexT','number',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('RealT','[arbitrary]',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('ImagT','[arbitrary]',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('IndexG','number',PLANCK_INT32,3,1,1,0)
    handle%columns(4)=fts_column('RealG','[arbitrary]',PLANCK_FLOAT32,3,2,1,0)
    handle%columns(5)=fts_column('ImagG','[arbitrary]',PLANCK_FLOAT32,3,3,1,0)
    handle%columns(6)=fts_column('IndexC','number',PLANCK_INT32,4,1,1,0)
    handle%columns(7)=fts_column('RealC','[arbitrary]',PLANCK_FLOAT32,4,2,1,0)
    handle%columns(8)=fts_column('ImagC','[arbitrary]',PLANCK_FLOAT32,4,3,1,0)
    return
  endif
  if (handle%type=='alm.LS_alm_pol_dp') then
    allocate(handle%columns(0:8))
    handle%columns(0)=fts_column('IndexT','number',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('RealT','[arbitrary]',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('ImagT','[arbitrary]',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('IndexG','number',PLANCK_INT32,3,1,1,0)
    handle%columns(4)=fts_column('RealG','[arbitrary]',PLANCK_FLOAT64,3,2,1,0)
    handle%columns(5)=fts_column('ImagG','[arbitrary]',PLANCK_FLOAT64,3,3,1,0)
    handle%columns(6)=fts_column('IndexC','number',PLANCK_INT32,4,1,1,0)
    handle%columns(7)=fts_column('RealC','[arbitrary]',PLANCK_FLOAT64,4,2,1,0)
    handle%columns(8)=fts_column('ImagC','[arbitrary]',PLANCK_FLOAT64,4,3,1,0)
    return
  endif
  if (handle%type=='toi.LS_toi') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('signal','K(Antenna)',PLANCK_FLOAT32,2,1,1,0)
    return
  endif
  if (handle%type=='toi.LS_timestamps') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('timestamp','microseconds',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='gap.LS_gaplist') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('gap_start','seconds',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('gap_end','seconds',PLANCK_FLOAT64,2,2,1,0)
    return
  endif
  if (handle%type=='pointing.LS_detpoint') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('theta','rad',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('phi','rad',PLANCK_FLOAT32,2,2,1,0)
    return
  endif
  if (handle%type=='pointing.LS_detpoint_dp') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('theta','rad',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('phi','rad',PLANCK_FLOAT64,2,2,1,0)
    return
  endif
  if (handle%type=='pointing.LS_detpoint_with_orientation') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('theta','rad',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('phi','rad',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('psi','rad',PLANCK_FLOAT32,2,3,1,0)
    return
  endif
  if (handle%type=='pointing.LS_detpoint_with_orientation_dp') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('theta','rad',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('phi','rad',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('psi','rad',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='focalplane.LS_focalplanedb') then
    allocate(handle%columns(0:22))
    handle%columns(0)=fts_column('detector','none',PLANCK_STRING,2,1,8,0)
    handle%columns(1)=fts_column('phi_uv','deg',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('theta_uv','deg',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('psi_uv','deg',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('psi_pol','deg',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('epsilon','number',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('nu_cen','Hz',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('nu_min','Hz',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('nu_max','Hz',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('f_knee','Hz',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('alpha','none',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('f_min','Hz',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('f_samp','Hz',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('tau_bol','s',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('tau_int','s',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('nread','number',PLANCK_INT32,2,16,1,0)
    handle%columns(16)=fts_column('beamfwhm','deg',PLANCK_FLOAT64,2,17,1,0)
    handle%columns(17)=fts_column('ellipticity','none',PLANCK_FLOAT64,2,18,1,0)
    handle%columns(18)=fts_column('psi_ell','deg',PLANCK_FLOAT64,2,19,1,0)
    handle%columns(19)=fts_column('net_rj','K(Antenna)sqrt(s)',PLANCK_FLOAT64,2,20,1,0)
    handle%columns(20)=fts_column('sldp_x','number',PLANCK_FLOAT64,2,21,1,0)
    handle%columns(21)=fts_column('sldp_y','number',PLANCK_FLOAT64,2,22,1,0)
    handle%columns(22)=fts_column('sldp_z','number',PLANCK_FLOAT64,2,23,1,0)
    return
  endif
  if (handle%type=='ring.LS_rings') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('ringdata','K(Antenna)',PLANCK_FLOAT32,2,1,1,0)
    return
  endif
  if (handle%type=='ringset.LS_ringset') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('ringsetdata','K(Antenna)',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('ringsets_present','number',PLANCK_INT32,3,1,1,0)
    return
  endif
  if (handle%type=='ringset.LS_ringset_dp') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('ringsetdata','K(Antenna)',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('ringsets_present','number',PLANCK_INT32,3,1,1,0)
    return
  endif
  if (handle%type=='satelliteinfo.LS_satinfo') then
    allocate(handle%columns(0:23))
    handle%columns(0)=fts_column('t_startpt_int','s',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('t_startpt_frac','s',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('satpos_x_startpt','m',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('satpos_y_startpt','m',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('satpos_z_startpt','m',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('theta_x_startpt','rad',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('phi_x_startpt','rad',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('theta_y_startpt','rad',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('phi_y_startpt','rad',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('theta_z_startpt','rad',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('phi_z_startpt','rad',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('rate_rot_z_startpt','rad/s',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('t_endpt_int','s',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('t_endpt_frac','s',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('satpos_x_endpt','m',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('satpos_y_endpt','m',PLANCK_FLOAT64,2,16,1,0)
    handle%columns(16)=fts_column('satpos_z_endpt','m',PLANCK_FLOAT64,2,17,1,0)
    handle%columns(17)=fts_column('nsamples','number',PLANCK_INT32,2,18,1,0)
    handle%columns(18)=fts_column('theta_x','rad',PLANCK_FLOAT64,3,1,1,0)
    handle%columns(19)=fts_column('phi_x','rad',PLANCK_FLOAT64,3,2,1,0)
    handle%columns(20)=fts_column('theta_y','rad',PLANCK_FLOAT64,3,3,1,0)
    handle%columns(21)=fts_column('phi_y','rad',PLANCK_FLOAT64,3,4,1,0)
    handle%columns(22)=fts_column('theta_z','rad',PLANCK_FLOAT64,3,5,1,0)
    handle%columns(23)=fts_column('phi_z','rad',PLANCK_FLOAT64,3,6,1,0)
    return
  endif
  if (handle%type=='quat.LS_satpt_quat') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('quat_w','[arbitrary]',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('quat_x','[arbitrary]',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('quat_y','[arbitrary]',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('quat_z','[arbitrary]',PLANCK_FLOAT64,2,4,1,0)
    return
  endif
  if (handle%type=='catalog.LS_pointsource_catalog') then
    allocate(handle%columns(0:11))
    handle%columns(0)=fts_column('theta_ecl','rad',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('phi_ecl','rad',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('flux_30GHz','Jy',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('flux_44GHz','Jy',PLANCK_FLOAT32,2,4,1,0)
    handle%columns(4)=fts_column('flux_70GHz','Jy',PLANCK_FLOAT32,2,5,1,0)
    handle%columns(5)=fts_column('flux_100GHz','Jy',PLANCK_FLOAT32,2,6,1,0)
    handle%columns(6)=fts_column('flux_143GHz','Jy',PLANCK_FLOAT32,2,7,1,0)
    handle%columns(7)=fts_column('flux_217GHz','Jy',PLANCK_FLOAT32,2,8,1,0)
    handle%columns(8)=fts_column('flux_353GHz','Jy',PLANCK_FLOAT32,2,9,1,0)
    handle%columns(9)=fts_column('flux_545GHz','Jy',PLANCK_FLOAT32,2,10,1,0)
    handle%columns(10)=fts_column('flux_857GHz','Jy',PLANCK_FLOAT32,2,11,1,0)
    handle%columns(11)=fts_column('source_name','none',PLANCK_STRING,2,12,32,0)
    return
  endif
  if (handle%type=='catalog.LS_pointsource_catalog_new') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('theta_ecl','rad',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('phi_ecl','rad',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('flux','K_RJ*sr',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('source_name','none',PLANCK_STRING,2,4,32,0)
    return
  endif
  if (handle%type=='catalog.LS_pointsource_catalog_pol_new') then
    allocate(handle%columns(0:5))
    handle%columns(0)=fts_column('theta_ecl','rad',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('phi_ecl','rad',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('flux','K_RJ*sr',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('source_name','none',PLANCK_STRING,2,4,32,0)
    handle%columns(4)=fts_column('polangle','deg',PLANCK_FLOAT32,2,5,1,0)
    handle%columns(5)=fts_column('polpercent','percent',PLANCK_FLOAT32,2,6,1,0)
    return
  endif
  if (handle%type=='catalog.LS_pointsource_catalog_pol') then
    allocate(handle%columns(0:29))
    handle%columns(0)=fts_column('theta_ecl','rad',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('phi_ecl','rad',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('flux_30GHz','Jy',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('polangle_30GHz','deg',PLANCK_FLOAT32,2,4,1,0)
    handle%columns(4)=fts_column('polpercent_30GHz','percent',PLANCK_FLOAT32,2,5,1,0)
    handle%columns(5)=fts_column('flux_44GHz','Jy',PLANCK_FLOAT32,2,6,1,0)
    handle%columns(6)=fts_column('polangle_44GHz','deg',PLANCK_FLOAT32,2,7,1,0)
    handle%columns(7)=fts_column('polpercent_44GHz','percent',PLANCK_FLOAT32,2,8,1,0)
    handle%columns(8)=fts_column('flux_70GHz','Jy',PLANCK_FLOAT32,2,9,1,0)
    handle%columns(9)=fts_column('polangle_70GHz','deg',PLANCK_FLOAT32,2,10,1,0)
    handle%columns(10)=fts_column('polpercent_70GHz','percent',PLANCK_FLOAT32,2,11,1,0)
    handle%columns(11)=fts_column('flux_100GHz','Jy',PLANCK_FLOAT32,2,12,1,0)
    handle%columns(12)=fts_column('polangle_100GHz','deg',PLANCK_FLOAT32,2,13,1,0)
    handle%columns(13)=fts_column('polpercent_100GHz','percent',PLANCK_FLOAT32,2,14,1,0)
    handle%columns(14)=fts_column('flux_143GHz','Jy',PLANCK_FLOAT32,2,15,1,0)
    handle%columns(15)=fts_column('polangle_143GHz','deg',PLANCK_FLOAT32,2,16,1,0)
    handle%columns(16)=fts_column('polpercent_143GHz','percent',PLANCK_FLOAT32,2,17,1,0)
    handle%columns(17)=fts_column('flux_217GHz','Jy',PLANCK_FLOAT32,2,18,1,0)
    handle%columns(18)=fts_column('polangle_217GHz','deg',PLANCK_FLOAT32,2,19,1,0)
    handle%columns(19)=fts_column('polpercent_217GHz','percent',PLANCK_FLOAT32,2,20,1,0)
    handle%columns(20)=fts_column('flux_353GHz','Jy',PLANCK_FLOAT32,2,21,1,0)
    handle%columns(21)=fts_column('polangle_353GHz','deg',PLANCK_FLOAT32,2,22,1,0)
    handle%columns(22)=fts_column('polpercent_353GHz','percent',PLANCK_FLOAT32,2,23,1,0)
    handle%columns(23)=fts_column('flux_545GHz','Jy',PLANCK_FLOAT32,2,24,1,0)
    handle%columns(24)=fts_column('polangle_545GHz','deg',PLANCK_FLOAT32,2,25,1,0)
    handle%columns(25)=fts_column('polpercent_545GHz','percent',PLANCK_FLOAT32,2,26,1,0)
    handle%columns(26)=fts_column('flux_857GHz','Jy',PLANCK_FLOAT32,2,27,1,0)
    handle%columns(27)=fts_column('polangle_857GHz','deg',PLANCK_FLOAT32,2,28,1,0)
    handle%columns(28)=fts_column('polpercent_857GHz','percent',PLANCK_FLOAT32,2,29,1,0)
    handle%columns(29)=fts_column('source_name','none',PLANCK_STRING,2,30,32,0)
    return
  endif
  if (handle%type=='catalog.LS_planet_data') then
    allocate(handle%columns(0:16))
    handle%columns(0)=fts_column('t_startpt_int','s',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('t_startpt_frac','s',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('theta_startpt','rad',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('phi_startpt','rad',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('d_sun_startpt','m',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('d_sat_startpt','m',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('angradius_startpt','arcmin',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('theta_sunsat_startpt','rad',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('t_endpt_int','s',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('t_endpt_frac','s',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('theta_endpt','rad',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('phi_endpt','rad',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('d_sun_endpt','m',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('d_sat_endpt','m',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('angradius_endpt','arcmin',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('theta_sunsat_endpt','rad',PLANCK_FLOAT64,2,16,1,0)
    handle%columns(16)=fts_column('planet_name','none',PLANCK_STRING,3,1,32,0)
    return
  endif
  if (handle%type=='module.psc.LS_pointsource_hits') then
    allocate(handle%columns(0:9))
    handle%columns(0)=fts_column('source_name','none',PLANCK_STRING,2,1,32,0)
    handle%columns(1)=fts_column('max_period','number',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('max_antenna_temp','K(Antenna)',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('max_intensity','Jy',PLANCK_FLOAT32,2,4,1,0)
    handle%columns(4)=fts_column('max_theta','deg',PLANCK_FLOAT32,2,5,1,0)
    handle%columns(5)=fts_column('max_phi','deg',PLANCK_FLOAT32,2,6,1,0)
    handle%columns(6)=fts_column('number_of_hits','number',PLANCK_INT32,2,7,1,0)
    handle%columns(7)=fts_column('periods_of_hits','number',PLANCK_INT32,2,8,1,0)
    handle%columns(8)=fts_column('contiguous_periods','number',PLANCK_INT32,2,9,1,0)
    handle%columns(9)=fts_column('Times_of_hits','none',PLANCK_STRING,2,10,128,0)
    return
  endif
  if (handle%type=='table.LS_detector_response') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('frequency','Hz',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('response','fraction',PLANCK_FLOAT64,3,1,1,0)
    return
  endif
  if (handle%type=='table.LS_toi_index_table') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('first_sample','number',PLANCK_INT64,2,1,1,0)
    handle%columns(1)=fts_column('last_sample','number',PLANCK_INT64,2,2,1,0)
    return
  endif
  if (handle%type=='ephemeris.LS_ephemeris') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('bodies','none',PLANCK_STRING,2,1,32,0)
    handle%columns(1)=fts_column('scalars','none',PLANCK_STRING,3,1,32,0)
    handle%columns(2)=fts_column('arrays','none',PLANCK_STRING,4,1,32,0)
    handle%columns(3)=fts_column('data','[arbitrary]',PLANCK_FLOAT64,5,1,1,0)
    return
  endif
  if (handle%type=='sat.LS_satpoint_real') then
    allocate(handle%columns(0:10))
    handle%columns(0)=fts_column('quat_time','seconds',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('quat_w','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('quat_x','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('quat_y','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('quat_z','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('quat_flag','none',PLANCK_INT32,2,6,1,0)
    handle%columns(6)=fts_column('pp_name','none',PLANCK_STRING,3,1,8,0)
    handle%columns(7)=fts_column('pp_tstart','seconds',PLANCK_FLOAT64,3,2,1,0)
    handle%columns(8)=fts_column('pp_tend','seconds',PLANCK_FLOAT64,3,3,1,0)
    handle%columns(9)=fts_column('pp_ifirst','number',PLANCK_INT64,3,4,1,0)
    handle%columns(10)=fts_column('pp_nquat','number',PLANCK_INT32,3,5,1,0)
    return
  endif
  if (handle%type=='sat.LS_tiltAngles') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('tilt1','arcmin',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('tilt2','arcmin',PLANCK_FLOAT64,2,2,1,0)
    return
  endif
  if (handle%type=='mission.LS_samples_real') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('pp_name','none',PLANCK_STRING,2,1,8,0)
    handle%columns(1)=fts_column('pp_nsamp','number',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('sample_index','number',PLANCK_INT64,3,1,1,0)
    handle%columns(3)=fts_column('sample_time','seconds',PLANCK_FLOAT64,3,2,1,0)
    return
  endif
  if (handle%type=='proc.List') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('objName','none',PLANCK_STRING,2,1,1999,0)
    return
  endif
  if (handle%type=='helsinki.map3d') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('pixel','number',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('ipsi','number',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('hits','number',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('tod','[arbitrary]',PLANCK_FLOAT32,2,4,1,0)
    return
  endif
  if (handle%type=='ctr.LS_ctr') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('CTR','nanosec',PLANCK_INT64,2,1,1,0)
    return
  endif
  if (handle%type=='toi.LS_temperature') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('temperatures','K',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='toi.science.LFI_Data') then
    allocate(handle%columns(0:7))
    handle%columns(0)=fts_column('obt_RC','none',PLANCK_INT16,2,1,1,0)
    handle%columns(1)=fts_column('sampleOBT','2^-16s',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('sampleSCET','microsec_since_Jan_1_1958',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('sky_adu','ADU',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('sky_volt','Volt',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('load_adu','ADU',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('load_volt','Volt',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('qualityFlag','none',PLANCK_INT32,2,8,1,0)
    return
  endif
  if (handle%type=='toi.LS_TimeTOD') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('sampleSCET','microseconds',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='toi.science.AVR') then
    allocate(handle%columns(0:6))
    handle%columns(0)=fts_column('sampleOBT','none',PLANCK_INT64,2,1,1,0)
    handle%columns(1)=fts_column('sampleSCET','none',PLANCK_INT64,2,2,1,0)
    handle%columns(2)=fts_column('sky_adu','ADU',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('sky_volt','Volt',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('load_adu','ADU',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('load_volt','Volt',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('qualityFlag','none',PLANCK_INT32,2,7,1,0)
    return
  endif
  if (handle%type=='toi.science.LFI_DataDiff') then
    allocate(handle%columns(0:4))
    handle%columns(0)=fts_column('obt_RC','none',PLANCK_INT16,2,1,1,0)
    handle%columns(1)=fts_column('sampleOBT','fracsec',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('sampleSCET','us',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('skyLoad','Volt',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('qualityFlag','none',PLANCK_INT32,2,5,1,0)
    return
  endif
  if (handle%type=='toi.attitude.HighFrequency') then
    allocate(handle%columns(0:26))
    handle%columns(0)=fts_column('dataStructureTypeID','none',PLANCK_STRING,2,1,1,0)
    handle%columns(1)=fts_column('pointingID','none',PLANCK_STRING,2,2,8,0)
    handle%columns(2)=fts_column('startTimeStablePointing_RC','none',PLANCK_INT16,2,3,1,0)
    handle%columns(3)=fts_column('startTimeStablePointing','2**-16sec',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('timeDataValid_RC','none',PLANCK_INT16,2,5,1,0)
    handle%columns(5)=fts_column('timeDataValid','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('ACMS','none',PLANCK_STRING,2,7,1,0)
    handle%columns(7)=fts_column('quaternionX','none',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('quaternionY','none',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('quaternionZ','none',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('quaternionS','none',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('eclipticLongitude','deg',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('eclipticLatitude','deg',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('spinPhaseAngle','deg',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('spinRate','deg/sec',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('solarAspect','deg',PLANCK_FLOAT32,2,16,1,0)
    handle%columns(16)=fts_column('lonDifference','deg',PLANCK_FLOAT64,2,17,1,0)
    handle%columns(17)=fts_column('latDifference','deg',PLANCK_FLOAT64,2,18,1,0)
    handle%columns(18)=fts_column('tiltAngle1','deg',PLANCK_FLOAT64,2,19,1,0)
    handle%columns(19)=fts_column('tiltAngle2','deg',PLANCK_FLOAT64,2,20,1,0)
    handle%columns(20)=fts_column('nutationAngle','deg',PLANCK_FLOAT64,2,21,1,0)
    handle%columns(21)=fts_column('bodyNutationPhaseAngle','deg',PLANCK_FLOAT64,2,22,1,0)
    handle%columns(22)=fts_column('inertialNutationPhaseAngle','deg',PLANCK_FLOAT64,2,23,1,0)
    handle%columns(23)=fts_column('slewOBT_RC','none',PLANCK_INT16,2,24,1,0)
    handle%columns(24)=fts_column('slewOBT','none',PLANCK_FLOAT64,2,25,1,0)
    handle%columns(25)=fts_column('astrFlag','none',PLANCK_STRING,2,26,1,0)
    handle%columns(26)=fts_column('astrQuality','none',PLANCK_FLOAT32,2,27,1,0)
    return
  endif
  if (handle%type=='toi.attitude.SpinPeriod') then
    allocate(handle%columns(0:15))
    handle%columns(0)=fts_column('dataStructureTypeID','none',PLANCK_STRING,2,1,1,0)
    handle%columns(1)=fts_column('pointintID','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('startTimeStablePointing','none',PLANCK_INT64,2,3,1,0)
    handle%columns(3)=fts_column('timeDataValid','none',PLANCK_INT64,2,4,1,0)
    handle%columns(4)=fts_column('eclipticLongitude','deg',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('eclipticLatitude','deg',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('spinPhaseAngle','deg',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('spinRate','deg/sec',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('solarAspect','deg',PLANCK_FLOAT32,2,9,1,0)
    handle%columns(9)=fts_column('lonDifference','deg',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('latDifference','deg',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('tiltAngle1','deg',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('tiltAngle2','deg',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('nutationAngle','deg',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('bodyNatationPhaseANgle','deg',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('intertialNutationPhaseANgle','deg',PLANCK_FLOAT64,2,16,1,0)
    return
  endif
  if (handle%type=='toi.attitude.ObservationFrequency') then
    allocate(handle%columns(0:20))
    handle%columns(0)=fts_column('dataStructureTypeID','none',PLANCK_STRING,2,1,1,0)
    handle%columns(1)=fts_column('pointingID','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('startTimeStablePointing','none',PLANCK_INT64,2,3,1,0)
    handle%columns(3)=fts_column('timeDataValid','none',PLANCK_INT64,2,4,1,0)
    handle%columns(4)=fts_column('eclipticLongitude','deg',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('eclipticLatitude','deg',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('spinPhaseAngle','deg',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('spinRate','deg/sec',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('solarAspect','deg',PLANCK_FLOAT32,2,9,1,0)
    handle%columns(9)=fts_column('lonDifference','deg',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('latDifference','deg',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('tiltAngle1','deg',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('tiltAngle2','deg',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('nutationAngle','deg',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('bodyNutationPhaseAngle','deg',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('inertialNutationPhaseAngle','deg',PLANCK_FLOAT64,2,16,1,0)
    handle%columns(16)=fts_column('obsDuration','sec',PLANCK_FLOAT32,2,17,1,0)
    handle%columns(17)=fts_column('nutSpinRate','none',PLANCK_FLOAT64,2,18,1,0)
    handle%columns(18)=fts_column('imbalance','none',PLANCK_FLOAT64,2,19,1,0)
    handle%columns(19)=fts_column('azimuthAngle','deg',PLANCK_FLOAT64,2,20,1,0)
    handle%columns(20)=fts_column('nutTimeConstant','s',PLANCK_FLOAT64,2,21,1,0)
    return
  endif
  if (handle%type=='toi.LFI_Detpoint') then
    allocate(handle%columns(0:4))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('sampleOBT','fracsec',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('sampleSCET','us',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('theta','rad',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('phi','rad',PLANCK_FLOAT64,2,5,1,0)
    return
  endif
  if (handle%type=='toi.LFI_Detpoint_pol') then
    allocate(handle%columns(0:5))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('sampleOBT','fracsec',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('sampleSCET','us',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('theta','rad',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('phi','rad',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('psi','rad',PLANCK_FLOAT64,2,6,1,0)
    return
  endif
  if (handle%type=='map.LFI_Map') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('I_Stokes','none',PLANCK_FLOAT32,2,1,1,0)
    return
  endif
  if (handle%type=='map.LFI_Map_pol') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('I_Stokes','none',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('Q_Stokes','none',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('U_Stokes','none',PLANCK_FLOAT32,2,3,1,0)
    return
  endif
  if (handle%type=='indb.LFI_instrumentdb') then
    allocate(handle%columns(0:15))
    handle%columns(0)=fts_column('detector','none',PLANCK_STRING,2,1,9,0)
    handle%columns(1)=fts_column('theta_uv','deg',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('phi_uv','deg',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('psi_uv','deg',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('psi_pol','deg',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('beamfwhm','deg',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('ellipticity','none',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('psi_ell','deg',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('nu_cen','Hz',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('nu_min','Hz',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('nu_max','Hz',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('f_knee','Hz',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('alpha','none',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('f_samp','Hz',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('tau_int','s',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('NET_RJ','K',PLANCK_FLOAT64,2,16,1,0)
    return
  endif
  if (handle%type=='polar.base') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('baseline','K',PLANCK_FLOAT64,2,2,1,0)
    return
  endif
  if (handle%type=='polar.matrix') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('cct','none',PLANCK_FLOAT32,2,1,1,0)
    return
  endif
  if (handle%type=='polar.matrix_dp') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('cct','none',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='polar.matrix_full') then
    allocate(handle%columns(0:6))
    handle%columns(0)=fts_column('cct','none',PLANCK_FLOAT32,2,1,1,0)
    handle%columns(1)=fts_column('cc11','none',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('cc12','none',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('cc13','none',PLANCK_FLOAT32,2,4,1,0)
    handle%columns(4)=fts_column('cc21','none',PLANCK_FLOAT32,2,5,1,0)
    handle%columns(5)=fts_column('cc22','none',PLANCK_FLOAT32,2,6,1,0)
    handle%columns(6)=fts_column('cc23','none',PLANCK_FLOAT32,2,7,1,0)
    return
  endif
  if (handle%type=='polar.matrix_full_dp') then
    allocate(handle%columns(0:6))
    handle%columns(0)=fts_column('cct','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('cc11','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('cc12','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('cc13','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('cc21','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('cc22','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('cc23','none',PLANCK_FLOAT64,2,7,1,0)
    return
  endif
  if (handle%type=='polar.hitmap_pol') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('Hits_unpol','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('Hits_pol','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('Det','none',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('ScalDet','none',PLANCK_FLOAT32,2,4,1,0)
    return
  endif
  if (handle%type=='polar.hitmap_mask') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('Hits','none',PLANCK_INT32,2,1,1,0)
    return
  endif
  if (handle%type=='rparam.R_table') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('r','none',PLANCK_FLOAT64,2,2,1,0)
    return
  endif
  if (handle%type=='noise_ps.LFI_noise_ps') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('freq','Hz',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('pwr','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('sigma','none',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='noise_ps.LFI_noise_filt') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('filter','1/K2',PLANCK_FLOAT64,2,1,1,0)
    return
  endif
  if (handle%type=='noise_ps.parameters') then
    allocate(handle%columns(0:5))
    handle%columns(0)=fts_column('alpha','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('s_alpha','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('fknee','Hz',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('s_fknee','Hz',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('sigma2','K2',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('s_sigma2','K2',PLANCK_FLOAT64,2,6,1,0)
    return
  endif
  if (handle%type=='flag.LFI_flag') then
    allocate(handle%columns(0:0))
    handle%columns(0)=fts_column('flag','none',PLANCK_INT32,2,1,1,0)
    return
  endif
  if (handle%type=='calib.LFI_Gain_table') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('gain','V/K',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('base','V',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='xy.LFI_beam_xy') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('x','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('y','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('distance','A.U.',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('signal','none',PLANCK_FLOAT64,2,4,1,0)
    return
  endif
  if (handle%type=='map.ifca.LFI_ps_catalog') then
    allocate(handle%columns(0:9))
    handle%columns(0)=fts_column('theta','deg',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('phi','deg',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('estim_ampl_fit','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('estim_ampl_direct','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('chi2','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('VecCoeff','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('estim_ampl_map','none',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('n_sigmas','none',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('amplification','none',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('probability','none',PLANCK_FLOAT64,2,10,1,0)
    return
  endif
  if (handle%type=='map.LFI_patch') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('x','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('y','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('value','none',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='map.LFI_SZ_catalog') then
    allocate(handle%columns(0:9))
    handle%columns(0)=fts_column('x','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('y','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('Yc','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('Delta_Yc','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('Rc','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('Delta_Rc','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('Flux_int','none',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('Delta_flux_int','none',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('theta','deg',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('phi','deg',PLANCK_FLOAT64,2,10,1,0)
    return
  endif
  if (handle%type=='map.ZLE_Micromap') then
    allocate(handle%columns(0:15))
    handle%columns(0)=fts_column('I_SPIN','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('t_startpt_int','sec',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('t_startpt_frac','sec',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('t_endpt_int','sec',PLANCK_INT32,2,4,1,0)
    handle%columns(4)=fts_column('t_endpt_frac','sec',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('SPIN_X','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('SPIN_Y','none',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('SPIN_Z','none',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('OBS_X','au',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('OBS_Y','au',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('OBS_Z','au',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('SUN_X','au',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('SUN_Y','au',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('SUN_Z','au',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('IDX0','none',PLANCK_INT32,2,15,1,0)
    handle%columns(15)=fts_column('ZLE','MJy/sterad',PLANCK_FLOAT64,3,1,1,0)
    return
  endif
  if (handle%type=='limo.Spikes') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('frequency','Hz',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('amplitude','WN',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('phase','rad',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='limo.RaaNonLinearParam') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('nonLinearity','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('linearGain','V/K',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('noiseT','K',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='limo.TransferFunction') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('amps4K','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('phi4K','rad',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('ampsFem','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('phiFem','rad',PLANCK_FLOAT64,2,4,1,0)
    return
  endif
  if (handle%type=='noise.Covariance_Matrix') then
    allocate(handle%columns(0:5))
    handle%columns(0)=fts_column('II','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('IQ','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('IU','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('QQ','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('QU','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('UU','none',PLANCK_FLOAT64,2,6,1,0)
    return
  endif
  if (handle%type=='beam.BeamWF') then
    allocate(handle%columns(0:36))
    handle%columns(0)=fts_column('multipole','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('TT_TT','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('EE_EE','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('BB_BB','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('TE_TE','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('TB_TB','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('EB_EB','none',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('TT_EE','none',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('TT_BB','none',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('TT_TE','none',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('TT_TB','none',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('TT_EB','none',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('EE_TT','none',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('EE_BB','none',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('EE_TE','none',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('EE_TB','none',PLANCK_FLOAT64,2,16,1,0)
    handle%columns(16)=fts_column('EE_EB','none',PLANCK_FLOAT64,2,17,1,0)
    handle%columns(17)=fts_column('BB_TT','none',PLANCK_FLOAT64,2,18,1,0)
    handle%columns(18)=fts_column('BB_EE','none',PLANCK_FLOAT64,2,19,1,0)
    handle%columns(19)=fts_column('BB_TE','none',PLANCK_FLOAT64,2,20,1,0)
    handle%columns(20)=fts_column('BB_TB','none',PLANCK_FLOAT64,2,21,1,0)
    handle%columns(21)=fts_column('BB_EB','none',PLANCK_FLOAT64,2,22,1,0)
    handle%columns(22)=fts_column('TE_TT','none',PLANCK_FLOAT64,2,23,1,0)
    handle%columns(23)=fts_column('TE_EE','none',PLANCK_FLOAT64,2,24,1,0)
    handle%columns(24)=fts_column('TE_BB','none',PLANCK_FLOAT64,2,25,1,0)
    handle%columns(25)=fts_column('TE_TB','none',PLANCK_FLOAT64,2,26,1,0)
    handle%columns(26)=fts_column('TE_EB','none',PLANCK_FLOAT64,2,27,1,0)
    handle%columns(27)=fts_column('TB_TT','none',PLANCK_FLOAT64,2,28,1,0)
    handle%columns(28)=fts_column('TB_EE','none',PLANCK_FLOAT64,2,29,1,0)
    handle%columns(29)=fts_column('TB_BB','none',PLANCK_FLOAT64,2,30,1,0)
    handle%columns(30)=fts_column('TB_TE','none',PLANCK_FLOAT64,2,31,1,0)
    handle%columns(31)=fts_column('TB_EB','none',PLANCK_FLOAT64,2,32,1,0)
    handle%columns(32)=fts_column('EB_TT','none',PLANCK_FLOAT64,2,33,1,0)
    handle%columns(33)=fts_column('EB_EE','none',PLANCK_FLOAT64,2,34,1,0)
    handle%columns(34)=fts_column('EB_BB','none',PLANCK_FLOAT64,2,35,1,0)
    handle%columns(35)=fts_column('EB_TE','none',PLANCK_FLOAT64,2,36,1,0)
    handle%columns(36)=fts_column('EB_TB','none',PLANCK_FLOAT64,2,37,1,0)
    return
  endif
  if (handle%type=='map.HitsMap') then
    allocate(handle%columns(0:12))
    handle%columns(0)=fts_column('Hits_total','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('Hits_01','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('Hits_02','none',PLANCK_INT32,2,3,1,0)
    handle%columns(3)=fts_column('Hits_03','none',PLANCK_INT32,2,4,1,0)
    handle%columns(4)=fts_column('Hits_04','none',PLANCK_INT32,2,5,1,0)
    handle%columns(5)=fts_column('Hits_05','none',PLANCK_INT32,2,6,1,0)
    handle%columns(6)=fts_column('Hits_06','none',PLANCK_INT32,2,7,1,0)
    handle%columns(7)=fts_column('Hits_07','none',PLANCK_INT32,2,8,1,0)
    handle%columns(8)=fts_column('Hits_08','none',PLANCK_INT32,2,9,1,0)
    handle%columns(9)=fts_column('Hits_09','none',PLANCK_INT32,2,10,1,0)
    handle%columns(10)=fts_column('Hits_10','none',PLANCK_INT32,2,11,1,0)
    handle%columns(11)=fts_column('Hits_11','none',PLANCK_INT32,2,12,1,0)
    handle%columns(12)=fts_column('Hits_12','none',PLANCK_INT32,2,13,1,0)
    return
  endif
  if (handle%type=='polar.LFI_Base') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('sampleSCET','microsec_since_01011958',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('baseline','K',PLANCK_FLOAT64,2,3,1,0)
    return
  endif
  if (handle%type=='healpix.hierarchicalMap') then
    allocate(handle%columns(0:1))
    handle%columns(0)=fts_column('geometry','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('data','[arbitrary]',PLANCK_FLOAT32,3,1,1,0)
    return
  endif
  if (handle%type=='tod_binned_2d') then
    allocate(handle%columns(0:2))
    handle%columns(0)=fts_column('pixel','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('weight','none',PLANCK_FLOAT32,2,2,1,0)
    handle%columns(2)=fts_column('signal','none',PLANCK_FLOAT32,2,3,1,0)
    return
  endif
  if (handle%type=='tod_binned_3d') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('pixel','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('ipsi','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('weight','none',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('signal','none',PLANCK_FLOAT32,2,4,1,0)
    return
  endif
  if (handle%type=='tod_binned_2d_perPID') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('pixel','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('weight','none',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('signal','none',PLANCK_FLOAT32,2,4,1,0)
    return
  endif
  if (handle%type=='tod_binned_3d_perPID') then
    allocate(handle%columns(0:9))
    handle%columns(0)=fts_column('pixel','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('ipsi','none',PLANCK_INT16,2,2,1,0)
    handle%columns(2)=fts_column('weight','none',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('signal1','none',PLANCK_FLOAT32,2,4,1,0)
    handle%columns(4)=fts_column('signal2','none',PLANCK_FLOAT32,2,5,1,0)
    handle%columns(5)=fts_column('pointingID','none',PLANCK_INT32,3,1,1,0)
    handle%columns(6)=fts_column('sample_offset','none',PLANCK_INT64,3,2,1,0)
    handle%columns(7)=fts_column('nsamples_PID','none',PLANCK_INT32,3,3,1,0)
    handle%columns(8)=fts_column('start_SCET','mus',PLANCK_FLOAT64,3,4,1,0)
    handle%columns(9)=fts_column('end_SCET','mus',PLANCK_FLOAT64,3,5,1,0)
    return
  endif
  if (handle%type=='map_4D') then
    allocate(handle%columns(0:10))
    handle%columns(0)=fts_column('pixel','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('ipsi','none',PLANCK_INT16,2,2,1,0)
    handle%columns(2)=fts_column('weight','none',PLANCK_FLOAT32,2,3,1,0)
    handle%columns(3)=fts_column('signal1','none',PLANCK_FLOAT32,2,4,1,0)
    handle%columns(4)=fts_column('signal2','none',PLANCK_FLOAT32,2,5,1,0)
    handle%columns(5)=fts_column('pointingID','none',PLANCK_INT32,3,1,1,0)
    handle%columns(6)=fts_column('sample_offset','none',PLANCK_INT64,3,2,1,0)
    handle%columns(7)=fts_column('nsamples_PID','none',PLANCK_INT32,3,3,1,0)
    handle%columns(8)=fts_column('start_SCET','mus',PLANCK_FLOAT64,3,4,1,0)
    handle%columns(9)=fts_column('end_SCET','mus',PLANCK_FLOAT64,3,5,1,0)
    handle%columns(10)=fts_column('searchtable','none',PLANCK_INT32,4,1,1,0)
    return
  endif
  if (handle%type=='table.elina_sampleinfo') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('tot_samples','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('start_time_SCET','mus',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('end_time_SCET','mus',PLANCK_FLOAT64,2,4,1,0)
    return
  endif
  if (handle%type=='calib.LFI_Gain_table_with_errors') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('pointingID','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('gain','V/K',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('base','V',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('gain_error','1/K',PLANCK_FLOAT64,2,4,1,0)
    return
  endif
  if (handle%type=='sparse_covariance') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('index2','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('covariance','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('index1','none',PLANCK_INT32,3,1,1,0)
    handle%columns(3)=fts_column('ncols','none',PLANCK_INT32,3,2,1,0)
    return
  endif
  if (handle%type=='deconv_matrix_unpol') then
    allocate(handle%columns(0:4))
    handle%columns(0)=fts_column('index2','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('TT_re','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('TT_im','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('index1','none',PLANCK_INT32,3,1,1,0)
    handle%columns(4)=fts_column('ncols','none',PLANCK_INT32,3,2,1,0)
    return
  endif
  if (handle%type=='dense_matrix') then
    allocate(handle%columns(0:3))
    handle%columns(0)=fts_column('TT_re','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('TT_im','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('index1','none',PLANCK_INT32,3,1,1,0)
    handle%columns(3)=fts_column('index2','none',PLANCK_INT32,4,1,1,0)
    return
  endif
  if (handle%type=='dense_matrix_pol') then
    allocate(handle%columns(0:19))
    handle%columns(0)=fts_column('TT_re','none',PLANCK_FLOAT64,2,1,1,0)
    handle%columns(1)=fts_column('TT_im','none',PLANCK_FLOAT64,2,2,1,0)
    handle%columns(2)=fts_column('TG_re','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('TG_im','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('TC_re','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('TC_im','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('GT_re','none',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('GT_im','none',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('GG_re','none',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('GG_im','none',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('GC_re','none',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('GC_im','none',PLANCK_FLOAT64,2,12,1,0)
    handle%columns(12)=fts_column('CT_re','none',PLANCK_FLOAT64,2,13,1,0)
    handle%columns(13)=fts_column('CT_im','none',PLANCK_FLOAT64,2,14,1,0)
    handle%columns(14)=fts_column('CG_re','none',PLANCK_FLOAT64,2,15,1,0)
    handle%columns(15)=fts_column('CG_im','none',PLANCK_FLOAT64,2,16,1,0)
    handle%columns(16)=fts_column('CC_re','none',PLANCK_FLOAT64,2,17,1,0)
    handle%columns(17)=fts_column('CC_im','none',PLANCK_FLOAT64,2,18,1,0)
    handle%columns(18)=fts_column('index1','none',PLANCK_INT32,3,1,1,0)
    handle%columns(19)=fts_column('index2','none',PLANCK_INT32,4,1,1,0)
    return
  endif
  if (handle%type=='pixel_covariance_pol') then
    allocate(handle%columns(0:12))
    handle%columns(0)=fts_column('ipix1','none',PLANCK_INT32,2,1,1,0)
    handle%columns(1)=fts_column('ipix2','none',PLANCK_INT32,2,2,1,0)
    handle%columns(2)=fts_column('II','none',PLANCK_FLOAT64,2,3,1,0)
    handle%columns(3)=fts_column('IQ','none',PLANCK_FLOAT64,2,4,1,0)
    handle%columns(4)=fts_column('IU','none',PLANCK_FLOAT64,2,5,1,0)
    handle%columns(5)=fts_column('QI','none',PLANCK_FLOAT64,2,6,1,0)
    handle%columns(6)=fts_column('QQ','none',PLANCK_FLOAT64,2,7,1,0)
    handle%columns(7)=fts_column('QU','none',PLANCK_FLOAT64,2,8,1,0)
    handle%columns(8)=fts_column('UI','none',PLANCK_FLOAT64,2,9,1,0)
    handle%columns(9)=fts_column('UQ','none',PLANCK_FLOAT64,2,10,1,0)
    handle%columns(10)=fts_column('UU','none',PLANCK_FLOAT64,2,11,1,0)
    handle%columns(11)=fts_column('index1','none',PLANCK_INT32,3,1,1,0)
    handle%columns(12)=fts_column('ncols','none',PLANCK_INT32,3,2,1,0)
    return
  endif
  call fatal_error('type '''//trim(handle%type)//''' not found in DDL')
end subroutine
