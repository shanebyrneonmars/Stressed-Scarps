def process_mcd_output_for_DISORT(ls,fname):
    tname, cname, rname, cols, rows, vals = read_mcd_output(fname)
    tname2, cname2, rname2, cols2, rows2, vals2 = read_mcd_output(fname.replace('PTDI','TAU'))

    H_LYR             = cols                     # Layer edges for the atmospheric profile
    pressure_profile  = interpn((cols,rows), vals[0,:,:], np.transpose(np.stack((H_LYR,np.ones(len(H_LYR))*lsubs))))
    TEMPER            = interpn((cols,rows), vals[1,:,:], np.transpose(np.stack((H_LYR,np.ones(len(H_LYR))*lsubs))))
    column_density    = -(H_LYR[:-1] - H_LYR[1:]) * ((pressure_profile[:-1] + pressure_profile[1:])/2.0) / ((TEMPER[:-1] + TEMPER[1:])/2.0) / 1.380649e-23

    altitude_midpoint = (H_LYR[:-1] + H_LYR[1:]) / 2.0              # mid point of layers
    dust_profile      = interpn((cols,rows), vals[2,:,:], np.transpose(np.stack((altitude_midpoint,np.ones(len(altitude_midpoint))*lsubs))))
    ice_profile       = interpn((cols,rows), vals[3,:,:], np.transpose(np.stack((altitude_midpoint,np.ones(len(altitude_midpoint))*lsubs))))

    vals2  = np.mean(vals2,axis=1)                # Average over all local times
    taud   = np.interp(lsubs,rows2,vals2[0,:])    # It should already be extinction rather than absorption and normalized by the surface pressure.
    taui   = np.interp(lsubs,rows2,vals2[1,:])    # This is actually column mass of ice, not optical depth.
    taud   = taud * 0.5                           # Convert optical depth from visible to 9.3 microns.
    taui   = taui * (3.0*2.62)/(4.0*920.0*2e-6)   # Convert to extinction optical depth. Qext for water ice at 9.3 microns is 2.62, density is 920 kg/m^3, and ice particle radius is 2 microns.

    return np.flip(H_LYR), np.flip(pressure_profile), np.flip(TEMPER), np.flip(column_density), np.flip(dust_profile), np.flip(ice_profile), taud, taui


def read_mcd_output(fname):
    # Read the entire file into a string
    with open(fname, 'r') as f:
        lines = f.readlines()

    #
    # Strip off the comments and split the table away from the row values
    #
    barelines       = [line for line in lines if not line.startswith(('#', '--'))]
    split_barelines = [line.split('||') for line in barelines]
    NL              = len(lines)     # Number of lines in the file
    NBL             = len(barelines) # Number of lines in the file without comments

    #
    # query the file for the table names, column names, and row names
    # Perform some sanity checks to check the file is structured as expected
    #
    tname   = [line for line in lines if "### Columns 2+ are " in line]
    tname   = [line.replace("### Columns 2+ are ", '', 1) for line in tname]

    colname = [line for line in lines if "### Line 1 is " in line]
    colname = [line.replace("### Line 1 is ", '', 1) for line in colname]

    rowname = [line for line in lines if "### Column 1 is " in line]
    rowname = [line.replace("### Column 1 is ", '', 1) for line in rowname]

    if not (len(tname) == len(colname) == len(rowname)):
        raise Exception("Headers in the file are not standard MCD output")
    NT      = len(tname)
    tname   = [line[:-1] for line in tname]
    print(NT, " tables found in the file:", fname)
    for i in range(NT):
        print("Table ", i+1, tname[i])

    if not (all(element == colname[0] for element in colname)):
        raise Exception("Columns in each Table have different names.")
    cname = colname[0][:-1]
    print("Columns are: ", cname)

    if not (all(element == rowname[0] for element in rowname)):
        raise Exception("Rows in each Table have different names.")
    rname = rowname[0][:-1]
    print("Rows are: ", rname)

    #
    # Extract column values, check they're the same for each table
    #
    colval = [line for line in lines if ("||" in line) and ("--" in line)]
    if not all(element == colval[0] for element in colval):
        raise Exception("Columns values in each Table are not the same.")
    colval = colval[0].split('||')[1][:-1]  # Extract the column values as a single string
    cols   = np.array(colval.split()).astype(float)       # Split columns into an array
    NC     = len(cols)                      # Number of columns in the file 


    #
    # Extract Row values, check they're the same for each table
    #
    if not (NL % NT == 0) or not (NBL % NT ==0):
        raise Exception("Numbers of lines per table seems to differ")
    NR     = int(NBL/NT)       # Number of rows per Table
    rowval = np.array([float(line[0]) for line in split_barelines]).reshape(NT, NR)
    for i in range(NT):
        if not all(element == rowval[0,i] for element in rowval[:,i]):
            raise Exception("Row values in each Table are not the same.")
    rows   = rowval[0,:]

    #
    # Extract tables and store in a 3D array
    #

    vals_string = [line[1] for line in split_barelines]
    vals = np.zeros([NT, NC, NR])
    for i in range(NT):
        for j in range(NR):
            vals[i,:,j] = np.array(vals_string[i*NR+j].split()) # for each row in each table, split the columns into an array and store


    return tname, cname, rname, cols, rows, vals


def FindSlopeFluxes(savename,slope,aspect):

    # Calculates the direct and diffuse fluxes as a function of solar incidence angle
    # Also returns the direct flux with the zero atmosphere approx. (geometry only)
    # Requires the standard DISORT results pickle file and the slope and aspect (in radians)

    with open(savename, 'rb') as f:
        data = pickle.load(f)
    if isinstance(data, dict):
        for key in data:
            globals()[key] = data[key]
    del f,key,data

    NUMU = uu_surf_inc.shape[0]
    NAZ  = uu_surf_inc.shape[1]
    NINC = uu_surf_inc.shape[2]
    NALB = uu_surf_inc.shape[3]

    UMU_up = 0.5*(np.roll(UMU, 1,axis=0) + UMU)
    UMU_dn = 0.5*(np.roll(UMU,-1,axis=0) + UMU)
    UMU_up[0]  = -1.0
    UMU_dn[-1] = 1.0
    inc_flat = np.zeros([NUMU,NAZ])
    inc_up   = np.zeros([NUMU,NAZ])
    inc_dn   = np.zeros([NUMU,NAZ])
    for i in range(NAZ):
        inc_flat[:,i] = np.arccos(-UMU)
        inc_up[:,i]   = np.arccos(-UMU_up)
        inc_dn[:,i]   = np.arccos(-UMU_dn)

    phi_flat = np.zeros([NUMU,NAZ])
    for i in range(NUMU):
        phi_flat[i,:] = np.deg2rad(PHI)

    dphi = 2.0*np.pi/NAZ
    dinc = np.abs(inc_up-inc_dn)
    sa   = dphi*dinc*0.5*(np.sin(inc_up) + np.sin(inc_dn))

    cosi = np.cos(slope)*np.cos(inc_flat) + np.sin(slope)*np.sin(inc_flat)*np.cos(phi_flat-aspect)
    cosi = np.where(cosi < 0.0, 0.0, cosi)

    diffuse = np.zeros([NINC,NALB])
    for i in range(NINC):
        for k in range(0,NALB):
            diffuse[i,k] = np.sum(uu_surf_inc[:,:,i,k]*sa*cosi)

    
    direct  = 0.0
    atmless = 0.0
    if PLANK == False:
        direct = np.zeros([NINC,NALB])
        cosi_direct = np.cos(slope)*UMU0_i + np.sin(slope)*np.sqrt(1.0-UMU0_i**2)*np.cos(np.deg2rad(PHI0)-aspect)
        cosi_direct = np.where(cosi_direct < 0.0, 0.0, cosi_direct)
        atmless = totalvisflux * cosi_direct

        for k in range(0,NALB):
            direct[:,k] = rfldir_inc[:,k] * cosi_direct / UMU0_i

    return UMU0_i, BTEMP_i, ALB_i, direct, diffuse, atmless



def MartianSolarDay(NSTP,LsubS,Lat):

    # Calculates NSTP values of the Solar Azimuth and Incidence angles (in radians) for one Martian day
    # Requires scalars for the solar longitude (in degrees) and Latitude (in degrees)

    sin_dec  = np.sin(np.deg2rad(LsubS))*np.sin(np.deg2rad(25.1))
    cos_dec  = np.sqrt(1.0-sin_dec**2)
    sin_lat  = np.sin(np.deg2rad(Lat))
    cos_lat  = np.sqrt(1.0-sin_lat**2)
    hr       = np.linspace(-np.pi+0.5*2.0*np.pi/NSTP, np.pi-0.5*2.0*np.pi/NSTP, NSTP)
    cos_hr   = np.cos(hr)

    cos_solarINC   = sin_lat*sin_dec + cos_lat*cos_dec*cos_hr
    sin_solarINC   = np.sqrt(1.0-cos_solarINC**2)
    solarAZ  = np.arccos( (sin_dec-sin_lat*cos_solarINC)/(cos_lat*sin_solarINC) )
    solarAZ  = np.where(hr>0, 2.0*np.pi - solarAZ, solarAZ)   

    solarINC = np.arccos(cos_solarINC)

    return solarAZ, solarINC





def MakeFluxLookupTable(savename, slope, aspect):
    # Calculates the direct and diffuse fluxes at many solar incidence and azimuth angles for a specific slope/aspect
    # Requires the standard DISORT results pickle file

    with open(savename, 'rb') as f:
        data = pickle.load(f)
    if isinstance(data, dict):
        for key in data:
            globals()[key] = data[key]
    del f,key,data

    NUMU = uu_surf_inc.shape[0]
    NAZ  = uu_surf_inc.shape[1]
    NINC = uu_surf_inc.shape[2]
    NALB = uu_surf_inc.shape[3]

    UMU_up = 0.5*(np.roll(UMU, 1,axis=0) + UMU)
    UMU_dn = 0.5*(np.roll(UMU,-1,axis=0) + UMU)
    UMU_up[0]  = -1.0
    UMU_dn[-1] = 1.0
    inc_flat = np.zeros([NUMU,NAZ])
    inc_up   = np.zeros([NUMU,NAZ])
    inc_dn   = np.zeros([NUMU,NAZ])
    for i in range(NAZ):
        inc_flat[:,i] = np.arccos(-UMU)
        inc_up[:,i]   = np.arccos(-UMU_up)
        inc_dn[:,i]   = np.arccos(-UMU_dn)

    phi_flat = np.zeros([NUMU,NAZ])
    for i in range(NUMU):
        phi_flat[i,:] = np.deg2rad(PHI)

    dphi = 2.0*np.pi/NAZ
    dinc = np.abs(inc_up-inc_dn)
    sa   = dphi*dinc*0.5*(np.sin(inc_up) + np.sin(inc_dn))

    cosi = np.cos(slope)*np.cos(inc_flat) + np.sin(slope)*np.sin(inc_flat)*np.cos(phi_flat-aspect)
    cosi = np.where(cosi < 0.0, 0.0, cosi)

    cosisa = cosi*sa

    direct  = np.zeros([360,NINC])
    diffuse = np.zeros([360,NINC,NALB])
    atmless = np.zeros([360,NINC])

    solarAZ_i       = np.linspace(0,359,360) * np.pi/180.0
    delta_solarAZ   = np.round((solarAZ_i*180.0/np.pi - PHI0) / (360.0/NAZ)).astype(int)        # Nearest neighbor azimuthal angle

    for k,solarAZ in enumerate(solarAZ_i): 
        uu_AZ_corrected = np.roll(uu_surf_inc, delta_solarAZ[k], axis=1)                                        # Roll the sky radiance array in azimuth to match the solar azimuth
        cosi_direct     = np.cos(slope)*UMU0_i + np.sin(slope)*np.sqrt(1.0-UMU0_i**2)*np.cos(solarAZ - aspect)  # Evaluate incidence angles on sloping surface for that solar azimuth
        cosi_direct     = np.where(cosi_direct < 0.0, 0.0, cosi_direct)                                         # zero out any solar positions that can't be seen by the sloping facet
        direct[k,:]     = rfldir_inc[:,0] * cosi_direct / UMU0_i                                                # Adjust all direct fluxes for updated incidence angles
        atmless[k,:]    = totalvisflux * cosi_direct                                                            # Figure out the direct flux with no atmosphere involved
               
        for i in range(0,NALB):
            for j in range(0,NINC):
                diffuse[k,j,i] = np.sum(uu_AZ_corrected[:,:,j,i] * cosisa) # Scale the radiances by cos(i) and solid-angle to get fluxes and sum over the sky

    sname2 = savename[:-4]+'_S'+"{:02d}".format(np.round(slope*180.0/np.pi).astype(int))+'A'+"{:03d}".format(np.round(aspect*180.0/np.pi).astype(int))+'.pkl'

    with open(sname2, 'wb') as f:
        pickle.dump({
            'uu_surf_inc': uu_surf_inc,
            'rfldir_inc': rfldir_inc,
            'rfldn_inc': rfldn_inc,
            'flup_inc': flup_inc,
            'PHI': PHI,
            'UMU': UMU,
            'UMU0_i': UMU0_i,
            'PHI0': PHI0,
            'ALB_i': ALB_i,
            'totalvisflux': totalvisflux,
            'slope': slope,
            'aspect': aspect,
            'solarAZ_i': solarAZ_i,
            'direct': direct,
            'diffuse': diffuse,
            'atmless': atmless
        }, f)

    return direct, diffuse, atmless



