"""Script to create Hazts_<source_zone> files from GEM
Subduction zone report parameters
Currently hard coded to read the GEM subduction zone database
and a mapping to the GAR source zones
"""
import os, sys

pathname = '/short/n74/GAR_2015/Tsunami_v2/DATA/FOR_MODEL'
final_source_path = os.path.join(pathname, 'source_zones')
austptha_gem_map = '/short/n74/GAR_2015/Tsunami_v2/DATA/FOR_MODEL/ancillary/aust_gem_mapping.csv'
gem_tab_path = '/short/n74/GAR_2015/Tsunami_v2/DATA/FOR_MODEL/ancillary'

source_zones = os.walk(final_source_path).next()[1]

# Create dictionary of mapping                                                               
mapping = {}
spareslips = {}
sparebs = {}
sparemmaxs = {}
sparemmaxs_pref = {}
sparemmaxs_min = {}
f_in = open(austptha_gem_map, 'r')
header = f_in.readline()
for line in f_in.readlines():
    row = line.split(',')
    mapping[row[0]]=row[1].replace(' ','').replace('-','')
    spareslips[row[0]]=row[8]
    sparemmaxs[row[0]]=row[9]
    if row[9] != '':
        sparemmaxs_pref[row[0]]=str(float(row[9])-0.1)
        sparemmaxs_min[row[0]]=str(float(row[9])-0.2)
    sparebs[row[0]]=row[10].rstrip('\r\n')
print mapping
f_in.close()

# Create dictionary of gem slip rates
gem_slips = {}
gem_coupling_pref = {}
gem_slips_pref_coupling = {}
gem_loc_file = os.path.join(gem_tab_path, 'GEM_TABLE_LOCATION.txt')
f_in = open(gem_loc_file, 'r')
header = f_in.readline()
for line in f_in.readlines():
    row = line.split(' ')
    name = row[0] + row[1]                                                                   
    name = name.replace(' ','').replace('-','')
    slip_rate = max(float(row[5]), float(row[9]))
    gem_slips[name] = slip_rate
#    gem_coupling_pref[name] = float(row[13])
#    gem_slips_pref_coupling[name] = slip_rate*gem_coupling_pref[name]
    slip_diff = abs(float(row[5]) - float(row[9]))
    if slip_diff > 15.0:
        print 'More than 15 mm/yr difference in slip rates at either end of zone %s, difference = %.2f mm/yr' % (name, slip_diff)
f_in.close()

# Get preferred coupling value
gem_geom_file = os.path.join(gem_tab_path, 'GEM_TABLE_SLAB_GEOMETRY.txt')
f_in = open(gem_geom_file, 'r')
header = f_in.readline()
for line in f_in.readlines():
    row = line.split(' ')
    name = row[0] + row[1]
    name = name.replace(' ','').replace('-','')
    gem_coupling_pref[name] = float(row[13])                                                                                               
    gem_slips_pref_coupling[name] = gem_slips[name]*gem_coupling_pref[name]

# Create dictionary of GEM b values, slip-rates weighted by
# coupling coefficient and maximum magnitudes (Mmax)
gem_bs = {}
gem_mmaxs = {}
gem_mmaxs_pref = {}
gem_mmaxs_min = {}
gem_coupling_min = {}
gem_coupling_max = {}
gem_slips_min_coupling = {}
gem_slips_max_coupling = {}
gem_haz_file = os.path.join(gem_tab_path, 'GEM_TABLE_HAZ.txt')
f_in = open(gem_haz_file, 'r')
header = f_in.readline()
for line in f_in.readlines():
    row = line.split(' ')
    name = row[0] + row[1]
    name = name.replace(' ','').replace('-','')
    coupling_min = float(row[2])
    gem_coupling_min[name] = coupling_min
    gem_slips_min_coupling[name] = gem_slips[name]*coupling_min
    coupling_max = float(row[3])
    gem_coupling_max[name] = coupling_max
    gem_slips_max_coupling[name] = gem_slips[name]*coupling_max
    max_mag_pref = row[4]
    gem_mmaxs_pref[name] = max_mag_pref
    max_mag = row[6]
    gem_mmaxs[name]=max_mag
    max_mag_min = row[5]
    gem_mmaxs_min[name] = max_mag_min
    b_value = float(row[7])
    gem_bs[name]=b_value
f_in.close()

# Combine into one dictionary
slip_rates_min_coupling = {}
slip_rates_max_coupling = {}
slip_rates_pref_coupling = {}
b_values = {}
mmaxs_pref = {}
mmaxs = {}
mmaxs_min = {}
for key in mapping.keys():
    print key
    if mapping[key] == 'NULL':
        slip_rates_pref_coupling[key] = float(spareslips[key])*0.5
        slip_rates_min_coupling[key] = float(spareslips[key])*0.3 #Assume min coupling of 0.3
        slip_rates_max_coupling[key] = float(spareslips[key])*0.7 # Assume max coupling of 0.7
        b_values[key] = float(sparebs[key])
        mmaxs[key] = sparemmaxs[key]
        mmaxs_pref[key] = sparemmaxs_pref[key]
        mmaxs_min[key] = sparemmaxs_min[key]
    else:
        gem_name = mapping[key]
        try:
            slip_pref_coupling = gem_slips_pref_coupling[gem_name]
            slip_min_coupling = gem_slips_min_coupling[gem_name]
            slip_max_coupling = gem_slips_max_coupling[gem_name]
            b_value = gem_bs[gem_name]
            mmax = gem_mmaxs[gem_name]
            mmax_pref = gem_mmaxs_pref[gem_name]
            mmax_min = gem_mmaxs_min[gem_name]
        except KeyError:
            try:
                gem_name = gem_name + '-'
                slip_pref_coupling = gem_slips_pref_coupling[gem_name]
                slip_min_coupling = gem_slips_min_coupling[gem_name]
                slip_max_coupling =gem_slips_max_coupling[gem_name]
                b_value = gem_bs[gem_name]
                mmax = gem_mmaxs[gem_name]
                mmax_pref = gem_mmaxs_pref[gem_name]
                mmax_min = gem_mmaxs_min[gem_name]
            except KeyError:
                print 'No mapping for %s to %s' % (gem_name, key)
                print mapping[key], type(mapping[key])
                print gem_name, type(gem_name)
                sys.exit()
        slip_rates_pref_coupling[key] = slip_pref_coupling
        slip_rates_min_coupling[key] = slip_min_coupling
        slip_rates_max_coupling[key] = slip_max_coupling
        b_values[key] = b_value
        mmaxs[key] = mmax
        mmaxs_pref[key] = mmax_pref
        mmaxs_min[key] = mmax_min

# Write out to Hazts file
for source_zone in source_zones:    
#    if source_zone != 'north_african':# and source_zone != 'cadiz' and source_zone != 'scaribbean' and source_zone != 'hellenic' and source_zone != 'antilles' and source_zone != 'cyprus' and source_zone != 'hispaniola' and source_zone != 'north_africa':  
 #       continue
    Hazts_filebasename = 'Hazts_' + source_zone + '.flt'
    Hazts_filename = os.path.join(final_source_path, source_zone, Hazts_filebasename)
    print Hazts_filename
    f_out = open(Hazts_filename,'w')
#b_value
#source_zone
#max_mag


    lines = range(35)

    lines[0] = '1 iCoor (0=(x,y), 1 =(long,lat))      ***started with slc-flt-sr09-dd.dat'
    lines[1] = '1 flts (no. of faults)'
    lines[2] = '   ' + source_zone + '  (fault system name)'
    lines[3] = '   1.0   Prob Activity'
    lines[4] = '   1     nSeg model'
    lines[5] = '      1 1.0         nflt segments, wt, Segmented model' # Assume no segments for now
    lines[6] = '         ' + source_zone
    lines[7] = '         10 1 0.1 1 0  source type, atten type, sampleStep (km), fltdir, syncron rup (1=yes)'
    lines[8] = '         1.0           aleatory seg wt'
    lines[9] = '         i_invall-' + source_zone
    lines[10] = '         0 0 0 0       along strike MIN MAX down dip MIN MAX'
    lines[11] = '         1             Number of dips'
    lines[12] = '         0             Dip variations'
    lines[13] = '         1.0           Weights for dip vars'
    lines[14] = '         1             Number of b-values'
    lines[15] = '         %.2f          b-values' % b_values[source_zone]
    lines[16] = '         1.0           Weights for b-values'
    lines[17] = '         3             Number of slip rates'
    lines[18] = '         %.2f %.2f %.2f         Slip rates (mm/yr)' % \
    (slip_rates_pref_coupling[source_zone], slip_rates_min_coupling[source_zone], \
         slip_rates_max_coupling[source_zone])
    lines[19] = '         0.6 0.2 0.2            Weights for slip rates'
    lines[20] = '         1             Number of recurrence models'
    lines[21] = '         1             Recur model, 0 = char, 1 = exp, 3 = max mag.'
    lines[22] = '         1.            Weights for recurrence models'
#    lines[23] = '         0.5 1.0 0.25  Delta_m1, delta_m2, for char mag rec model'
    lines[23] = '         1             Number of fault widths'
    lines[24] = '         50            Fault widths'
    lines[25] = '         1.0           Weights for fault widths'
    lines[26] = '         1             OverrideMagnitude'
    lines[27] = '         3             Number of maximum magnitudes'
    lines[28] = '         ' + mmaxs_pref[source_zone] + ' ' + mmaxs[source_zone] + ' ' + mmaxs_min[source_zone] + \
        '          Maximum magnitudes'
    lines[29] = '         0.7 0.2 0.1          Weights for maximum mags'
    lines[30] = '         7.8 0.1  5.0 3.0  1  1  0.0 (minmag, magstep, hxStep, hzStep, nRupArea, nRupWidth, minDepth)'
    lines[31] = '         1'
    lines[32] = '         -3.476  0.952  0.304             (rupture area consts a, b, sigma) Strasser et al 2010'
    lines[33] = '         -0.882  0.351  0.173             (rupture width consts a, b, sigma)'
    lines[34] = '         1    ftype (0 = s/s, 1 = rv, -1.0 = normal)'
    
    for line in lines:
#        print line
        f_out.write(line + '\n')
    f_out.close()











