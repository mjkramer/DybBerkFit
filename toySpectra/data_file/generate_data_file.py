#!/usr/bin/env python

#Feb. 21. 2013: removed absolute energy scale error, since it's double-counted with scinti_nl
#Feb. 26, 2013: removed accidental shape error, since this is precisely measured by data
#May  10, 2013: Add option to vary solar oscillation parameters
#June 26, 2013: Updated Li9 background treatment
#July 5, 2013: Updated amc background shape distortion. Also removed alpha-n shape distortion
base_filename = 'dyb_data_v1'

variations = [('randomizeSolarOscPars','solar_oscpars', '1'),
              ('randomizeReactorPower','reac_power', '1'),
              ('randomizeIavDistortion', 'iav', '1'),
              ('randomizeScintiNonLinear', 'scinti_nl', '1'),
              ('randomizeRelativeEnergyScale', 'rel_escale', '1'),
              ('randomizeResolution', 'resolutoin', '1'),
              ('randomizeCoreSpectra', 'core_spectra', '1'),
              ('randomizeDetectorEfficiency', 'det_eff', '1'),
              ('varyAccBg', 'vary_acc', '1'),
              ('varyAmcBg', 'vary_amc', '1'),
              ('varyFnBg', 'vary_fn', '1'),
              ('varyLi9Bg', 'vary_li9', '1'),
              ('varyAlnBg', 'vary_aln', '1' ),
              ('distortAmcBg', 'distort_amc', '0.15'),
              ('distortFnBg', 'distort_fn', '0.2'),
              ('distortLi9Bg', 'distort_li9', '../li9_spectrum/8he9li_distort_neutron100_alpha100_frac0.055_N250.root'),
              ('statisticalFluctuation', 'stat', '1' )]


base_file = open(base_filename+'_nominal.txt','r')
basefile_lines = base_file.readlines()

for par, opt, val in variations:
    output_file = open(base_filename+'_'+opt+'.txt','w')
    #for line in base_file.read().split('\n'):
    for line in basefile_lines:
        items = line.split()
        if len(items) > 1:
            if items[0] == par:
                output_file.write(par + ' ' + val + '\n')
            else:
                output_file.write(line)
        else:
            output_file.write(line)
    base_file.close()
    output_file.close()


output_file2 = open(base_filename+'_allsys.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        match = False
        for par, opt, val in variations:
            if items[0] == par and opt != 'stat':
                output_file2.write(par + ' ' + val + '\n')
                match = True
        if match == False:
            output_file2.write(line)
    else:
        output_file2.write(line)

        
output_file3 = open(base_filename+'_allsys_and_stat.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        match = False
        for par, opt, val in variations:
            if items[0] == par:
                output_file3.write(par + ' ' + val + '\n')
                match = True
        if match == False:
            output_file3.write(line)
    else:
        output_file3.write(line)

        

output_file4 = open(base_filename+'_bgsys.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        match = False
        for par, opt, val in variations:
            if 'vary' in items[0] or 'distort' in items[0]:
                if items[0] == par:
                    output_file4.write(par + ' ' + val + '\n')
                    match = True
        if match == False:
            output_file4.write(line)
    else:
        output_file4.write(line)

output_file5 = open(base_filename+'_sigsys.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        match = False
        for par, opt, val in variations:
            if 'randomize' in items[0]:
                if items[0] == par:
                    output_file5.write(par + ' ' + val + '\n')
                    match = True
        if match == False:
            output_file5.write(line)
    else:
        output_file5.write(line)

        
