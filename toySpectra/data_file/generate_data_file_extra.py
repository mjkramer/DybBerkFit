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
              ('distortLi9Bg', 'distort_li9', '../li9_spectrum/8he9li_distort_neutron100_alpha100_frac0.1_N250.root'),
              ('statisticalFluctuation', 'stat', '1' ),
              ('useBcwFluxUncertainty', 'bcwflux', '1' )]

reac_list = ["randomizeCoreSpectra",
             "randomizeReactorPower"]
             
det_list = ["randomizeIavDistortion",
            "randomizeScintiNonLinear",
            "randomizeRelativeEnergyScale",
            "randomizeResolution",
            "randomizeDetectorEfficiency"]

base_file = open(base_filename+'_nominal.txt','r')
basefile_lines = base_file.readlines()

output_file = open(base_filename+'_reactor_only.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        if items[0] in reac_list:
            for par, opt, val in variations:
                if par == items[0]:
                    output_file.write(par + ' ' + val + '\n')
        else:
            output_file.write(line)
    else:
        output_file.write(line)


output_file2 = open(base_filename+'_det_only.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        if items[0] in det_list:
            for par, opt, val in variations:
                if par == items[0]:
                    output_file2.write(par + ' ' + val + '\n')
        else:
            output_file2.write(line)
    else:
        output_file2.write(line)

                
output_file3 = open(base_filename+'_allsys_wo_reactor.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        if items[0] in reac_list:
            output_file3.write(line)
        else:
            match = False
            for par, opt, val in variations:
                if par == items[0] and opt != 'stat' and opt != 'bcwflux':
                    output_file3.write(par + ' ' + val + '\n')
                    match = True
            if match == False:
                output_file3.write(line)
                   
    else:
        output_file3.write(line)

output_file4 = open(base_filename+'_allsys_wo_det.txt','w')
for line in basefile_lines:
    items = line.split()
    if len(items) > 1:
        if items[0] in det_list:
            output_file4.write(line)
        else:
            match = False
            for par, opt, val in variations:
                if par == items[0] and opt != 'stat' and opt != 'bcwflux':
                    output_file4.write(par + ' ' + val + '\n')
                    match = True
            if match == False:
                output_file4.write(line)
                   
    else:
        output_file4.write(line)

# finally make data files with BCW flux


bcw_list = ['nominal',
            'allsys',
            'sigsys',
            'reactor_only',
            'allsys_wo_det']
# bcw_list = ['reactor_only',
#             'allsys_wo_det']

for bcw_opt in bcw_list:
    base_file_bcw = open(base_filename+'_'+bcw_opt+'.txt','r')
    basefile_bcw_lines = base_file_bcw.readlines()
    
    output_file_bcw = open(base_filename+'_'+bcw_opt+'_bcwflux.txt','w')
    for line in basefile_bcw_lines:
        # print line
        items = line.split()
        if len(items) > 1:
            if items[0] == 'useBcwFluxUncertainty':
                for par, opt, val in variations:
                    if par == items[0]:
                        output_file_bcw.write(par + ' ' + val + '\n')
            else:
                output_file_bcw.write(line)
        else:
            output_file_bcw.write(line)

    base_file_bcw.close()
    output_file_bcw.close()
