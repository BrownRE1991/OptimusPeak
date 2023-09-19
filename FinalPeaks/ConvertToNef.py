import sys
import os
import numpy
import math

def fitness(model):
    return(model[6])

dim = int(sys.argv[1])
a1 = sys.argv[2]
a2 = sys.argv[3]
#print(len(sys.argv))
if(dim > 2):
    if(len(sys.argv) > 4):
        a3 = sys.argv[4]
    else:
        print("Nucleus 3 missing")
if(dim > 3):
    if(len(sys.argv) > 5):
        a4 = sys.argv[5]
    else:
        print("Nucleus 4 missing")
        
atom1 = ""
atom2 = ""
atom3 = ""
atom4 = ""
        
if(a1 == "C"):
    atom1 = "13C"
if(a1 == "N"):
    atom1 = "15N"
if(a1 == "H"):
    atom1 = "1H"
if(a2 == "C"):
    atom2 = "13C"
if(a2 == "N"):
    atom2 = "15N"
if(a2 == "H"):
    atom2 = "1H"
if(dim > 2):
    if(a3 == "C"):
        atom3 = "13C"
    if(a3 == "N"):
        atom3 = "15N"
    if(a3 == "H"):
        atom3 = "1H"
if(dim > 3):
    if(a4 == "C"):
        atom4 = "13C"
    if(a4 == "N"):
        atom4 = "15N"
    if(a4 == "H"):
        atom4 = "1H"

        
files = os.listdir("./Models_Decon")
models = []
found = []
count = 0

for x in files:
    if(x != ".DS_Store"):
        starting = "./Models_Decon/" + x
        #print(starting)
        f = open(starting, "r")
        stuff = numpy.array(f.readlines())
        f.close()
        lines = []
        for x in stuff:
            lines.append((x.strip('\n')).split(","))

        dataset = numpy.array(lines)
        for x in dataset:
            models.append(x)
        
        
        
openfile = open("Results2.nef", "w")
toPrint = "data_default\n\n\n   save_nef_nmr_meta_data\n\n      _nef_nmr_meta_data.sf_category      nef_nmr_meta_data\n      _nef_nmr_meta_data.sf_framecode     nef_nmr_meta_data\n      _nef_nmr_meta_data.format_name      nmr_exchange_format\n      _nef_nmr_meta_data.format_version   1.1\n      _nef_nmr_meta_data.program_name     AnalysisAssign\n      _nef_nmr_meta_data.program_version  3.1.0\n      _nef_nmr_meta_data.creation_date    2022-11-08T16:01:48.432332\n      _nef_nmr_meta_data.uuid             AnalysisAssign-2022-11-08T16:01:48.432332-288545018\n      _nef_nmr_meta_data.coordinate_file_name  .\n\n      loop_\n         _nef_program_script.program_name\n         _nef_program_script.script_name\n         _nef_program_script.script\n\n         CcpNmr  exportProject  .\n      stop_\n   save_"
toPrint = toPrint + "\n\n\n   save_nef_molecular_system\n\n      _nef_molecular_system.sf_category   nef_molecular_system\n      _nef_molecular_system.sf_framecode  nef_molecular_system\n   save_\n"
#toPrint = toPrint + "\n\n   save_nef_chemical_shift_list_default\n\n      _nef_chemical_shift_list.sf_category   nef_chemical_shift_list\n      _nef_chemical_shift_list.sf_framecode  nef_chemical_shift_list_default\n      _nef_chemical_shift_list.ccpn_serial   1\n      _nef_chemical_shift_list.ccpn_auto_update  true\n      _nef_chemical_shift_list.ccpn_is_simulated  false\n      _nef_chemical_shift_list.ccpn_comment  ''\n   save_\n"
toPrint = toPrint + "\n\n   save_nef_nmr_spectrum_XP_AQP1UCN_2DCC_HORROR7ms_H2O_20181207_80ppm_GM80_80`137`\n\n      _nef_nmr_spectrum.sf_category                   nef_nmr_spectrum\n      _nef_nmr_spectrum.sf_framecode                  nef_nmr_spectrum_XP_AQP1UCN_2DCC_HORROR7ms_H2O_20181207_80ppm_GM80_80`137`\n      _nef_nmr_spectrum.num_dimensions                2\n      _nef_nmr_spectrum.chemical_shift_list           nef_chemical_shift_list_default\n      _nef_nmr_spectrum.experiment_classification     .\n      _nef_nmr_spectrum.experiment_type               XP_AQP1UCN_2DCC_HORROR7ms_H2O_20181207_80ppm_GM80_80\n      _nef_nmr_spectrum.ccpn_reference_experiment_dimensions  '(None, None)'\n      _nef_nmr_spectrum.ccpn_positive_contour_count   10\n      _nef_nmr_spectrum.ccpn_positive_contour_base    137545.0177\n      _nef_nmr_spectrum.ccpn_positive_contour_factor  1.41\n      _nef_nmr_spectrum.ccpn_positive_contour_colour  '#008080'\n      _nef_nmr_spectrum.ccpn_negative_contour_count   10\n      _nef_nmr_spectrum.ccpn_negative_contour_base    -137545.0177\n      _nef_nmr_spectrum.ccpn_negative_contour_factor  1.41\n"
toPrint = toPrint + "      _nef_nmr_spectrum.ccpn_negative_contour_colour  '#DA70D6'\n      _nef_nmr_spectrum.ccpn_slice_colour             '#008080'\n      _nef_nmr_spectrum.ccpn_spectrum_scale           1\n      _nef_nmr_spectrum.ccpn_spinning_rate            .\n      _nef_nmr_spectrum.ccpn_spectrum_comment         ''\n      _nef_nmr_spectrum.ccpn_spectrum_file_path       /Users/pixie/Desktop/OptmizationExperiments/GoodnessOfFitTests/Threshold/XP_AQP1UCN_2DCC_HORROR7ms_H2O_20181207_80ppm_GM80_80.ft2\n      _nef_nmr_spectrum.ccpn_file_type                NMRPipe\n      _nef_nmr_spectrum.ccpn_file_scale_factor        1\n      _nef_nmr_spectrum.ccpn_sample                   .\n      _nef_nmr_spectrum.ccpn_peaklist_serial          1\n      _nef_nmr_spectrum.ccpn_peaklist_comment         ''\n      _nef_nmr_spectrum.ccpn_peaklist_name            .\n      _nef_nmr_spectrum.ccpn_peaklist_is_simulated    false\n      _nef_nmr_spectrum.ccpn_peaklist_symbol_colour   '#7a7a7a'\n      _nef_nmr_spectrum.ccpn_peaklist_symbol_style    cross\n      _nef_nmr_spectrum.ccpn_peaklist_text_colour     '#7a7a7a'\n"
toPrint = toPrint + "      _nef_nmr_spectrum.ccpn_file_header_size         512\n      _nef_nmr_spectrum.ccpn_file_number_type         float\n      _nef_nmr_spectrum.ccpn_file_is_big_endian       false\n      _nef_nmr_spectrum.ccpn_file_byte_number         4\n      _nef_nmr_spectrum.ccpn_file_block_header_size   0\n\n      loop_\n         _nef_spectrum_dimension.dimension_id\n         _nef_spectrum_dimension.axis_unit\n         _nef_spectrum_dimension.axis_code\n         _nef_spectrum_dimension.spectrometer_frequency\n         _nef_spectrum_dimension.spectral_width\n         _nef_spectrum_dimension.value_first_point\n         _nef_spectrum_dimension.folding\n         _nef_spectrum_dimension.absolute_peak_positions\n         _nef_spectrum_dimension.is_acquisition\n         _nef_spectrum_dimension.ccpn_axis_code\n\n         1  ppm"
"  " + atom1 + "   499.8859863  3.422740928  9.805181798  circular  true  true   " + a1 + "\n         2  ppm  " + atom2 + "  50.6590004   29.60974335  133.2798701  circular  true  false  "+ atom2 + "\n      stop_\n\n"
toPrint = toPrint + "      loop_\n         _ccpn_spectrum_dimension.dimension_id\n         _ccpn_spectrum_dimension.point_count\n         _ccpn_spectrum_dimension.reference_point\n         _ccpn_spectrum_dimension.total_point_count\n         _ccpn_spectrum_dimension.assignment_tolerance\n         _ccpn_spectrum_dimension.lower_aliasing_limit\n         _ccpn_spectrum_dimension.higher_aliasing_limit\n         _ccpn_spectrum_dimension.measurement_type\n         _ccpn_spectrum_dimension.phase_0\n         _ccpn_spectrum_dimension.phase_1\n         _ccpn_spectrum_dimension.window_function\n         _ccpn_spectrum_dimension.lorentzian_broadening\n         _ccpn_spectrum_dimension.gaussian_broadening\n         _ccpn_spectrum_dimension.sine_window_shift\n\n         1  292  .  .  0.03  6.388301728  9.811042656  Shift  -67.59999847  -36.40000153  .  .  .  .\n         2  512  .  .  0.4   103.6990426  133.3087859  Shift  0             0             .  .  .  .\n      stop_\n\n      loop_\n"
toPrint = toPrint + "         _nef_spectrum_dimension_transfer.dimension_1\n         _nef_spectrum_dimension_transfer.dimension_2\n         _nef_spectrum_dimension_transfer.transfer_type\n         _nef_spectrum_dimension_transfer.is_indirect\n\n      stop_\n\n      loop_\n         _nef_peak.index\n         _nef_peak.peak_id\n         _nef_peak.volume\n         _nef_peak.volume_uncertainty\n         _nef_peak.height\n         _nef_peak.height_uncertainty\n         _nef_peak.position_1\n         _nef_peak.position_uncertainty_1\n         _nef_peak.position_2\n         _nef_peak.position_uncertainty_2\n         _nef_peak.chain_code_1\n         _nef_peak.sequence_code_1\n         _nef_peak.residue_name_1\n         _nef_peak.atom_name_1\n         _nef_peak.chain_code_2\n         _nef_peak.sequence_code_2\n         _nef_peak.residue_name_2\n         _nef_peak.atom_name_2\n         _nef_peak.ccpn_figure_of_merit\n         _nef_peak.ccpn_linked_integral\n         _nef_peak.ccpn_annotation\n         _nef_peak.ccpn_comment\n         _nef_peak.ccpn_peak_list_serial\n"

openfile.write(toPrint)

         # 1  1  .  .  2126418.75   .  7.731134715  .  111.3880317  .  .  .  .  .  .  .  .  .  1  .  .  .  1
         # 2  2  .  .  1584258.75   .  8.048839614  .  111.9611235  .  .  .  .  .  .  .  .  .  1  .  .  .  1
         # 3  3  .  .  1632458.75   .  8.171685508  .  112.5819729  .  .  .  .  .  .  .  .  .  1  .  .  .  1
         # 4  4  .  .  1390872.875  .  8.548695321  .  112.2476694  .  .  .  .  .  .  .  .  .  1  .  .  .  1
         # 5  5  .  .  1854757.125  .  8.616472366  .  111.2447588  .  .  .  .  .  .  .  .  .  1  .  .  .  1
         # 6  6  .  .  2380915.75   .  8.743554326  .  110.8626976  .  .  .  .  .  .  .  .  .  1  .  .  .  1
         # 7  7  .  .  1507927.25   .  8.463974015  .  110.2418482  .  .  .  .  .  .  .  .  .  1  .  .  .  1


#add a \n at beginning of line instead of end
models.sort(key=fitness)
#print(models)
toPrint2 = ""
for x in range(1, len(models)+1):
    if(dim == 2):
        toPrint2 = "\n         " + str(x) + "  " + str(x) + "  " + str(float(models[x-1][7])) + "  .  " + str(round(float(models[x-1][4]), 3)) + "  .  " + models[x-1][1] + "  .  " + models[x-1][0] + "  .  .  .  .  .  .  .  .  .  1  .  .  .  1"
    if(dim == 3):
        toPrint2 = "\n         " + str(x) + "  " + str(x) + "  " + str(float(models[x-1][9])) + "  .  " + str(round(float(models[x-1][6]), 3)) + "  .  " + models[x-1][2] + "  .  " + models[x-1][1] +  "  .  " + models[x-1][0] + "  .  .  .  .  .  .  .  .  .  .  .  .  .  1  .  .  .  1"
#     if(x < 10):
#         toPrint2 = "     "
#     if((x > 9) and (x < 100)):
#         toPrint2 = "    "
#     if(x > 99 | x < 1000):
#         toPrint2 = "   "
#     if(len(str(round(float(models[x-1][1]), 3))) == 3):
#         toPrint2 = toPrint2 + str(x) + "   " + str(round(float(models[x-1][1]), 3)) + str(0) + str(0)
#     if(len(str(round(float(models[x-1][1]), 3))) == 4):
#         toPrint2 = toPrint2 + str(x) + "   " + str(round(float(models[x-1][1]), 3)) + str(0)
#     if(len(str(round(float(models[x-1][1]), 3))) == 5):
#         toPrint2 = toPrint2 + str(x) + "   " + str(round(float(models[x-1][1]), 3))
#     if(len(str(round(float(models[x-1][1]), 3))) == 6):
#         toPrint2 = toPrint2 + str(x) + "  " + str(round(float(models[x-1][1]), 3))
#     if(len(str(round(float(models[x-1][1]), 3))) == 7):
#         toPrint2 = toPrint2 + str(x) + " " + str(round(float(models[x-1][1]), 3))
#
#     if(len(str(round(float(models[x-1][0]), 3))) == 6):
#         toPrint2 = toPrint2 + " " +  str(round(float(models[x-1][0]), 3)) + str(0)
#     else:
#         toPrint2 = toPrint2 + " " +  str(round(float(models[x-1][0]), 3))
#     toPrint2 = toPrint2 + " " +  "2 U          0.000e+00  0.00e+00 -   0 0 0 0"
#     toPrint2 = toPrint2 + "\n";
    #print(toPrint2)
    openfile.write(toPrint2)
    
toPrint = "      stop_\n\n      loop_\n         _ccpn_spectrum_reference_substances.serial\n         _ccpn_spectrum_reference_substances.name\n         _ccpn_spectrum_reference_substances.labelling\n\n      stop_\n   save_\n\n\n   save_ccpn_assignment\n\n      _ccpn_assignment.sf_category   ccpn_assignment\n      _ccpn_assignment.sf_framecode  ccpn_assignment\n\n      loop_\n         _nmr_chain.short_name\n         _nmr_chain.serial\n         _nmr_chain.label\n         _nmr_chain.is_connected\n         _nmr_chain.comment\n\n"
toPrint = toPrint + "         @-  1  @-  false  'Default NmrChain, used for ResonanceGroups not in other chains. Cannot be deleted or renamed.'\n      stop_\n\n      loop_\n         _nmr_residue.chain_code\n         _nmr_residue.sequence_code\n         _nmr_residue.residue_name\n         _nmr_residue.serial\n         _nmr_residue.comment\n\n      stop_\n\n      loop_\n         _nmr_atom.chain_code\n         _nmr_atom.sequence_code\n         _nmr_atom.serial\n         _nmr_atom.name\n         _nmr_atom.isotope_code\n         _nmr_atom.comment\n\n      stop_\n   save_\n"
toPrint = toPrint + "\n\n   save_ccpn_additional_data\n\n      _ccpn_additional_data.sf_category   ccpn_additional_data\n      _ccpn_additional_data.sf_framecode  ccpn_additional_data\n\n      loop_\n         _ccpn_internal_data.ccpn_object_pid\n         _ccpn_internal_data.internal_data_string\n\n         PeakList:sbw_15Nsbw_hsqc.1             \n;{\n  \"_ccpNmrV3internal\": {\n    \"meritEnabled\": false\n  }\n}\n;\n         Spectrum:sbw_15Nsbw_hsqc               \n;{\n  \"_ccpNmrV3internal\": {\n    \"_dataStore\": \"[\\n  [\\n    \\\"_metadata\\\",\\n    {\\n      \\\"jsonVersion\\\": 3.0,\\n      \\\"className\\\": \\\"DataStore\\\",\\n      \\\"classVersion\\\": 1.0,\\n      \\\"user\\\": \\\"pixie\\\",\\n      \\\"lastPath\\\": \\\"undefined\\\",\\n      \\\"timestamp\\\": \\\"Tue Nov  8 16:00:06 2022\\\"\\n    }\\n  ],\\n  [\\n    \\\"_path\\\",\\n    \\\"/Users/pixie/Desktop/OptmizationExperiments/GoodnessOfFitTests/Threshold/sbw_15Nsbw_hsqc.pipe\\\"\\n  ],\\n  [\\n    \\\"dataFormat\\\",\\n    \\\"NMRPipe\\\"\n  ],\\n  [\\n    \\\"useBuffer\\\",\\n    false\\n  ],\\n  [\\n    \\\"apiDataStoreName\\\",\\n    null\\n  ],\\n  [\\n    \\\"apiDataStoreDir\\\",\\n    null\\n  ],\\n  [\\n    \\\"apiDataStorePath\\\",\\n    null\\n  ],\\n  [\\n    \\\"pathRedirections\\\",\\n    {\\n      \\\"$DATA\\\": \\\"/Users/pixie\\\",\\n      \\\"$ALONGSIDE\\\": \\\"/var/folders/0x/0l2f0y556yb8y93bz73ml9vc0000gn/T/CcpnProject_5ipvu0ul.ccpn\\\",\\n      \\\"$INSIDE\\\": \\\"/var/folders/0x/0l2f0y556yb8y93bz73ml9vc0000gn/T/CcpnProject_5ipvu0ul.ccpn/default\\\"\\n    }\\n  ]\\n]\",    \"_peakPickerParameters\": \"[\n  [\n    \"_metadata\",\n    {\n      \"jsonVersion\": 3.0,\n      \"className\": \"PeakPickerNd\",\n      \"classVersion\": 1.0,\n      \"user\": \"pixie\",\n      \"lastPath\": \"undefined\",\n      \"timestamp\": \"Tue Nov  8 16:00:06 2022\"\n    }\n  ],\n  [\n    \"dimensionCount\",\n    2\n  ],\n  [\n    \"pointExtension\",\n    1\n  ],\n  [\n    \"autoFit\",\n    false\n  ],\n  [\n    \"dropFactor\",\n    0.1\n  ],\n  [\n    \"fitMethod\",\n    \"parabolic\"\n  ],\n  [\n    \"positiveThreshold\",\n    10000.0\n  ],\n  [\n    \"negativeThreshold\",\n    -10000.0\n  ],\n  [\n    \"noise\",\n    null\n  ],\n  [\n    \"minimumLineWidth\",\n    []\n  ],\n  [\n    \"checkAllAdjacent\",\n    true\n  ],\n  [\n    \"singularMode\",\n    true\n  ],\n  [\n    \"halfBoxFindPeaksWidth\",\n    4\n  ],\n  [\n    \"halfBoxSearchWidth\",\n    4\n  ],\n  [\n    \"halfBoxFitWidth\",\n    4\n  ],\n  [\n    \"searchBoxDoFit\",\n    true\n  ],\n  [\n    \"setLineWidths\",\n    true\n  ],\n  [\n    \"searchBoxMode\",\n    true\n  ],\n  [\n    \"searchBoxWidths\",\n    {}\n  ]\n]\",\n    \"_preferredAxisOrdering\": [\n      0,\n      1\n    ],\n    \"negativeNoiseLevel\": -97549.6579384697\n  }\n}\n;\n         SpectrumReference:sbw_15Nsbw_hsqc.1.1  \n;{\n  \"_ccpNmrV3internal\": {\n    \"dimensionType\": \"Frequency\",\n    \"foldingMode\": \"circular\"\n  }\n}\n;\n         SpectrumReference:sbw_15Nsbw_hsqc.2.1  \n;{\n  \"_ccpNmrV3internal\": {\n    \"dimensionType\": \"Frequency\",\n    \"foldingMode\": \"circular\"\n  }\n}\n;\n      stop_\n   save_\n\n\n"
toPrint = toPrint + "# End of data_default\n"
    
#so far leaves linewidths and volumes blank. Can add in later if I feel like it.
openfile.write(toPrint)
openfile.close()