import ROOT as r
import os
from ROOT import gROOT

gROOT.LoadMacro('printSizes.C')
gROOT.SetBatch(1)

verbose = True

def incolumn(text, width = 20):

    n = len(text)
    text2 = text + (width-n)*" " 
    return text2


### Load tree
"""
_file = r.TFile("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO/output_14.root")
_tree = _file.Get("Events")
tree_size = r.sizeOnDisk(_tree)
print(">>> Tree with {0} entries and with size {1} bytes".format(_tree.GetEntries(), tree_size))
"""
"""
_files = []
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_0.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_2.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_3.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_4.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_5.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_6.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_7.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_8.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_9.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_10.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_11.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_12.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_13.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_14.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_15.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_17.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_18.root")
_files.append("/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/RECO_new/output_19.root")
"""

_path = "/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/DisplacedSUSY_HS/RECO_filter1_v2/"
_files = []
for name in os.listdir(_path):

    _files.append(_path + name)

#tree_size = r.sizeOnDisk(_tree)
#print(">>> Tree with {0} entries and with size {1} bytes".format(_tree.GetEntries(), tree_size))


### Branch definition
branches = []
# muons
branches.append( ["recoMuons_displacedMuons__RECO.", 0.0, 0.0] )
# value maps (PF iso)
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueCharged03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueCharged04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueChargedAll03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueChargedAll04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueGamma03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueGamma04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueGammaHighThreshold03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueGammaHighThreshold04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueNeutral03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueNeutral04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueNeutralHighThreshold03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValueNeutralHighThreshold04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValuePU03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFMeanDRIsoValuePU04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueCharged03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueCharged04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueChargedAll03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueChargedAll04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueGamma03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueGamma04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueGammaHighThreshold03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueGammaHighThreshold04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueNeutral03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueNeutral04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueNeutralHighThreshold03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValueNeutralHighThreshold04_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValuePU03_RECO.", 0.0, 0.0] )
branches.append( ["doubleedmValueMap_displacedMuons_dispMuPFSumDRIsoValuePU04_RECO.", 0.0, 0.0] )
# time branches
branches.append( ["recoMuonTimeExtraedmValueMap_displacedMuons_combined_RECO.", 0.0, 0.0] )
branches.append( ["recoMuonTimeExtraedmValueMap_displacedMuons_csc_RECO.", 0.0, 0.0] )
branches.append( ["recoMuonTimeExtraedmValueMap_displacedMuons_dt_RECO.", 0.0, 0.0] )
# Muon reduced track extras
branches.append( ["recoTrackExtras_displacedMuonReducedTrackExtras__RECO.", 0.0, 0.0] )
branches.append( ["TrackingRecHitsOwned_displacedMuonReducedTrackExtras__RECO.", 0.0, 0.0] )
branches.append( ["recoTrackExtrasedmAssociation_displacedMuonReducedTrackExtras__RECO.", 0.0, 0.0] )
branches.append( ["SiPixelClusteredmNewDetSetVector_displacedMuonReducedTrackExtras__RECO.", 0.0, 0.0] )
branches.append( ["SiStripClusteredmNewDetSetVector_displacedMuonReducedTrackExtras__RECO.", 0.0, 0.0] )

if verbose: ">> Reading size of {0} branches".format(len(branches))

# Compute the size:
tree_size = 0.0
total_size = 0.0
nEntries = 0
for _file in _files:
    _tfile = r.TFile(_file)
    _tree = _tfile.Get("Events")
    nEntries += _tree.GetEntries()
    tree_size += r.sizeOnDisk(_tree)
    print("> Reading tree with entries: {0}".format(_tree.GetEntries()))
    for branch in branches:
        branch[1] += r.sizeOnDisk(_tree.GetBranch(branch[0]), True)
        #branch[2] += r.sizeOnDisk(_tree.GetBranch(branch[0]), True)/tree_size
        total_size += r.sizeOnDisk(_tree.GetBranch(branch[0]), True)
    _tfile.Close()

print(">>> Tree with {0} entries and with size {1} bytes".format(nEntries, tree_size))

_out = open('info.txt', 'w')
_out.write(140*"*"+"\n")
_out.write(">>> Tree with {0} entries and with size {1} bytes\n".format(nEntries, tree_size))
_out.write(140*"*"+"\n")
_out.write("--- Total space taken by displaced branches is {0} bytes which is the {1}% ---\n".format(total_size, total_size/tree_size*100.))
_out.write(140*"*"+"\n")
line = incolumn("Branch", 90)
line += "\t"
line += incolumn("Bytes")
line += "\t"
line += incolumn("Rel. Size in %")
line += "\n"
line.expandtabs(20)
_out.write(line)
_out.write(140*"*"+"\n")
for branch in branches:
    print(">> Branch {0} takes {1} bytes in disk, {2}% fraction of total size".format(branch[0], branch[1], branch[1]/tree_size*100.))
    line = incolumn(branch[0], 90)
    line += "\t"
    line += incolumn(str(branch[1]))
    line += "\t"
    line += incolumn(str(branch[1]/tree_size*100.))
    line += "\n"
    line.expandtabs(20)
    _out.write(line)
_out.write(140*"*"+"\n")
_out.close()

print("-- The total space taken by displaced branches is {0} bytes, fraction {1}% --".format(total_size, total_size/tree_size*100.))

