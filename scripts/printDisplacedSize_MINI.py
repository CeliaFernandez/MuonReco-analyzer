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


### Load files
_path = "/eos/user/f/fernance/Muon-POG/displacedCollection-implementation/11834.21/step4_220503_ptfinal/"
_files = []
for name in os.listdir(_path):

    _files.append(_path + name)


### Branch definition
branches = []
# muons
branches.append( ["patMuons_slimmedDisplacedMuons__PAT.", 0.0, 0.0] )
# segments
branches.append( ["CSCDetIdCSCSegmentsOwnedRangeMap_slimmedDisplacedMuons__PAT.", 0.0, 0.0] )
branches.append( ["DTChamberIdDTRecSegment4DsOwnedRangeMap_slimmedDisplacedMuons__PAT.", 0.0, 0.0] )
# Muon reduced track extras
branches.append( ["recoTrackExtras_slimmedDisplacedMuonTrackExtras__PAT.", 0.0, 0.0] )
branches.append( ["SiPixelClusteredmNewDetSetVector_slimmedDisplacedMuonTrackExtras__PAT.", 0.0, 0.0] )
branches.append( ["SiStripClusteredmNewDetSetVector_slimmedDisplacedMuonTrackExtras__PAT.", 0.0, 0.0] )
branches.append( ["TrackingRecHitsOwned_slimmedDisplacedMuonTrackExtras__PAT.", 0.0, 0.0] )

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

_out = open('pat_info.txt', 'w')
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

