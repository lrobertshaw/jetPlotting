# Load libraries
print("Loading libraries")
import sys
import uproot
import awkward as ak
import numpy as np
import pickle
from collections import namedtuple
print("Libraries loaded")

print("Parsing cmd line args")
path, branch = str(sys.argv[1]), str(sys.argv[2])

print(f"Loading {branch} data")
readData = lambda path: uproot.open(path)['Events']
Candidate = namedtuple("Candidate", "pt eta phi")

data = readData(path=path)

def getJetMass(data, branch, event, r=0.8, pTcut = 0.0):
    jetProps = {}

    pt, eta, phi = data[branch+"Jets_pt"].array()[event], data[branch+"Jets_eta"].array()[event], data[branch+"Jets_phi"].array()[event]
    nJets = len(pt)
    axes = tuple( map( Candidate, pt, eta, phi ) )

    for jetIdx in range(nJets):
        if axes[jetIdx].pt < pTcut: break    # Exit out of loop if pt of current jet doesnt meet threshold

        # Get jet axes
        jetProps[f"jet{jetIdx}"] = {"axis":axes[jetIdx]}

        # get jet cands
        jetProps[f"jet{jetIdx}"]["cands"] = []            
        candIdx = 0
        while True:
            try:
                pt, eta, phi = data[f"{branch}Jets_dau{candIdx}_pt"].array()[event], data[f"{branch}Jets_dau{candIdx}_eta"].array()[event], data[f"{branch}Jets_dau{candIdx}_phi"].array()[event]
                jetProps[f"jet{jetIdx}"]["cands"].append( Candidate(pt = pt[jetIdx], eta = eta[jetIdx], phi = phi[jetIdx]) ) if pt[jetIdx] != -1 else None
                candIdx += 1
            except: break

        # get jet seed
        pts, etas, phis = zip(*jetProps[f"jet{jetIdx}"]["cands"])
        jetProps[f"jet{jetIdx}"]["seed"] = jetProps[f"jet{jetIdx}"]["cands"][pts.index(max(pts))]

        # Calculate jet mass
        axisPt, axisEta, axisPhi = jetProps[f"jet{jetIdx}"]["axis"].pt, jetProps[f"jet{jetIdx}"]["axis"].eta, jetProps[f"jet{jetIdx}"]["axis"].phi
        seedPt, seedEta, seedPhi = jetProps[f"jet{jetIdx}"]["seed"].pt, jetProps[f"jet{jetIdx}"]["seed"].eta, jetProps[f"jet{jetIdx}"]["seed"].phi

        p1_tot = 0
        p1x_tot = 0
        p1y_tot = 0
        p1z_tot = 0

        constituents = jetProps[f"jet{jetIdx}"]["cands"]
        for constituent in constituents:
            constitPt, constitEta, constitPhi = constituent.pt, constituent.eta, constituent.phi
            if constitPhi - seedPhi > r  :  constitPhi -= 2*np.pi #phiDist = constitPhi - axisPhi - 2*np.pi
            if constitPhi - seedPhi < -r :  constitPhi += 2*np.pi #phiDist = constitPhi - axisPhi + 2*np.pi

            etaDist, phiDist = constitEta - axisEta, constitPhi - axisPhi
            dist2 = etaDist**2 + phiDist**2

            p1_tot  += constitPt
            p1x_tot += constitPt * etaDist
            p1y_tot += constitPt * phiDist
            p1z_tot += constitPt * (1 - (dist2 / 2))

        jetProps[f"jet{jetIdx}"]["mass"] = p1_tot**2 - p1x_tot**2 - p1y_tot**2 - p1z_tot**2

    return jetProps


nEvents = len(data[f"nscPuppiJets"].array())
branch = "scPuppi"
datList = np.array([])
ptCut = 100.0
pTmass = namedtuple("pTmass", "pt mass")

print(f"Calculating jet mass for all jets with pT > {ptCut} in {nEvents} events")
for event in range(nEvents):
    for m in getJetMass(data=data, branch=branch, event=event, pTcut=ptCut).values():
        datList = np.append(datList, pTmass(pt = m["axis"].pt, mass = m["mass"]) )

    if event%1000==0 and event!=0: np.save(f"sigOutput/datListUpTo{event}.npy", datList)
    if event%500==0 and event!=0: print(f"Processed and saved up to {event} events")

print(f"Jet masses calculated, total masses: {len(datList)}")
print(f"Saving list of jet masses to datList100GeV.npy")
np.save("sigOutput/datList100GeV.npy", datList)