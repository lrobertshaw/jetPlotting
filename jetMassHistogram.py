# Load libraries
print("Loading libraries")
import sys
import uproot
import awkward as ak
import numpy as np
import pickle
import vector
import matplotlib.pyplot as plt
print("Libraries loaded")

print("Loading data")
global DATA, BRANCH
DATA = uproot.open( str(sys.argv[1]) )['Events']
BRANCH = str(sys.argv[2])
print("Data loaded")


def plotJetMasses(masses):
    #binning = np.linspace( 0, 150, 5 )
    plt.figure(figsize=(10,6))
    plt.hist(masses, bins=int(np.ceil(len(masses) / 5)), color="black")

    plt.axvline(125.4 / 2, label="H mass / 2", linestyle="--", color="red")
    plt.axvline(91.2 / 2, label="Z mass / 2", linestyle="--", color="green")
    plt.axvline(80.4 / 2, label="W mass / 2", linestyle="--", color="blue")

    plt.title(f"Histogram of jet masses. nJets = {len(masses)}")
    plt.xlabel("Mass (GeV)")
    plt.legend(loc="upper right")
    plt.savefig(f"jetMassPlots/hist{len(masses)}.png")


def getJetMassVector():

    countData = DATA[f"n{BRANCH}Jets"].array()
    nEvents = len(countData)
    masses = np.array([])
    print(f"Total number of events: {nEvents}")
    for event in range(nEvents):
        nJets = countData[event]
        for jetIdx in range(nJets):
            candIdx = 0
            jet_vector = vector.obj(pt=0, phi=0, eta=0, m=0)
            while True:
                try:
                    pt, phi, eta, mass = DATA[f"{BRANCH}Jets_dau{candIdx}_pt"].array()[event][jetIdx], DATA[f"{BRANCH}Jets_dau{candIdx}_phi"].array()[event][jetIdx], DATA[f"{BRANCH}Jets_dau{candIdx}_eta"].array()[event][jetIdx], DATA[f"{BRANCH}Jets_dau{candIdx}_mass"].array()[event][jetIdx]
                    if pt != -1:
                        jet_vector += vector.obj(pt = pt, phi = phi, eta = eta, mass = mass)
                    else:
                        break
                    candIdx += 1
                except: break
            masses = np.append(masses, jet_vector.m)

        if (event+1) % 1000 == 0:
            print(f"Processed {event+1} events")
        # if (event+1) % 250 == 0:
        #     plotJetMasses(masses=masses)
        #     np.save(f"jetMassData/data{len(masses)}.npy", masses)

    np.save(f"{BRANCH}Masses.npy", masses)
    return 0


print("Calculating jet masses")
getJetMassVector()


# def getJetMassCalc(DATA, BRANCH, event, r=0.8, pTcut = 0.0):
#     jetProps = {}

#     pt, phi, eta, mass = DATA[BRANCH+"Jets_pt"].array()[event], DATA[BRANCH+"Jets_phi"].array()[event], DATA[BRANCH+"Jets_eta"].array()[event], DATA[BRANCH+"Jets_mass"].array()[event]
#     nJets = len(pt)
#     axes = tuple( map( Candidate, pt, phi, eta, mass ) )

#     for jetIdx in range(nJets):
#         if axes[jetIdx].pt < pTcut: break    # Exit out of loop if pt of current jet doesnt meet threshold

#         # Get jet axes
#         jetProps[f"jet{jetIdx}"] = {"axis":axes[jetIdx]}

#         # get jet cands
#         jetProps[f"jet{jetIdx}"]["cands"] = []
#         candIdx = 0
#         while True:
#             try:
#                 pt, phi, eta, mass = DATA[f"{BRANCH}Jets_dau{candIdx}_pt"].array()[event], DATA[f"{BRANCH}Jets_dau{candIdx}_phi"].array()[event], DATA[f"{BRANCH}Jets_dau{candIdx}_eta"].array()[event], DATA[f"{BRANCH}Jets_dau{candIdx}_mass"].array()[event]
#                 jetProps[f"jet{jetIdx}"]["cands"].append( vector.obj(pt=pt[jetIdx], phi=phi[jetIdx], eta=eta[jetIdx], m=mass[jetIdx]) ) if pt[jetIdx] != -1 else None
#                 candIdx += 1
#             except: break

#         # get jet seed
#         pts, phis, etas, masses = zip(*jetProps[f"jet{jetIdx}"]["cands"])
#         jetProps[f"jet{jetIdx}"]["seed"] = jetProps[f"jet{jetIdx}"]["cands"][pts.index(max(pts))]

#         # Calculate jet mass
#         axisPt, axisEta, axisPhi = jetProps[f"jet{jetIdx}"]["axis"].pt, jetProps[f"jet{jetIdx}"]["axis"].eta, jetProps[f"jet{jetIdx}"]["axis"].phi
#         seedPt, seedEta, seedPhi = jetProps[f"jet{jetIdx}"]["seed"].pt, jetProps[f"jet{jetIdx}"]["seed"].eta, jetProps[f"jet{jetIdx}"]["seed"].phi

#         p1_tot = 0
#         p1x_tot = 0
#         p1y_tot = 0
#         p1z_tot = 0

#         constituents = jetProps[f"jet{jetIdx}"]["cands"]
#         for constituent in constituents:
#             constitPt, constitEta, constitPhi = constituent.pt, constituent.eta, constituent.phi
#             if constitPhi - seedPhi > r  :  constitPhi -= 2*np.pi #phiDist = constitPhi - axisPhi - 2*np.pi
#             if constitPhi - seedPhi < -r :  constitPhi += 2*np.pi #phiDist = constitPhi - axisPhi + 2*np.pi

#             etaDist, phiDist = constitEta - axisEta, constitPhi - axisPhi
#             dist2 = etaDist**2 + phiDist**2

#             p1_tot  += constitPt
#             p1x_tot += constitPt * etaDist
#             p1y_tot += constitPt * phiDist
#             p1z_tot += constitPt * (1 - (dist2 / 2))

#         jetProps[f"jet{jetIdx}"]["mass"] = p1_tot**2 - p1x_tot**2 - p1y_tot**2 - p1z_tot**2

#     return jetProps