#include <vector>
#include <list>
#include "utility.h"
#include "../CTab.h"
#include "../StdResH.h"
#include "../ResBlk.h"
#include "pdb++.h"
#include "../PDBrec.h"
#include "../AtomPositions.h"

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

void recordSkipInfo(bool skipH, std::vector<std::string>& fixNotes,
   const PDBrec& theHatom, const PDBrec& heavyAtom,
   std::list< std::pair<PDBrec*, Point3d> >& nearr, const char * msg);

//All defined in reduce.cpp
extern bool DoOnlyAltA; 
extern float NonMetalBumpBias;
extern float MetalBumpBias;
extern bool Verbose;

bool okToPlaceHydHere(const PDBrec& theHatom, const atomPlacementPlan& pp,
					  const PDBrec& a1, const PDBrec& a2, AtomPositions& xyz,
					  bool& doNotAdjustSC, std::vector<std::string>& fixNotes ) {

	std::list< std::pair<PDBrec*, Point3d> > emptySet;

	const Point3d& theHpos = theHatom.loc();

	if (pp.hasFeature(ROTATEFLAG) || pp.hasFeature(NH3FLAG)) {
		const PDBrec& heavyAtom = a1;
		if (heavyAtom.elem().atno() != 6) { // OH, SH and NH3, but not CH3
			const double halfbondlen = heavyAtom.elem().covRad();
			const double  maxbondlen = halfbondlen
				+ ElementInfo::StdElemTbl().maxCovalentRadius();
			std::list< std::pair<PDBrec*, Point3d> > sym_nearr_list = xyz.neighbors( heavyAtom.loc(),
				halfbondlen + 0.1,
				maxbondlen  + 0.25);
			std::list< std::pair<PDBrec*, Point3d> > c2batoms;
			int countOfBonds = 0;
			PDBrec* rec = NULL;
			for (std::list< std::pair<PDBrec*,Point3d> >::const_iterator nearr = sym_nearr_list.begin(); nearr != sym_nearr_list.end(); ++nearr) {
				rec = nearr->first;

				if (interactingConfs(heavyAtom, *rec, DoOnlyAltA)
					&& (! rec->elem().isHydrogen())
					&& (! rec->isWater()) ) {
					const double actual = distance2(heavyAtom.loc(),
						nearr->second);
					const double expected = heavyAtom.elem().covRad()
						+ rec->elem().covRad();
					if ((actual >= (expected - 0.55))
						&&  (actual <= (expected + 0.25))) {

						c2batoms.push_front(*nearr);

						if (++countOfBonds >= 2) {
							// we have an R-X-R or R-X-X-R bond

							if (pp.hasFeature(NH3FLAG)) {
								bool skipthisH = FALSE;
								PDBrec* nbs = NULL;
								for(std::list< std::pair<PDBrec*, Point3d> >::const_iterator it = c2batoms.begin(); it != c2batoms.end(); ++it) {
									nbs = it->first;
									const double hRdist = distance2(theHpos,
										it->second);
									const double bumpLim = nbs->elem().explRad()
										+ NonMetalBumpBias;
									if (hRdist < bumpLim) { skipthisH = TRUE; }

									if ((bumpLim - hRdist) < 0.02) {
										cerr << "WARNING:" << heavyAtom.recName()
											<< " H atom too close to "
											<< nbs->recName() << " by "
											<< (bumpLim - hRdist) << "A" << endl;
										cerr << "        if warning inappropriate, adjust the cuttoffs using "
											<< "-METALBump or -NONMETALBump." << endl;
									}
								}

								if (skipthisH) {
									recordSkipInfo(TRUE, fixNotes, theHatom, heavyAtom, c2batoms, "(NH2R)");
									return FALSE; // don't add hyd.
								}
							}
							else { // C-S-S-C, C-O-C, C-O-P, etc.
								//recordSkipInfo(TRUE, fixNotes, theHatom, heavyAtom, c2batoms, "(RXR or RXXR)"); //2.21.mod_dcr
								return FALSE; // don't add hyd.
							}
						}
					}
				}
			}
		}
	}

	// *** special case for N/Q metal coordination ***
	if (pp.hasFeature(BONDBUMPFLAG)
		&& strstr(":ASN:asn:GLN:gln:ASX:asx:GLX:glx:", theHatom.resname())) {
		const PDBrec& nqoxygen = a2;
		const double halfbondlen = nqoxygen.elem().covRad();
		const double  maxbondlen = halfbondlen
			+ ElementInfo::StdElemTbl().maxCovalentRadius();


		std::list< std::pair<PDBrec*, Point3d> > sym_nearr_list = xyz.neighbors( nqoxygen.loc(),
			halfbondlen + 0.1,
			maxbondlen  + 0.25);

		PDBrec* rec = NULL;
		for (std::list< std::pair<PDBrec*, Point3d> >::const_iterator nearr = sym_nearr_list.begin(); nearr != sym_nearr_list.end(); ++nearr) {
			rec = nearr->first;
			if (interactingConfs(nqoxygen, *rec, DoOnlyAltA)
				&& rec->hasProp(METALIC_ATOM) ) {
				const double actual = distance2(nqoxygen.loc(), nearr->second);
				const double expected = nqoxygen.elem().covRad() + rec->elem().covRad();

				if ((actual >= (expected - 0.65))
					&& (actual <= (expected + 0.25)) ) {
					// we have determined that the Oxygen is "bonded" to a metal

					emptySet.push_front(std::make_pair(rec,rec->loc()));
					recordSkipInfo(FALSE, fixNotes, theHatom, nqoxygen, emptySet, "(metal ligand)");
					emptySet.pop_front();
					doNotAdjustSC = TRUE;
				}
			}
		}
	} // *** end of special case code for N/Q metal coord. ***

	if (pp.hasFeature(BONDBUMPFLAG) && (!doNotAdjustSC)) {
		const PDBrec& heavyAtom = a1;
		const double halfbondlen = heavyAtom.elem().covRad();
		const double  maxbondlen = halfbondlen
			+ ElementInfo::StdElemTbl().maxCovalentRadius();

		std::list< std::pair<PDBrec*, Point3d> > sym_nearr_list = xyz.neighbors( heavyAtom.loc(),
			halfbondlen + 0.1,
			maxbondlen  + 0.25);
		PDBrec* rec = NULL;
		for (std::list< std::pair<PDBrec*, Point3d> >::const_iterator nearr = sym_nearr_list.begin(); nearr != sym_nearr_list.end(); ++nearr) {
			rec = nearr->first;
			bool is_sym=(distanceSquared(rec->loc(), nearr->second) > 0.00001);  //Vishal: use TOL
			if (interactingConfs(heavyAtom, *rec, DoOnlyAltA)
				&& (! (sameres(theHatom, *rec, is_sym)
				&& StdResH::ResXtraInfo().match(theHatom.resname(), "AminoAcid") ))
				&& (! rec->elem().isHydrogen())
				&& (! rec->isWater()) ) {
				const double actual = distance2(heavyAtom.loc(), nearr->second);
				const double expected = heavyAtom.elem().covRad() + rec->elem().covRad();

				if (interactingConfs(theHatom, *rec, DoOnlyAltA)	// rec_is_not_sym not needed. Using as a place holder
					&& (actual >= (expected - 0.65))
					&& (actual <= (expected + 0.25)) ) {

					// this section determines how "bonded" groups place their hydrogens
					// especially as regards how metals determine the protonation pattern

					const double bias = rec->hasProp(METALIC_ATOM)
						? MetalBumpBias     // compensate for small "Metal VDW"
						: NonMetalBumpBias; // skip H very near non-Metal

					const double hRdist = distance2(theHpos, rec->loc());
					const double bumpLim = rec->elem().explRad() + bias;

					if (hRdist < bumpLim) {
						const char *msg = "(H bumps)";
						if ((bumpLim - hRdist) < 0.02) {
							cerr << "WARNING:" << heavyAtom.recName()
								<< " H atom too close to "
								<< rec->recName() << " by "
								<< (bumpLim - hRdist) << "A" << endl;
							cerr << "        you may need to adjust using "
								<< "-METALBump or -NONMETALBump." << endl;
						}
						if (pp.hasFeature(ROTATEFLAG)) { msg = "(short bond)"; }
						if ( rec->hasProp(METALIC_ATOM)
							&& strstr(":ASN:asn:GLN:gln:ASX:asx:GLX:glx:", theHatom.resname())
							&& (heavyAtom.elem().atno()==7)) { // N/Q Nitrogen metal ligand?
							cerr << "WARNING:" << heavyAtom.recName()
								<< " coordinates " << rec->recName()
								<< " (should consider flipping)" << endl;
						}
						emptySet.push_front(std::make_pair(rec,rec->loc()));
						recordSkipInfo(TRUE, fixNotes, theHatom, heavyAtom, emptySet, msg);
						emptySet.pop_front();
						return FALSE; // don't add hyd.
					}
				}
			}
		}
	}

	if (pp.hasFeature(IFNOPO4)) {
		const PDBrec& theoxygen = a1;
		double pobondlen = theoxygen.elem().covRad() +
			ElementInfo::StdElemTbl().element("P")->covRad();
		std::list< std::pair<PDBrec*, Point3d> > sym_nearr_list = xyz.neighbors( theoxygen.loc(),
			pobondlen - 0.25, pobondlen + 0.25);
		PDBrec* rec = NULL;
		for (std::list< std::pair<PDBrec*, Point3d> >::const_iterator nearr = sym_nearr_list.begin(); nearr != sym_nearr_list.end(); ++nearr) {
			rec = nearr->first;
			if (interactingConfs(theoxygen, *rec, DoOnlyAltA)
				&& rec->elem().atno() == 15 ) {
				return FALSE; // we have a OP bond -- don't add hyd.
			}
		}
	}
	return TRUE;
}

void recordSkipInfo(bool skipH, std::vector<std::string>& fixNotes,
					const PDBrec& theHatom, const PDBrec& heavyAtom,
					std::list< std::pair<PDBrec*, Point3d> >& nearr, const char * msg) {
	std::list<PDBrec*> conAtms;
	PDBrec* rec = NULL;
	for(std::list< std::pair<PDBrec*, Point3d> >::const_iterator it = nearr.begin(); it != nearr.end(); ++it) {
		rec = it->first;
		bool is_sym=(distanceSquared(rec->loc(), it->second) > 0.00001);  //Vishal: use TOL
		if (! (sameres(heavyAtom, *rec, is_sym)
			&& StdResH::ResXtraInfo().match(heavyAtom.resname(), "AminoAcid")) ) {
			conAtms.push_front(rec);
		}
	}
	std::string theOtherGuy = (!conAtms.empty()) ? ( (**(conAtms.begin())).recName()
		+ ( (conAtms.size() == 1) ? ":" : ":...") )
		: " cyclic :" ;

	std::string heavy = heavyAtom.atomname();

	if ((heavy == " O3'") || (heavy == " O3*") ||
		(heavy == " O5'") || (heavy == " O5*")) {
		for(std::list<PDBrec*>::const_iterator cat = conAtms.begin(); cat != conAtms.end(); ++cat) {
			std::string cat_elem = (*cat)->elemLabel();
			if (cat_elem == " P") {
				return; /* bypass -- do not warn about normal nucleic acid linkages */
			}
		}
	}

	fixNotes.push_back( (skipH ? "USER  MOD NoAdj-H:" : "USER  MOD NoAdj  :") +
		theHatom.recName() + ":" +
		heavyAtom.recName() + ":" +
		theOtherGuy + msg);

	if (Verbose) {
		if (skipH) {
			cerr << "SKIPPED H("<< theHatom.recName() <<"):"
				<< heavyAtom.recName() << " bonds";
		}
		else {
			cerr << "FIXED H("<< theHatom.recName() <<"):"
				<< heavyAtom.recName() << " bonds";
		}
		for(std::list<PDBrec*>::const_iterator it = conAtms.begin(); it != conAtms.end(); ++it) {
			cerr << "-" << (*it)->recName();
		}
		cerr << " " << msg << endl;
	}
}
